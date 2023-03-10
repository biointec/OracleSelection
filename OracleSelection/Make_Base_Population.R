## Genomic selection simulations
# This code runs iterations of long-term genomic selection simulations to 
## test the usefulness of updating a training population and how it impacts
## prediction accuracy, genetic variance, and selection potential
load("Genome/CAP_Genome_QTL100_h0.5_100_Iterations.RData")#------------------------------Set Path To Genome
cycle2save<-c(5,10,15,20)
#source hypred map, imported library as tar.gz file
library(hypred)
library("rrBLUP")
library(GSSimTPUpdate)
library(rrBLUP)
library(parallel)
library(stringr)
save.dir<-"Population"
# First and only argument is the pop makeup (MN, ND, or MNxND)
# Second argument is how the TP should be combined after each cycle (cumulative or window)
# Third argument is the number of QTL and the heritability
pop.makeup <- "MNxND" #wijst naar de datasets die gebruikt worden, MN en ND
tp.formation <- "window"
h2 <- Genome$data$heritability #heritability
tp.change <- "tails"

# Verify that the TP changes in the arguments are acceptable
if (!all(tp.change %in% c("best", "worst", "random", "nochange", "PEVmean", "CDmean", "tails")) )
  stop("The TP change arguments are not acceptable.")

# Load the datasets
data("CAP.haploids")
data("CAP.markers")

# Set the number of cores by detection

max.cores <- detectCores()
n.cores <- 10#----------------------------------------------------------------------------NumberOfCores
#if number of selected cores exceed the available number of cores
if (n.cores>max.cores){
  n.cores<-max.cores
}

# Other simulation parameters
n.QTL <- Genome$data$number_QTL
# Make sure heritability is numeric
h2 <- as.numeric(h2)

# How many cycles?
n.cycles =20 #----------------------------------------------------------------------------------------Cycles

# Number of phenotyping environments and reps
n.env = 3
n.rep = 1

# Minor allele frequency cut-off for markers
min.maf = 0.03

# Barley population genetics data
mutation.rate.snp = 0
mutation.rate.qtl = 0

# Selection intensity and the number of crosses
parents.sel.intensity = 100
n.crosses = 50

# The number of lines to add the TP after each cycle
tp.update.increment = 150
# Size of the TP to maintain - this is the same as the starting TP
tp.size <- nrow(CAP.haploids) / 2

# Parent selection and crossing parameters
ind.per.cross = 20
cycle.candidate.size = n.crosses * ind.per.cross

# Standardized selection intensity
std.sel.intensity = parents.sel.intensity / cycle.candidate.size

# Computation parameters
n.iterations = Genome$data$n.iterations

date <- format(Sys.time(), "%d%m%y-%H%M%S")

# Save the metadata to a list
metadata <- list(h2 = h2,
                 n.cycles = n.cycles,
                 n.QTL = n.QTL,
                 min.marker.maf = min.maf,
                 parents.sel.intensity = parents.sel.intensity,
                 n.env = n.env, 
                 n.rep = n.rep,
                 mutation.rate.qtl = mutation.rate.qtl,
                 mutation.rate.snp = mutation.rate.snp,
                 tp.update.increment = tp.update.increment,
                 tp.size = tp.size,
                 n.crosses = n.crosses,
                 ind.per.cross = ind.per.cross,
                 cycle.candidate.size = cycle.candidate.size,
                 std.sel.intensity = std.sel.intensity,
                 n.iterations = n.iterations,
                 pop.makeup = pop.makeup,
                 date = date)


#### Define genome characteristics ####
# Find the snps per chromsome   rs contains the name of markers, index the number of chromosome (factor), function  length
# tAPPLY will split cap.markers$rs into the different factors and give the length of each part
n.chr.snps <- Genome$data$n.chr.snps

# Find the chromsome lengths
#the cromosome length is calculated as the max distance we find a marker
chr.len <- Genome$data$chr.len

# Create a list of loci positions
genetic.map.list <- Genome$data$genetic.map.list

for (change in tp.change) {
  set.seed(1)
  # Split iterations into cores
  if (n.cores > 1) {
    iters.per.core <- split(x = 1:n.iterations, factor(cut(x = 1:n.iterations, breaks = n.cores)))
    names(iters.per.core) <- paste("set", seq_along(iters.per.core), sep = "")
  } else {
    iters.per.core <- 1:n.iterations
  }    
  
  # Apply the iterations over cores
  Simulated.populations <- mclapply(X = iters.per.core, FUN = function(iter.set) {
    
    # Create a vector of rep names
    reps <- str_c("rep", iter.set)
    
    # Iterate over reps
    sapply(X = reps, FUN = function(r) {
      
      rep.iter<-as.numeric(substring(r,4))
      # All code below this line is variable in each iteration of the simulation

      #### Define trait parameters ####
      hv.genome <- Genome$genome[[rep.iter]]
      
      TP.haploids.i <- CAP.haploids
      # Convert the gametes to genotypes
      TP.genos <- genotype.loci(haploid.genos = TP.haploids.i, genome = hv.genome)
      
      # Find the allele freq of all marker snps
      marker.af <- measure.af(genome = hv.genome, haploid.genos = TP.haploids.i)$snp
      
      loci.af<-measure.af(genome = hv.genome, haploid.genos = TP.haploids.i)$loci
      Fixed<-sum( loci.af==0 | loci.af == 1)
      # Calculate maf minor allele freq
      marker.maf <- sapply(marker.af, FUN = function(freq) min(freq, 1 - freq))
      
      # Determine which are below 
      markers.below.maf <- which(marker.maf < min.maf)
      
      
      # Calculate the frequency of the 1 allele in the base training population
      p_i <- apply(X = TP.genos + 1, MARGIN = 2, FUN = mean) / 2
      # Calculate the P matrix
      P = matrix(2 * (p_i - 0.5))
      # Calculate the normalization constant
      c = 2 * sum(p_i * (1 - p_i))
      
      ## Select the TP lines for use as parents
      # First separate MN and ND lines
      line.names <- row.names(TP.genos)
      ND.lines <- str_subset(string = line.names, pattern = "^ND")
      MN.lines <- setdiff(x = line.names, y = ND.lines)
      
      ## Set the inital variances for the heritability
      # True genetic variance
      TP.V_g <- genotypic.value(genome = hv.genome, haploid.genos = TP.haploids.i) %>%
        var()
      # Environmental variance (scale * 8 as in Bernardo 2015)
      V_E = TP.V_g * 8
      # Residual variance scaled to achieve desired h2
      V_e = n.rep * n.env * ((TP.V_g / h2) - TP.V_g)
      
      
      # Phenotype the training population
      
      TP.values <- phenotype.population(genome = hv.genome,
                                        haploid.genos = TP.haploids.i,
                                        V_E = V_E,
                                          V_e = V_e,
                                        n.env = n.env,
                                        n.rep = n.rep,
                                        run.anova = T)
      
      
      
      TP.phenos <- TP.values$mean.pheno.values
      
      # Next select the top MN and top ND
      top.MN.lines <- names(sort(TP.phenos[MN.lines,], decreasing = T)[1:parents.sel.intensity])
      top.ND.lines <- names(sort(TP.phenos[ND.lines,], decreasing = T)[1:parents.sel.intensity])
      
      # Create a list to store the line names
      pop.makeup.list <- list(
        MN = list(p1 = top.MN.lines, p2 = top.MN.lines),
        ND = list(p1 = top.ND.lines, p2 = top.ND.lines),
        MNxND = list(p1 = top.MN.lines[seq(parents.sel.intensity / 2)], p2 = top.ND.lines[seq(parents.sel.intensity / 2)])
      )
      
      # Set the parent gamete data input
      parent.lines.list <- pop.makeup.list[[pop.makeup]]
      parent.haploids <- select.haploids(haploid.genos = TP.haploids.i, line.names = unique(unlist(parent.lines.list)))
      
      
      # Set dummy variables for the phenos and genos
      TP.phenos.i <- TP.phenos
      TP.genos.i <- TP.genos
      
      # Create an initial data list
      simulation.population <- list()
      simulation.results <- list()
      # Loop over the number of cycles
      for (breeding.cycle in seq(n.cycles)) {

        ##### Start the Cycle Executions #####  

        ##### Step 1 - Crossing and inbreeding
        # Make a crossing block
       
        crossing.block.i <- make.crossing.block(parent1.lines = parent.lines.list$p1, 
                                                parent2.lines = parent.lines.list$p2, 
                                                n.crosses = n.crosses, 
                                                use.parents.once = T)
        
  
        # According to the crossing block, parents are crossed to form F1s, then
        ## the progeny are inbred to the F3 generation. Since each F1 plant is inbred
        ## individually, the resulting families consist of F1:3 lines
        library(dplyr)
        library(stringr)
        library(Rcpp)
     
        
        candidate.haploid.i <- make.population(genome = hv.genome, 
                                               parental.haploids = parent.haploids,
                                               crossing.block = crossing.block.i,
                                               N = ind.per.cross,
                                               cycle.number = breeding.cycle,
                                               generations = 2,
                                               pop.type = "inbred",
                                               mutation.rate.snp = mutation.rate.snp,
                                               mutation.rate.qtl = mutation.rate.qtl)
        

        
        #--------------------------------------------------------------------------------------------------------------
        ##### Step 2 - Genotype
        # Find the genotypes of the markers and QTL
        candidate.genos.i <- genotype.loci(haploid.genos = candidate.haploid.i, 
                                           genome = hv.genome, 
                                           include.QTL = T)
        # Just the marker genotypes
        candidate.marker.genos.i <- genotype.loci(haploid.genos = candidate.haploid.i, 
                                                  genome = hv.genome, 
                                                  include.QTL = F)
        
        
        ##### Step 3 - Genotypic Summary Statistics
        # Measure the frequency of the 1 allele in each of the haploids set
        candidate.af.i <- measure.af(genome = hv.genome, 
                                     haploid.genos = candidate.haploid.i)
        
        TP.af.i <- measure.af(genome = hv.genome, 
                              haploid.genos = TP.haploids.i)
        
        source("R/QTL_maxv2.R")
        qtl.out<-QTL.max(genome=hv.genome, haploid = candidate.haploid.i)
        if (breeding.cycle==1){
        source("R/QTL_maxv2_base.R")
        qtl.base<-QTL.max.base(genome=hv.genome, candidate.af.i = TP.af.i,haploid= TP.haploids.i)
        }
        ### LD measures
        # Candidates
        
        candidate.LD.window <- measure.LD(genome = hv.genome, 
                                          genos = candidate.genos.i, 
                                          Morgan.window = 0.25)
        
        # No window
        candidate.LD.genome <- measure.LD(genome = hv.genome,
                                          genos = candidate.genos.i)
        
        ### Measure the average relationship between the TP and the candidates
        # Assign M
        M <- rbind(TP.genos.i, candidate.marker.genos.i)
        # Subtract P to make Z (need to convert P into a repeated matrix)
        W = M - matrix(P, nrow(M), length(P), byrow = T)
        # Calculate the relationship matrix
        A = tcrossprod(W) / c
        
        # Subset the relationship matrix for the selection candidates
        A.sc <- A[row.names(candidate.marker.genos.i), row.names(candidate.marker.genos.i)]
        
        # Calculate the mean relationship between the TP and the candidates
        mu.relationship <- A[row.names(TP.genos.i), row.names(candidate.marker.genos.i)] %>%
          mean()
        
        ## Find the average inbreeding coefficient among the selection candidates
        candidate.inbreeding <- A.sc %>%
          diag() %>%
          - 1 %>%
          mean()
        
        
        ##### Step 4 - Prediction
        # Remove the markers with maf below the threshold (set at the start of the sim)
        ## also remove monomorphic markers

        
        markers.to.remove <- c( which(candidate.af.i$snp == 0), 
                                which(TP.af.i$snp == 0),
                                markers.below.maf ) %>%
          unique() %>%
          sort()

        
        # Filter the TP and candidate marker matrices for those markers
        TP.genos.use <- TP.genos.i[,-markers.to.remove]
        candidate.genos.use <- candidate.marker.genos.i[,-markers.to.remove]
        
        # Estimate marker effects
        predictions.out <- try(make.predictions(pheno.train = TP.phenos.i,
                                            geno.train = TP.genos.use,
                                            geno.pred = candidate.genos.use))
                                            
        # If it errors out, just return NA for this iteration
        if (class(predictions.out) == "try-error")
          return(NA)
        
        ##### Step 5 - Phenotype the population
        # We use the haploid genotypes from the F1:3 generation to measure the
        ## the phenotypes
        
        # Measure the phenotype and true genotypic values of all selection candidates
        candidate.values.i <- phenotype.population( genome = hv.genome,
                                                    haploid.genos = candidate.haploid.i,
                                                    V_E = V_E[1],
                                                    V_e = V_e[1],
                                                    n.env = n.env,
                                                    n.rep = n.rep )
        
        # Validate the predictions
        # Find the correlation between the GEBVs and the true genotypic value
        pred.accuracy.i <- cor(candidate.values.i$geno.values,
                               predictions.out$GEBV) %>%
          as.numeric()


        
        ##### Step 6 - Select the parents of the next generation --------------------------------------- 6 Select Parents
        
        # Make selections on the GEBVs
        # Select the top 100 based on GEBVs for parents of the next cycle
      #  source("select_population_SA.R")
      #  parent.selections.i<-select.populationSA(value.mat=predictions.out$GEBV,scoping) 
       
       
        parent.selections.i <- select.population(value.mat = predictions.out$GEBV, 
                                                sel.intensity = parents.sel.intensity, 
                                                selection = "best") 
                                   
        
        parent.lines.list <- list(p1 = parent.selections.i$lines.sel, 
                                  p2 = parent.selections.i$lines.sel)
        
        # The parents are selected and crossed at the F3 stage, so subset the haploid genotpyes from the F1:3
        parent.haploids <- select.haploids(haploid.genos = candidate.haploid.i,
                                          line.names = parent.selections.i$lines.sel)
        
        parent.values <- select.values(pheno.values.list = candidate.values.i, 
                                       line.names = parent.selections.i$lines.sel)
        
        ##### Step 7 - Update the TP --------------------------------------------------------------------7 Update TP
        
        # Skip this step if not called
        if (change != "nochange") {
          
          # If the TP change is best, worst, or random, simply subset the population.
          if (change %in% c("best", "worst", "random", "tails")) {
            
            TP.addition.list <- 
              list(TP.addition.lines = select.population(value.mat = predictions.out$GEBV,
                                                         sel.intensity = tp.update.increment,
                                                         selection = change)$lines.sel )
          }

          if (change %in% c("PEVmean", "CDmean")) {
            
            # Analyze using PEVmean or CDmean
            # We want to see what optimized TP is best for the parents, so we will optimize the training set based
            ## on the lines from the whole candidate set, including the parents
            phenotyped.index <- which( row.names(candidate.marker.genos.i) %in%
                                         row.names(A.sc) )
            unphenotyped.index <- which( parent.selections.i$lines.sel %in%
                                           row.names(A.sc) )

            # V_e is estimated from maximum likelihood
            V_e.i <- predictions.out$solve.out$Ve
            # V_a is estimated as the variance among marker effects * the number of markers
            V_a.i <- predictions.out$solve.out$Vu * ncol(TP.genos.use)
            
            # Run the optimization algorithm
            optimized.TP <- optimize.subset(method = change, A = A.sc,
                                            n.training = tp.update.increment,
                                            phenotyped.index = phenotyped.index,
                                            unphenotyped.index = unphenotyped.index,
                                            V_a = V_a.i, V_e = V_e.i, max.iter = 500)
            
            # Using the optimized index of "phenotyped" entries, determine
            # which candidates should be added to the TP
            TP.addition.list <- list(TP.addition.lines = row.names(A.sc)[optimized.TP$OTP],
                                     optimization.info = optimized.TP)
            
          } # Close the tp optimization algorithm if statement
          
          ### Collect information on the new TP lines
          
          # TP additions
          TP.addition.lines <- TP.addition.list$TP.addition.lines
          
          TP.addition.haploids <- select.haploids(haploid.genos = candidate.haploid.i,
                                                  line.names = TP.addition.lines)
          
          # Subset the geno matrix for these lines
          TP.addition.genos <- candidate.marker.genos.i[TP.addition.lines,]

          # Gather genotypic and phenotypic values of the TP additions
          TP.addition.values <- select.values(pheno.values.list = candidate.values.i, 
                                              line.names = TP.addition.lines)
          # Separate the phenotypes
          TP.addition.phenos <- TP.addition.values$mean.pheno.values
          
          # Measure the average relationship among these lines
          TP.addition.mu.relationship <- A[TP.addition.lines, TP.addition.lines] %>%
            .[upper.tri(.)] %>% 
            mean()
          
          # Measure inbreeding among the additions
          TP.addition.inbreeding <- A[TP.addition.lines, TP.addition.lines] %>%
            diag() %>% 
            - 1 %>%
            mean()
          
          # Combine the new data to the TP
          TP.phenos.i <- rbind(TP.phenos.i, TP.addition.phenos)
          TP.genos.i <- rbind(TP.genos.i, TP.addition.genos)
          TP.haploids.i <- rbind(TP.haploids.i, TP.addition.haploids)
          
          ## Measure the expected heterozygosity of the additions, using all
          ## loci including QTL
          TP.addition.list[["Exp.het"]] <- 
            measure.expected.het(genome = hv.genome, 
                                 haploid.genos = TP.addition.haploids)
          
          
          # If the TP formation calls for a sliding window, use only the ~750 most recent training individuals
          if (tp.formation == "window") {
            
            # If the breeding cycle * tp addition size is less than the starting tp size, 
            ## include the most recent 150 additions, then randomly sample from the
            ## remaining TP
            if ((breeding.cycle * tp.update.increment) < tp.size) {
              
              # Find the index of the most recent additions
              tp.recent.index <- tail(1:nrow(TP.genos.i), (breeding.cycle * tp.update.increment))
              # Randomly select among the remaining index to maintain the tp.size
              tp.random.index <- sort( sample(setdiff(1:nrow(TP.genos.i), tp.recent.index), size = (tp.size - length(tp.recent.index)) ) )
              # Combine
              tp.keep.index <- sort(c(tp.recent.index, tp.random.index))
              
            } else { # Otherwise just take the last tp.size individuals added to the TP
              tp.keep.index <- tail(1:nrow(TP.genos.i), tp.size)
              
            }
              # Set the TP.pheno and TP.genos
              TP.phenos.i <- as.matrix(TP.phenos.i[tp.keep.index,])
              TP.genos.i <- as.matrix(TP.genos.i[tp.keep.index,])
              TP.haploids.i <- select.haploids(haploid.genos = TP.haploids.i,
                                               line.names = row.names(TP.genos.i))
          }
          
        } else {
          
          TP.addition.list <- list(TP.addition.lines = NA,
                                   Exp.het = NA)
          
          TP.addition.inbreeding <- NA
          TP.addition.mu.relationship <- NA
          
          
        } # Close the tp.change if statement
                              
        
        print( paste("Cycle", breeding.cycle, "complete.") )
        
        cycle.name <- paste("cycle", breeding.cycle, sep = "")
        simulation.results[[cycle.name]] <- 
          list(MM.solve = predictions.out$solve.out,
               candidate.values = candidate.values.i,
               QTL_info = qtl.out,
               candidate.values = candidate.values.i,
               candidate.GEBV = predictions.out$GEBV,
               prediction.accuracy = pred.accuracy.i,
               parents = crossing.block.i,
               tp.update = TP.addition.list)
        

        #save populationat certain breeding cycles
        if (any(breeding.cycle==cycle2save)){
        
        cycle.name <- paste("cycle", breeding.cycle, sep = "")
        
        # Gather data for analysis
        simulation.population[[cycle.name]] <- list(
                Haploids=candidate.haploid.i,
                Genotype=candidate.marker.genos.i,
                Phenotype=candidate.values.i,
                QTL.info=qtl.out,
                TP=list(pheno.train = TP.phenos.i,
                        geno.train = TP.genos.i,
                        haploids.train=TP.haploids.i)
                 )
        }
        
      }          

      # Return a list
      list(population = simulation.population, simulation=simulation.results, genome = hv.genome)
      
      # Close the sapply over replications
    }, simplify = F)
    
    # End mclapply
  }, mc.cores = n.cores)
  
  # Save the tp.change data
  filename <- file.path(save.dir, paste("Population_", pop.makeup, "_", change, "_", 
                                        tp.formation, "_", date, ".RData", sep = "") )
  
  
  save(list = c("Simulated.populations", "change", "metadata"), file = filename)
  
  
  
} # Close the tp.change for loop


