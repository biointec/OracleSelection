## Genomic selection simulations
#-------------------------------------------------------------------------
#Copyright (c) 2016 Jeff Neyhart
#-------------------------------------------------------------------------

load("Genome/CAP_Genome_QTL100_h0.5_100_Iterations.RData")#------------------------------Set Path To Genome
source("R/OGM.R")


library(hypred) 
library(rrBLUP)
library(dplyr)
library(stringr)
library(Rcpp)
library(GSSimTPUpdate) 
library(parallel)
library(resample)
library(BGLR)
library(zoo)
library(GenomicMating)


  
save.dir<-"own_results"
# First and only argument is the pop makeup (MN, ND, or MNxND)
# Second argument is how the TP should be combined after each cycle (cumulative or window)
# Third argument is the number of QTL and the heritability
pop.makeup <- "MNxND" 
tp.formation <- "window"
h2 <- Genome$data$heritability #heritability
tp.change <- "tails"

# Verify that the TP changes in the arguments are acceptable
if (!all(tp.change %in% c("best", "worst", "random", "nochange", "PEVmean", "CDmean", "tails")) )
  stop("The TP change arguments are not acceptable.")

# Load the datasets
data("CAP.haploids")
data("CAP.markers")

scoping=3
# Set the number of cores by detection
max.cores <- detectCores()
n.cores <- 1#----------------------------------------------------------------------------NumberOfCores
#if number of selected cores exceed the available number of cores
if (n.cores>max.cores){
  n.cores<-max.cores
}

# Other simulation parameters
n.QTL <- Genome$data$number_QTL
# Make sure heritability is numeric
h2 <- as.numeric(h2)

# How many cycles?
n.cycles <-50 #----------------------------------------------------------------------------------------Cycles

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
n.crosses = parents.sel.intensity/2

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

# Iterate over the different tp.change types

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
  experiment.sub.results <- mclapply(X = iters.per.core, FUN = function(iter.set) {

    # Create a vector of rep names
    reps <- str_c("rep", iter.set)
    
    # Iterate over reps
    sapply(X = reps, FUN = function(r) {
      
    rep.iter<-as.numeric(substring(r,4))
    

    
    hv.genome <- Genome$genome[[rep.iter]]
      
      # All code below this line is variable in each iteration of the simulation

      #### Define trait parameters ####
      
      TP.haploids.i <- CAP.haploids
      # Convert the gametes to genotypes
      TP.genos <- genotype.loci(haploid.genos = TP.haploids.i, genome = hv.genome)
      
      #make a genetic map
      mask<-matrix(F,dim(CAP.haploids)[2],1)
      for(i in seq(dim(CAP.haploids)[2])){
        mask[i]<-any(CAP.markers$rs[i]==colnames(TP.genos))
      }
      
      map<-data.frame(CHROMOSOME=CAP.markers$chrom[mask],
                      POSITION=CAP.markers$pos[mask],
                      MARKER=as.character(CAP.markers$rs)[mask],
                      stringsAsFactors = FALSE)
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
      
      g1=TP.values
      
      
      G_G0<-genotypic.value(genome=hv.genome, haploid.genos =TP.haploids.i)  %>% 
        mean()
      G_G0_Top <-sort(genotypic.value(genome=hv.genome, haploid.genos=TP.haploids.i),decreasing=TRUE)[1:10]%>%
        mean()
      
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
      simulation.results <- list()
      
      # Loop over the number of cycles
      for (breeding.cycle in seq(n.cycles)) {

        ##### Start the Cycle Executions #####

        ##### Step 1 - Crossing and inbreeding
        # Make a crossing block
        if (breeding.cycle==1){
          
        crossing.block.i <- make.crossing.block(parent1.lines = parent.lines.list$p1, 
                                                parent2.lines = parent.lines.list$p2, 
                                                n.crosses = n.crosses, 
                                                use.parents.once = T)
        
        #source the function for next breeding cycles

        }else{
           
          
        # new crossing block construction
 
          crossing.block.i <- selection.pre$crossing.block.i

          
        }
        # According to the crossing block, parents are crossed to form F1s, then
        ## the progeny are inbred to the F3 generation. Since each F1 plant is inbred
        ## individually, the resulting families consist of F1:3 lines

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
        
        
        
        parents.fixed<-0
        for (i in seq(dim(candidate.marker.genos.i)[2])){
          if (all(candidate.marker.genos.i[,i]==-1) | all(candidate.marker.genos.i[,i]==1)){
            parents.fixed<-parents.fixed+1
          }
        }
        
        
        
        ##### Step 3 - Genotypic Summary Statistics
        # Measure the frequency of the 1 allele in each of the haploids set
        candidate.af.i <- measure.af(genome = hv.genome, 
                                     haploid.genos = candidate.haploid.i)
        
        TP.af.i <- measure.af(genome = hv.genome, 
                              haploid.genos = TP.haploids.i)
        
        source("R/QTL_maxv2.R")
        qtl.out<-QTL.max(genome=hv.genome,  haploid = candidate.haploid.i)
        
        if(breeding.cycle==1){
          source("R/QTL_maxv2_base.R")
          qtl_base<-QTL.max.base(genome=hv.genome, candidate.af.i=TP.af.i, haploid = TP.haploids.i)
        }
        ### LD measures
        # Candidates
        candidate.LD.window <- measure.LD(genome = hv.genome, 
                                          genos = candidate.genos.i, 
                                          Morgan.window = 0.25)
        
        # No window
        candidate.LD.genome <- measure.LD(genome = hv.genome,
                                          genos = candidate.genos.i)
        
        
        #if the candidate.LD.genome is NA, next calculations cannot be calculated.
        if (!any(is.na(candidate.LD.genome))){
        # For the whole genome, find the mean LD value across those
        ## max LD values per QTL
        candidate.mean.max.LD.genome <- apply(X = candidate.LD.genome, MARGIN = 1, FUN = function(qtl)
          max(qtl^2) ) %>%
          mean(na.rm = T)
        
        # Measure genomic LD on the TP
        TP.LD.genome <- measure.LD(genome = hv.genome, 
                                   genos = genotype.loci(haploid.genos = TP.haploids.i,
                                                         genome = hv.genome,
                                                         include.QTL = T) )
        
        # Measure the mean LD value across those
        ## max LD values per QTL in the TP
        TP.mean.max.LD.genome <- apply(X = TP.LD.genome, MARGIN = 1, FUN = function(qtl)
          max(qtl^2) ) %>%
          mean(na.rm = T)
        
        ## Persistance of LD phase
        # First find the common polymorphic QTL
        common.poly.QTL <- intersect( row.names(TP.LD.genome), row.names(candidate.LD.genome) )
        common.poly.markers <- intersect( colnames(TP.LD.genome), colnames(candidate.LD.genome) )
        
        # Subset the TP and candidates for those markers and QTL, then create a
        # data.frame
        TP.candidate.LD <- data.frame(TP = TP.LD.genome[common.poly.QTL, common.poly.markers] %>%
                                        as.vector(),
                                      candidates = candidate.LD.genome[common.poly.QTL, common.poly.markers] %>%
                                        as.vector())
        
        
        # Correlate
        TP.candidate.persistance.of.phase <- cor(TP.candidate.LD) %>% 
          .[upper.tri(.)]
        
        # Determine the number of loci used to calculate LD
        n.marker.LD <- min(ncol(TP.LD.genome), ncol(candidate.LD.genome))
        n.qtl.LD <- min(nrow(TP.LD.genome), nrow(candidate.LD.genome))
        
        
        # Create a list to save
        qtl.marker.LD.i <- list(sc.mean.max.genome = candidate.mean.max.LD.genome,
                                tp.mean.max.genome = TP.mean.max.LD.genome,
                                persistance.of.phase = TP.candidate.persistance.of.phase,
                                n.qtl.LD = n.qtl.LD,
                                n.marker.LD = n.marker.LD)
        }else{
          qtl.marker.LD.i<-NA #indicates that the qtl.marker.LD.i cannot be calculated.
        }
        
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
        
        if (breeding.cycle==1){
          Ne=NA
          Ft2=candidate.inbreeding
        }else{
          Ft1=Ft2
          Ft2=candidate.inbreeding
          delta_F=(Ft2-Ft1)/(1-Ft1)
          Ne=0.5/delta_F
        }
        
        ##### Step 4 - Prediction
        # Remove the markers with maf below the threshold (set at the start of the sim)
        ## also remove monomorphic markers

        
        # Estimate marker effects
        predictions.out <- try(make.predictions(pheno.train = TP.phenos.i,
                                            geno.train = TP.genos.i,
                                            geno.pred = candidate.marker.genos.i))
        
    
                            
        # If it errors out, just return NA for this iteration
        if (class(predictions.out) == "try-error")
          return(NA)
        
        ##### Step 5 - Phenotype the population
        # We use the haploid genotypes from the F1:3 generation to measure the
        ## the phenotypes
        
        # Measure the phenotype and true genotypic values of all selection candidates
        candidate.values.i <- phenotype.population( genome = hv.genome,
                                                    haploid.genos = candidate.haploid.i,
                                                    V_E = V_E,
                                                    V_e = V_e,
                                                    n.env = n.env,
                                                    n.rep = n.rep )
        
        # Validate the predictions
        # Find the correlation between the GEBVs and the true genotypic value
        pred.accuracy.i <- cor(candidate.values.i$geno.values,
                               predictions.out$GEBV) %>%
          as.numeric()


        
        ##### Step 6 - Select the parents of the next generation --------------------------------------- 6 Select Parents
        
        #H_UC_bridging.R
        selection.pre<-OGM_H(candidate.marker.genos.i,
                             candidate.haplotype=candidate.haploid.i,
                             Bhat=predictions.out$solve.out$u,
                             pheno.train,
                             p.sel=parents.sel.intensity,
                             presel=parents.sel.intensity*scoping,
                             A.sc,
                             parents.information=parent.selections.list,
                             trainingpanel=list(geno=TP.genos.i,pheno=TP.phenos.i),
                             map=map,
                             rep.iter,
                             genome=hv.genome)
        
        parent.selections.i<-selection.pre$parent.selections.i

        #remove the top-10 individuals
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
                              
        
        print( paste("Repeat", rep.iter,"Cycle", breeding.cycle, "complete.") )
        print( paste("GV of top-10:", round(mean(sort(candidate.values.i$geno.values,decreasing=T)[1:10])/42.8,2)))
        
        
        cycle.name <- paste("cycle", breeding.cycle, sep = "")
        
        # Gather data for analysis
        simulation.results[[cycle.name]] <- 
          list(geno.summary.stats = list(candidate.af = candidate.af.i,
                                         TP.af = TP.af.i,
                                         qtl.marker.LD = qtl.marker.LD.i,
                                         mu.TP.candidate.rel = mu.relationship),
               MM.solve = predictions.out$solve.out,
               Generatie_nul = g1,
               candidate.values = candidate.values.i,
               QTL_info = qtl.out,
               QTL_base=qtl_base,
               gen0=G_G0,
               gen0_top=G_G0_Top,
               CrossingBlock= crossing.block.i,
               candidate.values=candidate.values.i,
               candidate.GEBV=predictions.out$GEBV,
               parents = parent.lines.list,
               Ne=Ne,
               SR=scoping,
               tp.update = TP.addition.list,
               inbreeding = list(candidates = candidate.inbreeding,
                                 TP.additions = TP.addition.inbreeding),
               relationship = list(TP.candidates = mu.relationship,
                                   TP.additions = TP.addition.mu.relationship))
                                   
        
      } # Close the per-cycle loop
      
      # Return a list
      list(sim.results = simulation.results, genome = hv.genome)
      
      # Close the sapply over replications
    }, simplify = F)
  # End mclapply
  }, mc.cores = n.cores)
  
  # Save the tp.change data
  #filename <- file.path(save.dir, paste("simulation_results_", pop.makeup, "_", change, "_", 
  #                                      tp.formation,"_Scoping" ,scoping,"_", date, ".RData", sep = "") )
  filename <- file.path(save.dir, paste("OCS_v2_",n.iterations,"rep_", change,"_",date, ".RData", sep = "") )                      
  
  save(list = c("experiment.sub.results", "change", "metadata"), file = filename)
  
  

  } # Close the tp.change for loop


