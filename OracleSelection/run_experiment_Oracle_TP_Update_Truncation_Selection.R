## Genomic selection simulations
#-------------------------------------------------------------------------
#Copyright (c) 2016 Jeff Neyhart
#-------------------------------------------------------------------------

load("Genome/CAP_Genome_QTL100_h0.5_100_Iterations.RData")#------------------------------Set Path To Genome

#source hypred map, imported library as tar.gz file
library(hypred)
library(rrBLUP)
library(GSSimTPUpdate)
library(rrBLUP)
library(parallel)
library(stringr)

save.dir<-"own_results"
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
n.cycles =15 #----------------------------------------------------------------------------------------Cycles

# Number of phenotyping environments and reps
n.env = 3
n.rep = 1

# Minor allele frequency cut-off for markers
min.maf = 0.03

# Barley population genetics data
mutation.rate.snp = 0
mutation.rate.qtl = 0

# Selection intensity and the number of crosses
parents.sel.intensity = 100 #pas aan -------------------------------------------------------------------------------
n.crosses = parents.sel.intensity/2

# The number of lines to add the TP after each cycle
tp.update.increment = 50
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
  experiment.sub.results <- mclapply(X = iters.per.core, FUN = function(iter.set) {
    
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
      G_G0 <- genotypic.value(genome= hv.genome, haploid.genos = TP.haploids.i) %>% 
        mean()
      G_G0_top <- sort(genotypic.value(genome= hv.genome, haploid.genos = TP.haploids.i),decreasing=TRUE)[1:10]%>%
        mean()
      
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
      
      id.tp=sample(length(TP.phenos),parents.sel.intensity)
      # Set dummy variables for the phenos and genos
      TP.phenos.i <- as.matrix(TP.phenos[id.tp,])
      TP.genos.i <- as.matrix(TP.genos[ id.tp,])
      
      # Create an initial data list
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
        qtl.out<-QTL.max(genome=hv.genome,  haploid = candidate.haploid.i)
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
        
        # For the whole genome, find the mean LD value across those
        ## max LD values per QTL
        
        #if the candidate.LD.genome is NA, next calculations cannot be calculated.
        if (!any(is.na(candidate.LD.genome))){
          
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
        } else {
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
                                            geno.train = TP.genos.i,                     
                                            geno.pred = candidate.marker.genos.i))    
              
        
        
        # Measure the phenotype and true genotypic values of all selection candidates
        candidate.values.i <- phenotype.population( genome = hv.genome,
                                                    haploid.genos = candidate.haploid.i,
                                                    V_E = V_E,
                                                    V_e = V_e,
                                                    n.env = n.env,
                                                    n.rep = n.rep )
                                      
        # If it errors out, just return NA for this iteration
        if (class(predictions.out) == "try-error")
          return(NA)
        
       
        
        ##### Update the TP --------------------------------------------------------------------
        #pre-allocate
        best.accuracu<-0
        Nochange<-F
        added<-0
        removed<-0
        
        
        
        name.add<-c()
        name.remove<-c()
        #check the prediction performance of each individual in the breeding population
        if (length(unique(candidate.values.i$geno.values))>1){
          
          for (tp.sel in seq(tp.update.increment)){
            ind.sel<-0
            for (ind in seq(parents.sel.intensity*ind.per.cross/2)){
              geno.temp<-candidate.marker.genos.i[ind,]
              pheno.temp<-candidate.values.i$mean.pheno.values[ind]
              
              predictions.temp <- try(make.predictions(pheno.train = c(TP.phenos.i, pheno.temp),
                                                       geno.train = rbind(TP.genos.i,geno.temp),                    
                                                       geno.pred = candidate.marker.genos.i))     
              
              names.pred<-rownames(candidate.values.i$geno.values)
              names.pred<-names.pred[!(names.pred%in%c(rownames(TP.phenos.i), names.pred[ind.sel]))]
              
              pred.accuracy.temp <- cor(candidate.values.i$geno.values[names.pred,],
                                     predictions.temp$GEBV[names.pred,]) %>%
                as.numeric()
              
              if (pred.accuracy.temp> best.accuracu){
                best.accuracu<-pred.accuracy.temp
                ind.sel=ind
              }
            }
            
            #If an individual increases the prediction performance, add the individual with the greatest contribution
            if (ind.sel!=0){
              TP.phenos.i <-rbind(TP.phenos.i, as.matrix(candidate.values.i$mean.pheno.values[ind.sel,]))
              TP.genos.i <-as.matrix(rbind(TP.genos.i, candidate.marker.genos.i[ind.sel,]))
              row.names(TP.genos.i)[length(TP.phenos.i)]<-row.names(as.matrix(candidate.values.i$mean.pheno.values[ind.sel,]))
              added<-added+1
              name.add<-c(name.add, row.names(as.matrix(candidate.values.i$mean.pheno.values[ind.sel,])))
              Nochange<-F
            }else{
              if (Nochange){
                break
              }else{
                Nochange=T
              }
            }
            
            
            #Check if removing an individual from the breeding population increases the prediction performance
            ind.sel=0
            
            for (ind in seq(length(TP.phenos.i))){
              
              
              predictions.temp <- try(make.predictions(pheno.train = c(TP.phenos.i[-ind]),
                                                       geno.train = rbind(TP.genos.i[-ind,]),                    
                                                       geno.pred = candidate.marker.genos.i))     
              
              names.pred<-rownames(candidate.values.i$geno.values)
              names.pred<-names.pred[!(names.pred%in%c(rownames(TP.phenos.i)[-ind]))]
            
              pred.accuracy.temp <- cor(candidate.values.i$geno.values[names.pred,],
                                      predictions.temp$GEBV[names.pred,]) %>%
              as.numeric()
              
              if (pred.accuracy.temp>best.accuracu){
                best.accuracu<-pred.accuracy.temp
                ind.sel=ind
              }
            }
            
            #If the removal of an individual increases the prediction performance, remove that individual
            if (ind.sel!=0){
              name.remove<-c(name.remove, row.names(as.matrix(TP.phenos.i[ind.sel,])))
              TP.phenos.i <-as.matrix(TP.phenos.i[-ind.sel,])
              TP.genos.i <-as.matrix(TP.genos.i[-ind.sel,])
              Nochange<-F
              removed<-removed+1
            }else{
              if (Nochange){
                break
              }else{
                Nochange=T
              }
            }
          }
        }
        # Validate the predictions
        # Find the correlation between the GEBVs and the true genotypic value
        notinTP<-!rownames(candidate.marker.genos.i)%in%rownames(TP.phenos.i)
        predictions.out <- try(make.predictions(pheno.train = TP.phenos.i,
                                                geno.train = TP.genos.i,                    
                                                geno.pred = candidate.marker.genos.i))  
        
        
        
        pred.accuracy.i <- cor(candidate.values.i$geno.values[notinTP],
                               predictions.out$GEBV[notinTP]) %>%
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
        
        
        
        
        M.var<-matrix(0,dim(TP.genos.i)[2],1)
        for(i in seq(dim(TP.genos.i)[2])){
          M.var[i]<-var(TP.genos.i[,i])
        }
        
        mean.phenos<-mean(TP.phenos.i)
        sd.phenos<-sd(TP.phenos.i)
        
        
        
        
        TP.info=list(meanPheno=mean.phenos, sdPheno= sd.phenos, MVAR=M.var, name.add=name.add, name.remove=name.remove)
        

                  
        
        print( paste("Cycle", breeding.cycle, "complete.") )
        
        cycle.name <- paste("cycle", breeding.cycle, sep = "")
        
        # Gather data for analysis
        simulation.results[[cycle.name]] <- 
          list(MM.solve = predictions.out$solve.out,
               candidate.values = candidate.values.i,
               QTL_info = qtl.out,
               QTL_Gen0=qtl.base,
               Genotype_Gen0=G_G0,
               Genotype_Gen0_top=G_G0_top,
               candidate.values=candidate.values.i,
               candidate.GEBV=predictions.out$GEBV,
               prediction.accuracy = pred.accuracy.i,
               parents = parent.lines.list,
               lengthTP= list(Length= length(TP.phenos.i), Added=added, Removed=removed),
               RelationshipTP=mu.relationship,
               Inbreeding=candidate.inbreeding,
               TP.info=TP.info
               )
      }          
      
      # Return a list
      list(sim.results = simulation.results, genome = hv.genome)
      
      # Close the sapply over replications
    }, simplify = F)
    
    # End mclapply
  }, mc.cores = n.cores)
  
  # Save the tp.change data
  filename <- file.path(save.dir, paste("sterwise_forward_backward_", pop.makeup, "_", change, "_", 
                                        tp.formation, "_", date, ".RData", sep = "") )
  
  
  save(list = c("experiment.sub.results", "change", "metadata"), file = filename)
  
  
  
} # Close the tp.change for loop


