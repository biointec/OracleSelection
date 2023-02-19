built.bank<-function(genome    = hv.genome,
                     Haploids  = CAP.haploids,
                     Genotype  = TP.genos,
                     Phenotype = TP.phenos,
                     size      = 500){
    
    n.ind<-dim(Genotype)[1]
    n.parents<-size/ind.per.cross
    #make an empty crossing block (pre-allocatie)
    crossing.block.i <- matrix("",n.parents,2)
    
    line.names <- row.names(Genotype)
    ND.lines <- str_subset(string = line.names, pattern = "^ND")
    MN.lines <- setdiff(x = line.names, y = ND.lines)
    
    #split parent into two populations MN and ND lines
    MN.lines <- names(sort(Phenotype[MN.lines,], decreasing = T))
    ND.lines <- names(sort(Phenotype[ND.lines,], decreasing = T))
    
   # merkers<-matrix(T,dim(Genotype)[2],1)
   # allele1<-matrix(F,dim(Genotype)[2],1)
   # allele2<-matrix(F,dim(Genotype)[2],1)
    
    available<-seq(length(ND.lines))
    chosen<-c()
    
    #For each individual select the one that maximizes the S-Score (max variation)
    for (i in seq(n.parents)){
      
      p1<-MN.lines[((i-1)%%(length( MN.lines)))+1]
      score<-0
      for(j in available){
        score.temp<- sum(colVars(rbind(Genotype[p1,],Genotype[ND.lines[j],])))
        if (score.temp>score){
          score<-score.temp
          chosen<-j
        }
      }
      
      
      available<-available[-which(available==chosen)]
      
      if (length(available)==0){
        available<-seq(length(ND.lines))
      }
      
      
      p2<-ND.lines[chosen]
        
      crossing.block.i[i,1]<-p1
      crossing.block.i[i,2]<-p2
     
    }
    
  
    if( length(unique(c(crossing.block.i[,1], crossing.block.i[,2])))!=n.parents*2 ){
      stop("Non-unique parents have been selected")
    }
    
    #construct haploids based on the selected parents
    candidate.haploid.i <- make.population(genome = genome, 
                                           parental.haploids = Haploids,
                                           crossing.block = crossing.block.i,
                                           N = ind.per.cross,
                                           cycle.number = 0,
                                           generations = 2,
                                           pop.type = "inbred",
                                           mutation.rate.snp = mutation.rate.snp,
                                           mutation.rate.qtl = mutation.rate.qtl)
  

  return(candidate.haploid.i)
}
