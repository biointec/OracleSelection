## Calculate H-score
#-------------------------------------------------------------------------
#Copyright (c) 2019 Antoine Allier
#-------------------------------------------------------------------------

DeepScoping<-function(candidate.marker.genos.i,
                      candidate.haplotype,
                      GeneBank,
                      GeneBank.genotype,
         GeneBank.Phenotype,
         greedy,
         pre,
         sco,
         layers,
         predictions.out,
         parents.information=parent.selections.list,
         trainingpanel=list(geno=TP.genos.i,pheno=TP.phenos.i),
         Map=map,
         rep.iter,
         genome=hv.genome){

  Geno<-rbind(candidate.marker.genos.i, GeneBank.genotype)
 
  Bhat <- predictions.out$solve.out$u 
 
  GEBV <- Geno%*%Bhat
  
  # Assuming the Elites are the 10 best lines
  PopE <- names(GEBV[order(GEBV,decreasing = TRUE),])[seq(greedy)]
  Poptop<-names(GEBV[order(GEBV,decreasing = TRUE),])[seq(greedy*10)]
  
  # Donors are sampled among the remaining lines
  #PopD <- sample(setdiff(rownames(GEBV),PopE),50,replace = FALSE)
  PopD1 <- sample(rownames(GeneBank.genotype),length(rownames(GeneBank.genotype)),replace = FALSE) #blend in individuals of the genebank
  
  PopD2<-list()
  for (i in seq(layers)){
    if (i == 1){
      PopD2[[i]]<-names(GEBV[seq((greedy)*10+1,(greedy+pre+0*sco/4)*10),])
    }else{
      PopD2[[i]]<-names(GEBV[seq((greedy+pre+(i-2)*(sco/layers))*10+1,(greedy+pre+(i-1)*sco/layers)*10),])
    }
  }
  
  ###########
  ###########
  #Remove QTL markers from haplotype
  Haplo<-rbind(as.matrix(do.call("rbind", candidate.haplotype)),
               as.matrix(do.call("rbind",GeneBank)))
  
  positions<-c()
  begin=0
  for (i in seq(7)){
    positions<-c(positions, genome[[i]]@pos.add.qtl$ID+begin)
    begin<-begin+genome[[i]]@num.snp.chr+genome[[i]]@num.add.qtl.chr
  }
  
  Haplo<-Haplo[,-positions]
  #Get HEBVs
  ObjectHEBV <- GetHEBVmat(ObjectGeno = Haplo,
                           ObjectBeta = predictions.out$solve.out$u,
                           ObjectMap = Map,
                           WindowSize = WindowSizeUI,
                           StepSize = StepSizeUI)
  
  d1<-dim(ObjectHEBV$HEBV)
  HEBV<-c()
  for(i in seq(1,d1[1],by=2)){
    ind.hebv<-rbind(ObjectHEBV$HEBV[i,],ObjectHEBV$HEBV[i+1,])
    temp.hebv<-c()
    for (j in seq(d1[2])){
      temp.hebv<-c(temp.hebv,max(ind.hebv[,j]))
    }
    HEBV<-rbind(HEBV, temp.hebv)
  }
  
  line.names <- row.names(ObjectHEBV$HEBV) %>%
    str_replace(pattern = "\\.[0-9]$", replacement = "") %>%
    unique()
  
  rownames(HEBV)<-line.names
  colnames(HEBV)<- colnames(ObjectHEBV$HEBV)
  
  ObjectHEBV<-list(HEBV=HEBV,POSITION=ObjectHEBV$POSITION)
  ########
  ########
  #Calculate H-score for layer 0
  SelDonor1 <- c()
  SelDonor2 <- list()
  for (i in seq(layers+1)){
    SelDonor2[[i]]<-NULL
    if (i == layers+1){
      SelDonor2[[i]]<-0 #to create an empty list with NULL elements
    }
  }
  
  HSelDonor1 <- ComputeH(ObjectHEBVmat = ObjectHEBV$HEBV,# elites only
                         ListPopE = PopE,
                         ListPopD = NULL)
  
  HSelDonor2<-list()
  for(i in seq(layers)){
    HSelDonor2[[i]] <- ComputeH(ObjectHEBVmat = ObjectHEBV$HEBV,# elites only
                           ListPopE = PopE,
                           ListPopD = NULL)
  }
  
  invisible(lapply(1:(pre/2),function(iter){
    tmp <- ComputeH(ObjectHEBVmat = ObjectHEBV$HEBV,
                    ListPopE = c(SelDonor1,PopE),# PopE + selected donors
                    ListPopD = PopD1[!PopD1%in%c(SelDonor1)] # yet unselected donors
    )
    #Increment with selected donor
    SelDonor1 <<- c(SelDonor1,tmp[order(tmp$H,decreasing =TRUE),]$LINE[1])
    # Increment with H criterion value of elites + selected donors
    HSelDonor1 <<- rbind(HSelDonor1,tmp[order(tmp$H,decreasing =TRUE),][1,])
  }))
  
  #for each layer (1-4) calculate H-score
  for (i in seq(layers)){

  invisible(lapply(1:(sco/(layers*2)),function(iter){
    tmp <- ComputeH(ObjectHEBVmat = ObjectHEBV$HEBV,
                    ListPopE = c(SelDonor2[[i]],PopE),# PopE + selected donors
                    ListPopD = PopD2[[i]][!PopD2[[i]]%in%c(SelDonor2[[i]])] # yet unselected donors
    )
    #Increment with selected donor
    SelDonor2[[i]] <<- c(SelDonor2[[i]],tmp[order(tmp$H,decreasing =TRUE),]$LINE[1])
    # Increment with H criterion value of elites + selected donors
    HSelDonor2[[i]] <<- rbind(HSelDonor2[[i]],tmp[order(tmp$H,decreasing =TRUE),][1,])
  }))
  
  }

  # Assuming we have selected the donors:
  SelDonors4UC1 <- HSelDonor1$LINE[seq(2,pre/2+1)]
  SelDonors4UC2<-list()
  for (i in seq(layers)){
  SelDonors4UC2[[i]] <- HSelDonor2[[i]]$LINE[seq(2,sco/(2*layers)+1)]
  }

  
  
  PopTop<-list()
  for (i in seq(layers)){
    
    PopD2[[i]][1:(length(PopD2[[i]])-5)]
    pos<-c()
    for (j in seq(length( SelDonors4UC2[[i]]))){
      pos<-c(pos,which(PopD2[[i]]==SelDonors4UC2[[i]][j]))
    }
    PopTop[[i]] <- c(Poptop,PopD2[[i]][-pos])
  }
  
  #We have the Donor individual, now we will reselect the PopE population based on the greedy population with the eye on the preservation of the population
  #Build the crossing block for the layers 1-4
  cros.sco<-Crossing.Scoping(candidates_donor=SelDonors4UC2,
                   candidates_popE=PopTop,
                   poplayer=PopD2,
                   Z=Geno,
                   GEBV,
                   layers)
  
  #Build the crossing block for layer 0
  cros.pre<-Crossing.Scoping(candidates_donor=list(SelDonors4UC1),
                             candidates_popE=list(Poptop),
                             poplayer=PopD2,
                             Z=Geno,
                             GEBV,
                             layers=1)
  
  

  #Summarize the crossing information
  
  parent.sel<-parents.information$parent.selections.i
  cros.block<-parents.information$crossing.block.i
  
  block<-cbind(as.character(cros.pre[[1]]$p1), as.character(cros.pre[[1]]$p2))
  for(i in seq(layers)){
    block<-rbind(block, cbind(as.character(cros.sco[[i]]$p1), as.character(cros.sco[[i]]$p2)))
  }
  
  cros.block[seq((100-(pre+sco))/2+1,50),]<-block
  
  parent.sel$lines.sel<-c(parent.sel$lines.sel,block[,1],block[,2])
  
  
  
  parent.sel$value.sel<-rbind(parent.sel$value.sel, as.matrix(GeneBank.Phenotype$mean.pheno.values[cros.pre[[1]]$p1,]),as.matrix(GEBV[cros.pre[[1]]$p2,]))
  
  for(i in seq(layers)){
    parent.sel$value.sel<-rbind(parent.sel$value.sel, as.matrix(GEBV[cros.sco[[i]]$p1,]),as.matrix(GEBV[cros.sco[[i]]$p2,]))
    
  }
  
  return(list(parent.selections.i=parent.sel,crossing.block.i=cros.block))
  
}
