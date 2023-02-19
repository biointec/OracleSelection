#use a separate genebank for prebreeding

Prebreeding<-function(GEBV = predictions.out$GEBV,
                      Z = candidate.marker.genos.i,
                      parents.sel=parents.sel.intensity,
                      fraction_greedy = 50,
                      fraction_pre=20,
                      fraction_greedyII=20,
                      GeneBank_G= GeneBank.genotype,
                      GeneBank_P= GeneBank.Phenotype,
                      parents.information=parent.selections.list,
                      G,
                      UseBank=T){
  #get size
  n.pop<-dim(Z)[1]
  k.pop<-dim(Z)[2]
  
  #get names and values
  ind.names <- rownames(GEBV)
  ind.values<- GEBV
  bank.names<- rownames(GeneBank_G)
  
  #-------------------------------------------------------------------------------------------------------------------------------
  #The Greedy Selection, importing information
  #-------------------------------------------------------------------------------------------------------------------------------
  
  #Information of parents selected using a greedy approach
  id.top<-order(GEBV,decreasing =T)[seq(1,(fraction_greedy))]
  Z.top<-Z[id.top,]
  
  #Fixed markers in the top/breeding population
  Markers.topfixed<-seq(k.pop)[colVars(Z.top)==0] #(resample library required) Fixed markers in greedy selected population
  Markers.fixed<-seq(k.pop)[colVars(Z)==0] #fixed markers in breeding population
  
  #which markers are of interest to reintroduce in the population
  marker.mask<-matrix(FALSE,k.pop,1)
  marker.mask[Markers.topfixed]<-TRUE
  marker.mask[Markers.fixed]<-FALSE #is T if marker is fixed in the greedy population but not in the breeding population
  #Markers.interest are markers that are fixed in the top population but not in the breeding population
  Markers.interest<-which(marker.mask)
  #genotype of fixed markers in the greedy selected population
  Z.topfixed<-Z.top[1,Markers.interest] 
  


  
  #-------------------------------------------------------------------------------------------------------------------------------
  #Preselection
  #-------------------------------------------------------------------------------------------------------------------------------
  
  #preallocation
  selected<-c()
  
  
  #check thatfraction_pre is even
  if (fraction_pre %% 2 != 0){
    stop("fraction_pre should be an even number")
  }
  
  if(UseBank){
  
    n.bank<-dim(GeneBank_G)[1]
    id.available<-seq(n.bank)
    bolean.available<-matrix(TRUE,n.bank,1)
  
  
    for (round in seq(fraction_pre/2)){
      ind.diff<-matrix(0,n.pop,1)
      for (i in seq(n.bank)){
        ind.diff[i]<-sum(GeneBank_G[i,Markers.interest]==(Z.topfixed*-1))#/(length(Z.topfixed))
      }

      p1.available<-which.max(ind.diff[bolean.available])
      p1<-id.available[bolean.available][p1.available]
      p1.geno<-GeneBank_G[p1,]
      #only p1 cannot be selected twice as p1
      bolean.available[p1]<-FALSE
    
      #compare genotype with top candidates
      diff.markers<-abs((Z.top[(round-1)%%fraction_greedy+1,]-p1.geno))==2 #search were there are 2 alleles different with the top parents
      diff.genotype<-p1.geno[diff.markers]*(-1) #find an individual that has the complementary allelles to reproduce an ideal offspring
   
      ind.diff<-matrix(0,n.pop,1)
      for (i in seq(n.pop)){
        ind.diff[i]<-sum(Z[i,diff.markers]==(diff.genotype))/(length(diff.genotype))
      }
      p2<-which.max(ind.diff)
      p2.geno<-Z[p2,]
      selected<-rbind(selected,cbind(bank.names[p1],ind.names[p2]))
    
      #redefine the markers of interest
      Markers.selected<-p1.geno[Markers.interest]!=(Z.topfixed) | p2.geno[Markers.interest]!=(Z.topfixed)
      Markers.interest<-Markers.interest[!Markers.selected]
      #if all markers have been selected, start over
      if (length(Markers.interest)==0){
        Markers.interest<-which(marker.mask) #if all fixed markers are selected, restart, the available vector prevents the selection of the same p1 parents
      }
    Z.topfixed<-Z.top[1,Markers.interest] 
    }
  }
  
  
  #---------------------------------------------------------------------------------------------------------------------------------
  #The Greedy selection II
  #---------------------------------------------------------------------------------------------------------------------------------

  #select the candidates 
  
  if(UseBank){
    top.id<-order(GEBV[seq(n.pop-(fraction_greedyII+fraction_pre)*10+1,n.pop-fraction_greedyII*10)], decreasing=T)[1:(fraction_greedyII/2)]
    parents<-ind.names[seq(n.pop-(fraction_greedyII+fraction_pre)*10+1,n.pop)][top.id] 

    top.id<-order(GEBV[seq(n.pop-(fraction_greedyII)*10+1,n.pop)], decreasing=T)[1:(fraction_greedyII/2)]
    parents<-c(parents, ind.names[seq(n.pop-(fraction_greedyII)*10+1,n.pop)][top.id] )
    topper.id<-order(GEBV[parents,],decreasing=T)
    parents<-parents[topper.id]
  }else{
    top.id<-order(GEBV[seq(n.pop-(fraction_greedyII+fraction_pre)*10+1,n.pop-fraction_greedyII*10)], decreasing=T)[1:((fraction_greedyII+fraction_pre)/2)]
    parents<-ind.names[seq(n.pop-(fraction_greedyII+fraction_pre)*10+1,n.pop)][top.id] 
    
    top.id<-order(GEBV[seq(n.pop-(fraction_greedyII)*10+1,n.pop)], decreasing=T)[1:((fraction_greedyII+fraction_pre)/2)]
    parents<-c(parents, ind.names[seq(n.pop-(fraction_greedyII+fraction_pre)*10+1,n.pop)][top.id] )
    topper.id<-order(GEBV[parents,],decreasing=T)
    parents<-parents[topper.id]
    
  }

  
  #build a crossing block
  p1<-c()
  p2<-c()
  available<-seq(length(topper.id))
  
  for (i in seq(length(topper.id)/2)){
    p1<-c(p1,available[1])
    available<-available[-1]
    
    G_pop<-G[parents[p1[i]],parents[available]]
    
    id.max<-which.min(G_pop)
    p2<-c(p2,available[id.max])
    available<-available[-id.max]
  }
  
  p1<-parents[p1]
  p2<-parents[p2]
  
  parent.sel<-parents.information$parent.selections.i
  cros.block<-parents.information$crossing.block.i
  
  #add the selected to the already selected parents
  Parents<-c(selected[,1],selected[,2],p1,p2)
  parent.sel$lines.sel<-c(parent.sel$lines.sel,Parents)
  
  if (UseBank){
    parent.sel$value.sel<-as.matrix(c(parent.sel$value.sel, as.matrix(GeneBank_P$mean.pheno.values[Parents[seq(1,fraction_pre/2)],]), as.matrix(GEBV[Parents[seq((fraction_pre/2)+1,fraction_pre+fraction_greedyII)],])))
  }else{
    parent.sel$value.sel<-as.matrix(c(parent.sel$value.sel, as.matrix(GEBV[Parents,])))
  }
  rownames(parent.sel$value.sel)<-parent.sel$lines.sel
  
    cros.block[seq((parents.sel-(fraction_pre+fraction_greedyII))/2+1,50),1]<-c(selected[,1],p1)
    cros.block[seq((parents.sel-(fraction_pre+fraction_greedyII))/2+1,50),2]<-c(selected[,2],p2)
    
  #}else{
    
   # parent.sel$value.sel<-c(parent.sel$value.sel, as.matrix(GEBV[Parents,]))
  #  cros.block[seq((parents.sel-(fraction_greedyII))/2+1,50),1]<-c(selected[,1],p1)
  #  cros.block[seq((parents.sel-(fraction_greedyII))/2+1,50),2]<-c(selected[,2],p2)
 # }

  return(list(parent.selections.i=parent.sel,crossing.block.i=cros.block))
}