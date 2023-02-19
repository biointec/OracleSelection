#' Rapid breeding selection
#' @value.mat contains the gebv of all the individuals
#' @sel.intensity The number of individuals that have te be selected as parent
#' @genos The genotype of all the individuals
#' @fraction_greedy The fraction that need to be selected greedily
#' @G the relationship matrix calculated according to VanRaden
#' 
#' This method selects parents and return a list with the selected parents and crossing block
Elite.Selection<-function(value.mat = predictions.out$GEBV, 
                                   sel.intensity = parents.sel.intensity,
                                   genos = candidate.marker.genos.i,
                                   fraction_greedy,
                                   fraction_pre,
                                   G){
  if(parents.sel.intensity%%2!=0){
    stop("parent.sel.intensity should be an even number, parents are always selected as couple")
  }
  G_threshold<-1
  
  #dimention of population
  k_pop<-dim(genos)[2]
  
  #keep track of the selected alleles
  allele1<-matrix(0,k_pop,1)
  allele2<-matrix(0,k_pop,1)
  p_select<-matrix(TRUE,k_pop,1)
  
  #standirize the gebv
  gebv<-as.matrix(((value.mat-mean(value.mat))/sd(value.mat))[c(seq(1,fraction_greedy*10),seq((fraction_greedy+fraction_pre)*10+1,1000)),])
  ind.names<-rownames(gebv)
  genotype<-genos[ind.names,]
  n_pop<-length(gebv)[1]
  #Keep track of the population
  id_pop<-seq(n_pop)
  available_parents<-matrix(TRUE,n_pop,1)
  available_id<-id_pop[available_parents]
  
  
  #----------------------------------------------------------------------------------------------
  #select the top individuals  
  #----------------------------------------------------------------------------------------------  
  
  #Calculate the number of parents that need to be selected greedily
  #Parents are always selected per couple, so the number to select should be whole even numbers
  number_top<-round((fraction_greedy/100)*sel.intensity)
  if (number_top%%2 != 0){
    number_top = number_top+1
  }
  list.top<-order(gebv,decreasing=TRUE)
  parents.top<-c()
  top<-1
  
  if (number_top>0){
    for (i in seq(to=number_top, by=2)){
      while (!available_parents[list.top[top]]){
        top<-top+1 #If the next top parent has already been selected, go to the next top parent
      }
      #select the first parent
      P1<-list.top[top]
      parents.top<-c(parents.top, P1)
      #change status of the selected parent
      available_parents[P1]<-FALSE 
      #update top parent
      top<-top+1
      
      
      top.candidates<-list.top[seq(number_top)]
      top.av.candidates<-top.candidates[available_parents[top.candidates]]
      P2<-top.av.candidates[which.min(G[P1,top.av.candidates])]
      
      parents.top<-c(parents.top, P2)
      #change status of the selected parent
      available_parents[P2]<-FALSE 
      
      #update the selected alleles
      allele1[genotype[P1,]==-1]<-1
      allele2[genotype[P1,]==1]<-1
      allele1[genotype[P2,]==-1]<-1
      allele2[genotype[P2,]==1]<-1
      
      
      #Check that all parents are unique (this should be TRUE)
      if (length(unique(parents.top)) !=length(parents.top)){
        stop("Non unique parents selected on line 89, Rapid breeding selection")
      }
    }
    
    #update the selected alleles
    #p_select contain a TRUE if both alleles of a marker have not yet been selected in the parental selection
    p_select[allele1==1 & allele2==1]<-FALSE
    
    #available id (not yet selected parents id)
  }
 
  parents<-ind.names[parents.top]
  #check that all parents are unique
  if (length(parents)!= length(unique(parents))){
    stop("Non unique parents selected on line 188, Rapid breeding selection")
  }
  
  #processing the parents into a list
  parent.selections.i<-list(value.sel=as.matrix(value.mat[parents,]),lines.sel=parents)
  crossing.block.i<-cbind(Parent1=parents[seq(1,sel.intensity,2)],Parent2=parents[seq(2,sel.intensity,2)])
  return(list(parent.selections.i=parent.selections.i,crossing.block.i=crossing.block.i))
}
