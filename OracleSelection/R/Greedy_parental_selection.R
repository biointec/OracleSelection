#' Rapid breeding selection
#' @value.mat contains the gebv of all the individuals
#' @sel.intensity The number of individuals that have te be selected as parent
#' @genos The genotype of all the individuals
#' @fraction_greedy The fraction that need to be selected greedily
#' @G the relationship matrix calculated according to VanRaden
#' 
#' This method selects parents and return a list with the selected parents and crossing block
Rapid.parental.selection<-function(value.mat = predictions.out$GEBV, 
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
  gebv<-as.matrix((value.mat-mean(value.mat))/sd(value.mat))
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
  list.top<-order(gebv,decreasing=TRUE)
  
  
  
  
  
 
  
  #processing the parents into a list
  parent.selections.i<-list(value.sel=as.matrix(value.mat[parents,]),lines.sel=parents)
  crossing.block.i<-cbind(Parent1=parents[seq(1,sel.intensity,2)],Parent2=parents[seq(2,sel.intensity,2)])
  return(list(parent.selections.i=parent.selections.i,crossing.block.i=crossing.block.i))
}