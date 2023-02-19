Crossing.Scoping<-function(candidates_donor=SelDonors4UC2,
                           candidates_popE=Poptop,
                           poplayer=PopD2,
                           Z=Geno,
                           GEBV,
                           layers){
  
  
  #get size
  n.pop<-dim(Z)[1]
  k.pop<-dim(Z)[2]
  
  #pre allocate 
  marker_allele=matrix(0,2,k.pop)
  marker_seq=matrix(TRUE,1,k.pop)
  CrosBlock<-list()
  preselect=0
 
  #for each layer ...	Layer 4 is first selected, then layer 3, 2, ...
  for (Layer in seq(layers,1,-1)){
    PopL<-c()
    if (Layer< layers){
      for (i in seq(layers-Layer)){
        PopL<-c(PopL,poplayer[[Layer+i]])
      }
    }
    
    selected.names<-c(candidates_popE[[Layer]],PopL)
    
    #delete individuals that have already been selected in previous rounds
    if (length(CrosBlock)!=0){
      already.selected<-c(as.character(CrosBlock[[1]][,1]),as.character(CrosBlock[[1]][,2]))
    for (i in length(already.selected)){
      selected.names<-selected.names[!selected.names==already.selected[i]]
    }
  }
    
    
    
  ordr<-order(GEBV[selected.names,],decreasing=T)
  #get names and values
  ind.names <- selected.names[ordr]
  ind.values<- as.matrix(GEBV[selected.names,])[ordr,]  

    
  Parents_available<-ind.names
    
  #Increase the threshold of the pre-selection
  preselect<-preselect+100
  
  free.position<-seq(preselect)
  all_parents<-c()

  #select first parent------------------------------------------------------------
  name1<-candidates_donor[[Layer]][1]
  gen1<-Z[name1,]
  all_parents<-rbind(all_parents,Z[name1,])
  p1<- name1
  
  score<-0
  sel<-0
  #Select first parent of that layer
  if (Layer==layers){ #if this is your first layer
    
  for (i in free.position){
    temp.name<- Parents_available[i]
    temp.geno<- as.matrix(Z[temp.name,])
    temp.score<-sum(colVars(rbind(all_parents, t(temp.geno)))[marker_seq]) #S-Score
    if (temp.score == 0 ){ 
      temp.score<-sum(colVars(rbind(all_parents, t(temp.geno))))
    }
    
    if (temp.score>score){
      score<-temp.score
      sel<-i
    }
  }  
    chosen<- sel
  }else{
    for (i in free.position[-chosen]){
      temp.name<- Parents_available[i]
      temp.geno<- as.matrix(Z[temp.name,])
      temp.score<-sum(colVars(rbind(all_parents, t(temp.geno)))[marker_seq]) #S-Score
      if (temp.score == 0 ){ 
        temp.score<-sum(colVars(rbind(all_parents, t(temp.geno))))
      }
      
      if (temp.score>score){
        score<-temp.score
        sel<-i
      }
    }  
    chosen<-c(chosen,sel)
  }
  
  p2<-c(Parents_available[sel])
  name2<-Parents_available[sel]
  gen2<-Z[sel,]
  
  all_parents<-rbind(all_parents, Z[p2,])
  
  # Keep track of the selected alleles
  marker_allele[1,(gen1+gen2)==-2]<-marker_allele[1,(gen1+gen2)==-2]+1
  marker_allele[2,(gen1+gen2)== 2]<-marker_allele[2,(gen1+gen2)==2]+1
 
  marker_allele[,abs(gen1+gen2)== 1]<-marker_allele[,abs(gen1+gen2)== 1]+1
  marker_seq[marker_allele[1,]>0 & marker_allele[2,]>0]<- FALSE
  
  #Selected the remaining parents for that layer
  for (sel.round in seq(2,length(candidates_donor[[Layer]]))){
    
    name1<-candidates_donor[[Layer]][sel.round]
    gen1<-Z[name1,]
    all_parents<-rbind(all_parents, Z[name1,])
    
    
    p1<-c(p1, name1)
    geno1<- as.matrix(Z[name1,])
    score<-0
    sel<-0
    for (i in free.position[-chosen]){
      temp.name<-Parents_available[i]
      temp.geno<- as.matrix(Z[temp.name,])
      
      temp.score<-sum(colVars(rbind(all_parents, t(temp.geno)))[marker_seq]) #S-Score
      if (temp.score == 0 ){ 
        temp.score<-sum(colVars(rbind(all_parents, t(temp.geno))))
      }
      
      if (temp.score>score){
        score<-temp.score
        sel<-i
      }
    }
    
    chosen<-c(chosen, sel)
    p2<-c(p2, Parents_available[sel])
    name2<-Parents_available[sel]
    gen2<-Z[name2,]
    #Keep track of selected alleles
    marker_allele[1,(gen1+gen2)==-2]<-marker_allele[1,(gen1+gen2)==-2]+1
    marker_allele[2,(gen1+gen2)== 2]<-marker_allele[2,(gen1+gen2)==2]+1
   
    marker_allele[,abs(gen1+gen2)== 1]<-marker_allele[,abs(gen1+gen2)== 1]+1
    
    marker_seq[marker_allele[1,]>0 & marker_allele[2,]>0]<- FALSE
    
    all_parents<-rbind(all_parents, Z[Parents_available[sel],])
    score<-temp.score
    sel<-i
    
  }
  
  

  p.values<-as.matrix(c(GEBV[p1,],GEBV[p2,]))
  p.names<-rownames(p.values)
  
  #construct crossing block
  CrosBlock[[layers-Layer+1]]<-as.data.frame(cbind(p1,p2))
  }
  
  return(rev(CrosBlock))
}
