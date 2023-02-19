## Calculate OCS-score based on allier et al.
#-------------------------------------------------------------------------
#Copyright (c) 2019 Antoine Allier
#-------------------------------------------------------------------------

OGM_H<-function(candidate.marker.genos.i,
                        candidate.haplotype=candidate.haploid.i,
                        Bhat=predictions.out$solve.out$u,
                        pheno.train,
                        p.sel=100,
                        presel=300,
                        A.sc,
                        parents.information=parent.selections.list,
                        trainingpanel=list(geno=TP.genos.i,pheno=TP.phenos.i),
                        map=map,
                        rep.iter,
                        genome=hv.genome){
  
  Parent.genos <- genotype.loci(haploid.genos = candidate.haplotype, 
                                genome = genome, 
                                include.QTL = F)
  
  
  GEBV <-   Parent.genos %*% Bhat
  
  
  
  # Assuming the Elites are the 10 best lines
  PopE <- names(GEBV[order(GEBV,decreasing = TRUE),])[seq(p.sel/2)]
  
  preselect <- names(GEBV[order(GEBV,decreasing = TRUE),])[seq(presel)]

  PopD <- sample(setdiff(preselect,PopE),(length(preselect)-p.sel/2),replace = FALSE)
  
 
  
      
      UCtmp <-getGaSolutionsFrontier(
        Markers=Parent.genos[PopE,],
        Markers2=Parent.genos[PopD,],
        K=A.sc[c(PopE,PopD),c(PopE,PopD)],
        markereffects=Bhat,
        markermap=NULL,
        nmates=p.sel/2,
        npopGA=100,
        nitGA=100,
        mc.cores=1,
        mutprob=0.999,
        noself=TRUE,
        method=3,
        type=2L,
        generation=1L,
        plotiters=F)
      
      n.iters=length(UCtmp[[2]])
      
  cros<-UCtmp[[2]][[n.iters]]
  colnames(cros)<-c("Parent1","Parent2")
  
  selection<-c(cros[,1], cros[,2] )


  parent.sel<-c()
  parent.sel$lines.sel<-unique(selection)
  
  parent.sel$value.sel<-GEBV[unique(selection),]
  
  return(list(parent.selections.i=parent.sel,crossing.block.i=cros))
  
}
