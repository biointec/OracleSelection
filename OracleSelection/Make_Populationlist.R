load("Population/Population_H80_file.RData")
Population.list<-list()
n.set<-length(Simulated.populations)
for (i in seq(n.set)){
  set<-Simulated.populations[[i]]
  n.rep<-length(set)
  repnames<-names(set)
  for(j in seq(n.rep)){
    rep<-set[[j]]
    pop<-as.numeric(substr(repnames[j],4,nchar(repnames[j])))
    if(!is.na(rep)){
    Population.list[[pop]]<-list(Population=rep$population, genome=rep$genome)
    }else{
    Population.list[[pop]]<-list(Population=NA, genome=NA)
    }
  }
}
save("Population.list",file="Population/Population_DeepScoping_H08.RData")


###### Part 2: save the truncation selection
experiment.sub.results<-list()

for (i in seq(n.set)){
  set<-Simulated.populations[[i]]
  n.rep<-length(set)
  repnames<-names(set)
  for(j in seq(n.rep)){
    rep<-set[[j]]
    pop<-as.numeric(substr(repnames[j],4,nchar(repnames[j])))
    if(!any(is.na(rep))){
      sim.results<-rep$simulation
      experiment.sub.results[[paste0("set",i)]][[repnames[j]]]<-list(sim.results=sim.results,genome=rep$genome)
    }
  }
}

save("experiment.sub.results", file="Population/Truncation_selection_BC20_H08.RData")
