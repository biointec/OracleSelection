#make figure of Oracle TP update methods. Number of additions, deletions and TP size

library(dplyr)
library(ggplot2)
library(egg)
library(gridExtra)

#map<-"Data/TP_Scoping_Train"
map<-"Data/TP_Train_100"
if (file.exists(paste0(map,"/Results"))){
  unlink(paste0(map,"/Results"), recursive=TRUE)
}

files<-dir(map, full.names = T)
names<-dir(map, full.names = F)
results<-c()
level<-c()
fac<-c()
dir.create(paste0(map,"/Results"))

#load("own_results/Baseline_method_100reps.RData")
#baseline<-experiment.sub.results

iter<-0
startbc<-c(1)
endbc<-c(15)
select.exp<-c(7)
exp.name<-c("Oracle")


Useable_exp<-matrix(T,100,1)
#iter=0
#for (i in c(select.exp)){
#  iter=iter+1
#  load(files[i])
#  n.set<-length(experiment.sub.results)
#  for(s in seq(n.set)){
#    set <- experiment.sub.results[[s]]
#    if (!is.null(set)  & (class(set) != "try-error")){
#      nrep<-length(set)
#      repnames<-names(set)
#      for(r in seq(nrep)){ 
#        rep<-set[[r]]
#        if (!any(is.na(rep)) ){
#          rname<-as.numeric(substr(repnames[r],4,nchar(repnames[r])))
#          Useable_exp[ rname,iter]<-T
#        }
#      }
#    }
#  }
#}

#Useable_exp<-Useable_exp[,1] & Useable_exp[,2]& Useable_exp[,3]& Useable_exp[,4]
iter=0


for(expir in select.exp){
  iter<-iter+1
  load(files[expir])
  a<-names[expir]
  n.set<-length(experiment.sub.results)
  
  acc_fixed<-c()
  for(s in seq(n.set)){
    set <- experiment.sub.results[[s]]
    repnames<-names(set)
    if(!is.null(set)  & (class(set) != "try-error")){
    nrep<-length(set)
    for(r in seq(nrep)){ 
      rep<-set[[r]]
      rname<-as.numeric(substr(repnames[r],4,nchar(repnames[r])))
      if (!any(is.na(rep))){#& Useable_exp[rname]){
      genome<-rep$genome
      sim<-rep$sim.results
      
     
      out<-c()    
      geneticValue<-c()
      geneticValue.top<-c()
      TP_size<-c()
      add<-c()
      remove<-c()
      tot<-c()
     # for(br.cycle in seq(metadata$n.cycles)){
      for(br.cycle in seq(startbc[iter],endbc[iter])){
      #  if(br.cycle<=startbc[iter]){
        #  gen<-rep_baseline$sim.results[[br.cycle]]
      #  }else{
          cycle.name<-paste0("cycle",br.cycle)
          gen<-sim[[cycle.name]]
      #  }
          candidate.af.i<-gen$geno.summary.stats$candidate.af
          geneticValue<-cbind(geneticValue, gen$candidate.values$mu.g/gen$QTL_info[['maximum']])
          geneticValue.top<-cbind(geneticValue.top, mean(sort(gen$candidate.values$geno.values,decreasing = TRUE)[1:10])/gen$QTL_info[['maximum']])
        
          add<-cbind(add, gen$lengthTP$Added)
          remove<-cbind(remove, gen$lengthTP$Removed)
          tot<- cbind(tot, gen$lengthTP$Length)
        
         
          #out<-cbind(out, gen$QTL_top)
          out<-cbind(out, gen$QTL_info)
          
          
          
      }
      

      Data<-data.frame(t(out))
      
      Data$cycle<-seq(startbc[iter],endbc[iter])
     
     
      # nodige data in 1 column zetten met factor for ggplot legend
      #correct the breeding cycle for the preselected results
 
      cycle<-Data$cycle#+((iter-1)*5)
      
      m1<-cbind(cycle, t(add), t(remove), t(tot), iter)
      m<-rbind(m1)
      m<-data.frame(m)
      colnames(m)<-c("cycle","Add","Remove","TPsize","Legend")
      results<-rbind(results, m)
      }
     }
   }
  }
  level<-c(level, toString(iter))
  
  fac<-c(fac,exp.name[iter])

}#Close the for expir loop
results$Legend<-factor(results$Legend, labels = fac)


#bereken de gemiddelde per generatie en de std om deze dan te plotten
sum.result<-data.frame(cycle = seq(metadata$n.cycles))
  

result.cycle<-c()
 
for(b.c in seq(1,15)){
    results %>%
      filter(cycle == b.c) %>%
      group_by(Legend) %>%
      summarise(GV = mean(Add), SD = sd(Add))->temp.cycle
    temp.cycle$Type<-"Add"
    temp.cycle$cycle<-b.c
    result.cycle<-rbind(result.cycle, temp.cycle)
    
    results %>%
      filter(cycle == b.c) %>%
      group_by(Legend) %>%
      summarise(GV = mean(Remove), SD = sd(Remove))->temp.cycle
    temp.cycle$Type<-"Remove"
    temp.cycle$cycle<-b.c
    result.cycle<-rbind(result.cycle, temp.cycle)
    
    results %>%
      filter(cycle == b.c) %>%
      group_by(Legend) %>%
      summarise(GV = mean(TPsize), SD = sd(TPsize))->temp.cycle
    temp.cycle$Type<-"TP size"
    temp.cycle$cycle<-b.c
    result.cycle<-rbind(result.cycle, temp.cycle)
}

result.cycle$Type<-paste0(result.cycle$Type," (",result.cycle$Legend,")")

title<-"Greedy parental selection with prebreeding"

save("result.cycle", file="TP_Truncation_add.RData")
  
  
f1<-ggplot(data=result.cycle, aes(x=cycle, y=GV, colour=Type,linetype=Type)) +
    geom_line(size=1.2)  +
    #geom_errorbar(aes(ymin=GV-SD, ymax=GV+SD), width=.2,position=position_dodge(.9))+ 
    scale_linetype_manual(values=c(rep("twodash",length(select.exp)),rep("solid",length(select.exp)),rep("dashed",length(select.exp))))+
    scale_colour_manual(values=rep(c(1,2,3,4,5,6,7),2))+
    ylab('Genetic Value (GV)')+
    xlab('Breeding Cycle')+
    theme_bw()+
    scale_y_continuous(breaks = seq(0, 300, 50),expand=c(0,0),limits=c(0,250),sec.axis = dup_axis(name = NULL, labels = NULL))+
    scale_x_continuous(breaks = seq(0, 50, 5),expand=c(0,0),limits=c(1,15),sec.axis = dup_axis(name = NULL, labels = NULL))+
    theme(legend.title =  element_blank(),text = element_text(size=12))+
    guides(linetype=guide_legend(ncol=2))+
    theme(axis.text = element_text(colour="black"))+
    theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.25,"cm"))+
    theme(legend.key.width = unit(1,"cm"),legend.key.height = unit(0.8,"cm"),legend.text = element_text(size=12))+
    theme(axis.text.x =element_text(margin=unit(c(t = 0.5, r = 0, b = 0, l = 0), "cm")))+
    theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.5, b = 0, l = 0), "cm")))+
    theme(plot.margin=unit(c(1,1,0.5,0.5),"cm"))
    
  
  


titel<-"Effects of prebreeding on a greedy selection"
filename<-paste0(map,"/Results/",titel,".eps")
ggsave(filename, height=8 ,width= 10, f1,device="eps")





