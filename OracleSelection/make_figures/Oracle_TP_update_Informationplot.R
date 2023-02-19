library(dplyr)
library(ggplot2)
library(egg)
library(gridExtra)

map<-"Data/TP_Scoping_Train"
#map<-"own_results"
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
startbc<-c(1,1,1,1,1,1)
endbc<-c(15,15,15,15,15,15)
select.exp<-c(1,2,3,4,5,6)
exp.name<-c("Best", 'CDMean', 'PEVMean', "Random",'Tails',"Stepwise Selection")


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
      mpheno<-c()
      mvar<-c()
      dif<-c()
      mgv<-c()
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
        
          mpheno<-cbind(mpheno, gen$TP.info$meanPheno/gen$QTL_info[['maximum']])
          #mvar<-cbind(mvar, mean(gen$TP.info$MVAR))
          
          
          if(expir!=6){
            names.selected<-gen$tp.update$TP.addition.lines
            mgv<-cbind(mgv, mean(gen$candidate.values$geno.values[names.selected,])/gen$QTL_info[['maximum']])
            
            
          }else{
            dif<-cbind(dif, mean(gen$TP.info$Diff/gen$QTL_info[['maximum']]))
            names.selected<-names(gen$TP.info$Diff)
            mgv<-cbind(mgv, mean(gen$candidate.values$geno.values[names.selected,])/gen$QTL_info[['maximum']])
          }
          
         
          
          #out<-cbind(out, gen$QTL_top)
          out<-cbind(out, gen$QTL_info)
      }
      

      Data<-data.frame(t(out))
      
      Data$cycle<-seq(startbc[iter],endbc[iter])
     
     
      # nodige data in 1 column zetten met factor for ggplot legend
      #correct the breeding cycle for the preselected results
 
      cycle<-Data$cycle#+((iter-1)*5)
      
      m1<-cbind(cycle, t(mgv), t(geneticValue), t(mvar), t(dif), iter)
      m<-rbind(m1)
      m<-data.frame(m)
      colnames(m)<-c("cycle","mgv","meanGV","mvar","diff","Legend")
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
  

result.gv<-c()
result.var<-c()
result.ph<-c()
result.diff<-c()
 
for(b.c in seq(1,15)){
  
  results %>%
    filter(cycle == b.c) %>%
    group_by(Legend) %>%
    summarise(GV = mean(mgv), SD = sd(mgv, na.rm=T))->temp.cycle
  temp.cycle$Type<-"Mean GV TP"
  temp.cycle$cycle<-b.c
  result.gv<-rbind(result.gv, temp.cycle)
  
  
    results %>%
      filter(cycle == b.c) %>%
      group_by(Legend) %>%
      summarise(GV = mean(meanGV), SD = sd(meanGV))->temp.cycle
    temp.cycle$Type<-"Mean GV"
    temp.cycle$cycle<-b.c
    result.ph<-rbind(result.ph, temp.cycle)
    
    results %>%
      filter(cycle == b.c) %>%
      group_by(Legend) %>%
      summarise(GV = mean(mvar), SD = sd(mvar))->temp.cycle
    temp.cycle$Type<-"Mean Marker Variance"
    temp.cycle$cycle<-b.c
    result.var<-rbind(result.var, temp.cycle)
    
    results %>%
      filter(cycle == b.c) %>%
      group_by(Legend) %>%
      summarise(GV = mean(diff), SD = sd(diff))->temp.cycle
    temp.cycle$Type<-"Residual Noise"
    temp.cycle$cycle<-b.c
    result.diff<-rbind(result.diff, temp.cycle)
}

result.gv$Type<-paste0(result.gv$Type," (",result.gv$Legend,")")
result.var$Type<-paste0(result.var$Type," (",result.var$Legend,")")
result.ph$Type<-paste0(result.ph$Type," (",result.ph$Legend,")")
result.diff$Type<-paste0(result.diff$Type," (",result.diff$Legend,")")

title<-"Greedy parental selection with prebreeding"
result.cycle<-list(gv=result.gv, var=result.var, pheno=result.ph, diff=result.diff)
save(result.cycle, file="TP_info_scoping_100.RData")
  
  
f1<-ggplot(data=result.gv, aes(x=cycle, y=GV, colour=Type,linetype=Type)) +
    geom_line(size=1.2)  +
    scale_linetype_manual(values=c(rep("twodash",length(select.exp)),rep("solid",length(select.exp)),rep("dashed",length(select.exp))))+
    scale_colour_manual(values=rep(c(1,2,3,4,5,6,7),2))+
    ylab('Genetic Value (GV)')+
    xlab('Breeding Cycle')+
    theme_bw()+
    scale_y_continuous(breaks = seq(0, 1, 0.2),expand=c(0,0),limits=c(0,1),sec.axis = dup_axis(name = NULL, labels = NULL))+
    scale_x_continuous(breaks = seq(0, 50, 5),expand=c(0,0),limits=c(1,15),sec.axis = dup_axis(name = NULL, labels = NULL))+
    theme(legend.title =  element_blank(),text = element_text(size=12) , legend.position= c(0.55,0.8),legend.direction = "vertical")+
    guides(linetype=guide_legend(ncol=2))+
    theme(axis.text = element_text(colour="black"))+
    theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.25,"cm"))+
    theme(legend.key.width = unit(1,"cm"),legend.key.height = unit(0.8,"cm"),legend.text = element_text(size=12))+
    theme(axis.text.x =element_text(margin=unit(c(t = 0.5, r = 0, b = 0, l = 0), "cm")))+
    theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.5, b = 0, l = 0), "cm")))+
    theme(plot.margin=unit(c(1,1,0.5,0.5),"cm"))
    
f2<-ggplot(data=result.var, aes(x=cycle, y=GV, colour=Type,linetype=Type)) +
  geom_line(size=1.2)  +
  #geom_errorbar(aes(ymin=GV-SD, ymax=GV+SD), width=.2,position=position_dodge(.9))+ 
  scale_linetype_manual(values=c(rep("twodash",length(select.exp)),rep("solid",length(select.exp)),rep("dashed",length(select.exp))))+
  scale_colour_manual(values=rep(c(1,2,3,4,5,6,7),2))+
  ylab('Genetic Value (GV)')+
  xlab('Breeding Cycle')+
  theme_bw()+
  scale_y_continuous(breaks = seq(0, 1, 0.2),expand=c(0,0),limits=c(0,1),sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_x_continuous(breaks = seq(0, 50, 5),expand=c(0,0),limits=c(1,15),sec.axis = dup_axis(name = NULL, labels = NULL))+
  theme(legend.title =  element_blank(),text = element_text(size=12) , legend.position= c(0.5,0.8),legend.direction = "vertical")+
  guides(linetype=guide_legend(ncol=2))+
  theme(axis.text = element_text(colour="black"))+
  theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.25,"cm"))+
  theme(legend.key.width = unit(1,"cm"),legend.key.height = unit(0.8,"cm"),legend.text = element_text(size=12))+
  theme(axis.text.x =element_text(margin=unit(c(t = 0.5, r = 0, b = 0, l = 0), "cm")))+
  theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.5, b = 0, l = 0), "cm")))+
  theme(plot.margin=unit(c(1,1,0.5,0.5),"cm"))

f3<-ggplot(data=result.ph, aes(x=cycle, y=GV, colour=Type,linetype=Type)) +
  geom_line(size=1.2)  +
  #geom_errorbar(aes(ymin=GV-SD, ymax=GV+SD), width=.2,position=position_dodge(.9))+ 
  scale_linetype_manual(values=c(rep("twodash",length(select.exp)),rep("solid",length(select.exp)),rep("dashed",length(select.exp))))+
  scale_colour_manual(values=rep(c(1,2,3,4,5,6,7),2))+
  ylab('Genetic Value (GV)')+
  xlab('Breeding Cycle')+
  theme_bw()+
  scale_y_continuous(breaks = seq(0, 1, 0.2),expand=c(0,0),limits=c(0,1),sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_x_continuous(breaks = seq(0, 50, 5),expand=c(0,0),limits=c(1,15),sec.axis = dup_axis(name = NULL, labels = NULL))+
  theme(legend.title =  element_blank(),text = element_text(size=12) , legend.position= c(0.55,0.8),legend.direction = "vertical")+
  guides(linetype=guide_legend(ncol=2))+
  theme(axis.text = element_text(colour="black"))+
  theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.25,"cm"))+
  theme(legend.key.width = unit(1,"cm"),legend.key.height = unit(0.8,"cm"),legend.text = element_text(size=12))+
  theme(axis.text.x =element_text(margin=unit(c(t = 0.5, r = 0, b = 0, l = 0), "cm")))+
  theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.5, b = 0, l = 0), "cm")))+
  theme(plot.margin=unit(c(1,1,0.5,0.5),"cm"))
  

f4<-ggplot(data=result.diff, aes(x=cycle, y=GV, colour=Type,linetype=Type)) +
  geom_line(size=1.2)  +
  #geom_errorbar(aes(ymin=GV-SD, ymax=GV+SD), width=.2,position=position_dodge(.9))+ 
  scale_linetype_manual(values=c(rep("twodash",length(select.exp)),rep("solid",length(select.exp)),rep("dashed",length(select.exp))))+
  scale_colour_manual(values=rep(c(1,2,3,4,5,6,7),2))+
  ylab('Genetic Value (GV)')+
  xlab('Breeding Cycle')+
  theme_bw()+
  scale_y_continuous(breaks = seq(0, 1, 0.2),expand=c(0,0),limits=c(0,1),sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_x_continuous(breaks = seq(0, 50, 5),expand=c(0,0),limits=c(1,15),sec.axis = dup_axis(name = NULL, labels = NULL))+
  theme(legend.title =  element_blank(),text = element_text(size=12) , legend.position= c(0.55,0.23),legend.direction = "vertical")+
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





