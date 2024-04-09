library(ggplot2)
library(viridis)
library(cowplot)
library(ggpubr)

###################################################################
sample_and_NLL = readRDS("sample_and_NLL.RDS")
resample = readRDS("resample.RDS")
resample_sample = sample_and_NLL[resample,]
index= unique(resample)
###################################################################

prevalences0 = readRDS("Bprevalences_factual.RDS")


prevC0_0    = prevalences0[,1]
prevH0_0    = prevalences0[,5]
prevN0_0    = prevalences0[,9]
prevWSHN0_0 = prevalences0[,13]
prevS0_0    = prevalences0[,17]
prevW0_0    = prevalences0[,21]
prevWSH0_0  = prevalences0[,25]

prevC12_0    = prevalences0[,4]
prevH12_0    = prevalences0[,8]
prevN12_0    = prevalences0[,12]
prevWSHN12_0 = prevalences0[,16]
prevS12_0    = prevalences0[,20]
prevW12_0    = prevalences0[,24]
prevWSH12_0  = prevalences0[,28]


#Effectiveness
W0=1-prevW12_0/prevC12_0
H0=1-prevH12_0/prevC12_0
S0=1-prevS12_0/prevC12_0
WSH0=1-prevWSH12_0/prevC12_0
N0=1-prevN12_0/prevC12_0
WSHN0=1-prevWSHN12_0/prevC12_0

print(100*c(median((W0)),
        median((S0)),
        median((H0)),
        median((WSH0)),
        median((N0)),
        median((WSHN0))))

###################################################################
rel_rel_risks = function(file,color=NULL){
  
  prevalences1= readRDS(file)
  

  
  prevC0_1    = prevalences1[,1]
  prevH0_1    = prevalences1[,5]
  prevN0_1    = prevalences1[,9]
  prevWSHN0_1 = prevalences1[,13]
  prevS0_1    = prevalences1[,17]
  prevW0_1    = prevalences1[,21]
  prevWSH0_1  = prevalences1[,15]
  
  prevC12_1    = prevalences1[,4]
  prevH12_1    = prevalences1[,8]
  prevN12_1    = prevalences1[,12]
  prevWSHN12_1 = prevalences1[,16]
  prevS12_1    = prevalences1[,20]
  prevW12_1    = prevalences1[,24]
  prevWSH12_1  = prevalences1[,28]
  
  #Effectiveness
  W1=1-prevW12_1/prevC12_1
  H1=1-prevH12_1/prevC12_1
  S1=1-prevS12_1/prevC12_1
  WSH1=1-prevWSH12_1/prevC12_1
  N1=1-prevN12_1/prevC12_1
  WSHN1=1-prevWSHN12_1/prevC12_1
  
  #Disregard points where there is no disease in the counterfactual scenario
  index=which(prevC0_1==1E-10)
  plot_sample=setdiff(1:nrow(prevalences0),index)

  y = 100*c((W1[plot_sample])-(W0[plot_sample]),
          (S1[plot_sample])-(S0[plot_sample]),
          (H1[plot_sample])-(H0[plot_sample]),
          (WSH1[plot_sample])-(WSH0[plot_sample]),
          (N1[plot_sample])-(N0[plot_sample]),
          (WSHN1[plot_sample])-(WSHN0[plot_sample]))
  
  print(c(median(100*(W1)),
          median(100*(S1)),
          median(100*(H1)),
          median(100*(WSH1)),
          median(100*(N1)),
          median(100*(WSHN1))))
  print(c(median(100*(W1-W0)),
          median(100*(S1-S0)),
          median(100*(H1-H0)),
          median(100*(WSH1-WSH0)),
          median(100*(N1-N0)),
          median(100*(WSHN1-WSHN0))))
  
  Arm = rep(c("W","S","H","WSH","N","WSHN"),each=length(plot_sample))
  
  Color=NULL
  
  if(color=="water pathway"){Color=rep(resample_sample[plot_sample,2],6)}
  if(color=="fomite pathway"){Color=rep((1-resample_sample[plot_sample,2])*(resample_sample[plot_sample,3]),6)}
  if(color=="other pathway"){Color=rep((1-resample_sample[plot_sample,2])*(1-resample_sample[plot_sample,3]),6)}
  if(color=="water"){Color=rep(resample_sample[plot_sample,5],6)}
  if(color=="sanitation"){Color=rep(resample_sample[plot_sample,7],6)}
  if(color=="hygiene"){Color=rep(resample_sample[plot_sample,6],6)}
  if(color=="nutrition"){Color=rep(resample_sample[plot_sample,4],6)}
  
  DATA=data.frame(y,Arm,Color)
  DATA$Arm=factor(DATA$Arm,levels=c("W","S","H","WSH","N","WSHN"))
  DATA$Quantile = cut(DATA$Color, breaks=quantile(DATA$Color, probs=c(0:5/5)), labels=1:5,
                      include.lowest=TRUE)
  DATA$QuantileArm=apply(DATA[,c("Arm","Quantile")],1,paste,collapse="")
  DATA$QuantileArm=factor(DATA$QuantileArm,levels=c("W1","W2","W3","W4","W5",
                                                    "S1","S2","S3","S4","S5",
                                                    "H1","H2","H3","H4","H5",
                                                    "WSH1","WSH2","WSH3","WSH4","WSH5",
                                                    "N1","N2","N3","N4","N5",
                                                    "WSHN1","WSHN2","WSHN3","WSHN4","WSHN5"))
  
  return(DATA)
}


### Now plots counterfactual vs factual relative risk

DATA=rel_rel_risks("Bprevalences_no_conditions_with_var.RDS",color="water pathway")
write.csv(DATA,"Table_Fig3A.csv")
p1=ggplot(DATA, aes(x=Arm, y=y,fill=Arm))+
  ggtitle("No preexisting WASH conditions")+
  scale_fill_manual(guide="none",values=c("#117733", "#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(bg=Arm))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-100,10))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p1

DATA=rel_rel_risks("Bprevalences_2prev_with_var.RDS",color="water pathway")
write.csv(DATA,"Table_Fig3B.csv")
p2=ggplot(DATA, aes(x=Arm, y=y,fill=Arm))+
  ggtitle("Double baseline disease prevalence")+
  scale_fill_manual(guide="none",values=c("#117733", "#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(bg=Arm))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-100,10))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p2

DATA=rel_rel_risks("Bprevalences_full_adherence_with_var.RDS",color="water pathway")
write.csv(DATA,"Table_Fig3C.csv")
p3=ggplot(DATA, aes(x=Arm, y=y,fill=Arm))+
  ggtitle("Full compliance")+
  scale_fill_manual(guide="none",values=c("#117733", "#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(bg=Arm))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p3

DATA=rel_rel_risks("Bprevalences_50_less_other_with_var.RDS",color="water pathway")
write.csv(DATA,"Table_Fig3D.csv")
p4=ggplot(DATA, aes(x=Arm, y=y,fill=Arm))+
  ggtitle("Half of other pathway transmission can be intervened on")+
  scale_fill_manual(guide="none",values=c("#117733", "#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(bg=Arm))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p4

DATA=rel_rel_risks("Bprevalences_100efficacy_W_with_var.RDS",color="water pathway")
write.csv(DATA,"Table_Fig3E.csv")
p5=ggplot(DATA, aes(x=Arm, y=y,fill=Arm))+
  ggtitle("Double efficacy of water chlorination")+
  scale_fill_manual(guide="none",values=c("#117733", "#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(bg=Arm))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p5

DATA=rel_rel_risks("Bprevalences_100efficacy_S_with_var.RDS",color="water pathway")
write.csv(DATA,"Table_Fig3F.csv")
p6=ggplot(DATA, aes(x=Arm, y=y,fill=Arm))+
  ggtitle("Double efficacy of a latrine water seal")+
  scale_fill_manual(guide="none",values=c("#117733", "#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(bg=Arm))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p6

DATA=rel_rel_risks("Bprevalences_100efficacy_H_with_var.RDS",color="water pathway")
write.csv(DATA,"Table_Fig3G.csv")
p7=ggplot(DATA, aes(x=Arm, y=y,fill=Arm))+
  ggtitle("Double efficacy of handwashing with soap and water")+
  scale_fill_manual(guide="none",values=c("#117733", "#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(bg=Arm))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p7

DATA=rel_rel_risks("Bprevalences_100efficacy_N_with_var.RDS",color="water pathway")
write.csv(DATA,"Table_Fig3H.csv")
p8=ggplot(DATA, aes(x=Arm, y=y,fill=Arm))+
  ggtitle("Double efficacy of nutrition supplementation")+
  scale_fill_manual(guide="none",values=c("#117733", "#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(bg=Arm))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p8

plot_list= list(p1,p2,p3,p4,p5,p6,p7,p8)
plot_grid(plotlist=plot_list,labels="auto",nrow=4, ncol=2,align="v",axis="lr")
ggsave("Figure 3.pdf",height=8, width=6.0)


######### Coverage
DATA=rel_rel_risks("Bprevalences_coverage20_with_var.RDS",color="water pathway")
write.csv(DATA,"Table_Fig4A.csv")
p1=ggplot(DATA, aes(x=Arm, y=y,fill=Arm))+
  ggtitle("Increase community coverage to 20%")+
  scale_fill_manual(guide="none",values=c("#117733", "#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(bg=Arm))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p1

DATA=rel_rel_risks("Bprevalences_coverage20_with_var.RDS",color="water pathway")
write.csv(DATA,"Table_Fig4B.csv")
p2=ggplot(DATA, aes(x=Arm, y=y,fill=QuantileArm,Group=Arm))+
  ggtitle("by of the strength of the water pathway")+
  scale_fill_manual(guide="none",values=c("#1dc956","#169c43","#117733", "#0d5926", "#0a431d",
                                          "#7fccbf","#5bbead","#44AA99","#3a9283","#2c6d62",
                                          "#d2ecf9","#a6d9f2","#88CCEE","#63bce9","#36a9e2",
                                          "#f0e8c2","#e6d999","#DDCC77","#d6c25c","#ccb333",
                                          "#e6b3bb","#d98c99","#CC6677","#bf4055","#993344",
                                          "#cc7fbf","#be5bad","#AA4499","#923a83","#6d2c62"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=QuantileArm,bg=QuantileArm),position = position_dodge(width=0.9))+
  # stat_summary(fun="median",geom="line",linewidth=0.1,aes(group=Arm,x=dodge(Arm,width=0.9)))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p2

DATA=rel_rel_risks("Bprevalences_coverage20_with_var.RDS",color="fomite pathway")
write.csv(DATA,"Table_Fig4C.csv")
p3=ggplot(DATA, aes(x=Arm, y=y,fill=QuantileArm,Group=Arm))+
  ggtitle("by of the strength of the fomite pathway")+
  scale_fill_manual(guide="none",values=c("#1dc956","#169c43","#117733", "#0d5926", "#0a431d",
                                          "#7fccbf","#5bbead","#44AA99","#3a9283","#2c6d62",
                                          "#d2ecf9","#a6d9f2","#88CCEE","#63bce9","#36a9e2",
                                          "#f0e8c2","#e6d999","#DDCC77","#d6c25c","#ccb333",
                                          "#e6b3bb","#d98c99","#CC6677","#bf4055","#993344",
                                          "#cc7fbf","#be5bad","#AA4499","#923a83","#6d2c62"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=QuantileArm,bg=QuantileArm),position = position_dodge(width=0.9))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p3

DATA=rel_rel_risks("Bprevalences_coverage20_with_var.RDS",color="other pathway")
write.csv(DATA,"Table_Fig4D.csv")
p4=ggplot(DATA, aes(x=Arm, y=y,fill=QuantileArm,Group=Arm))+
  ggtitle("by of the strength of all other pathways")+
  scale_fill_manual(guide="none",values=c("#1dc956","#169c43","#117733", "#0d5926", "#0a431d",
                                          "#7fccbf","#5bbead","#44AA99","#3a9283","#2c6d62",
                                          "#d2ecf9","#a6d9f2","#88CCEE","#63bce9","#36a9e2",
                                          "#f0e8c2","#e6d999","#DDCC77","#d6c25c","#ccb333",
                                          "#e6b3bb","#d98c99","#CC6677","#bf4055","#993344",
                                          "#cc7fbf","#be5bad","#AA4499","#923a83","#6d2c62"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=QuantileArm,bg=QuantileArm),position = position_dodge(width=0.9))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p4

DATA=rel_rel_risks("Bprevalences_coverage20_with_var.RDS",color="water")
write.csv(DATA,"Table_Fig4E.csv")
p5=ggplot(DATA, aes(x=Arm, y=y,fill=QuantileArm,Group=Arm))+
  ggtitle("by of the efficacy of the water intervention")+
  scale_fill_manual(guide="none",values=c("#1dc956","#169c43","#117733", "#0d5926", "#0a431d",
                                          "#7fccbf","#5bbead","#44AA99","#3a9283","#2c6d62",
                                          "#d2ecf9","#a6d9f2","#88CCEE","#63bce9","#36a9e2",
                                          "#f0e8c2","#e6d999","#DDCC77","#d6c25c","#ccb333",
                                          "#e6b3bb","#d98c99","#CC6677","#bf4055","#993344",
                                          "#cc7fbf","#be5bad","#AA4499","#923a83","#6d2c62"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=QuantileArm,bg=QuantileArm),position = position_dodge(width=0.9))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p5

DATA=rel_rel_risks("Bprevalences_coverage20_with_var.RDS",color="sanitation")
write.csv(DATA,"Table_Fig4F.csv")
p6=ggplot(DATA, aes(x=Arm, y=y,fill=QuantileArm,Group=Arm))+
  ggtitle("by of the efficacy of the sanitation intervention")+
  scale_fill_manual(guide="none",values=c("#1dc956","#169c43","#117733", "#0d5926", "#0a431d",
                                          "#7fccbf","#5bbead","#44AA99","#3a9283","#2c6d62",
                                          "#d2ecf9","#a6d9f2","#88CCEE","#63bce9","#36a9e2",
                                          "#f0e8c2","#e6d999","#DDCC77","#d6c25c","#ccb333",
                                          "#e6b3bb","#d98c99","#CC6677","#bf4055","#993344",
                                          "#cc7fbf","#be5bad","#AA4499","#923a83","#6d2c62"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=QuantileArm,bg=QuantileArm),position = position_dodge(width=0.9))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p6

DATA=rel_rel_risks("Bprevalences_coverage20_with_var.RDS",color="hygiene")
write.csv(DATA,"Table_Fig4G.csv")
p7=ggplot(DATA, aes(x=Arm, y=y,fill=QuantileArm,Group=Arm))+
  ggtitle("by of the efficacy of the hygiene intervention")+
  scale_fill_manual(guide="none",values=c("#1dc956","#169c43","#117733", "#0d5926", "#0a431d",
                                          "#7fccbf","#5bbead","#44AA99","#3a9283","#2c6d62",
                                          "#d2ecf9","#a6d9f2","#88CCEE","#63bce9","#36a9e2",
                                          "#f0e8c2","#e6d999","#DDCC77","#d6c25c","#ccb333",
                                          "#e6b3bb","#d98c99","#CC6677","#bf4055","#993344",
                                          "#cc7fbf","#be5bad","#AA4499","#923a83","#6d2c62"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=QuantileArm,bg=QuantileArm),position = position_dodge(width=0.9))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p7

DATA=rel_rel_risks("Bprevalences_coverage20_with_var.RDS",color="nutrition")
write.csv(DATA,"Table_Fig4H.csv")
p8=ggplot(DATA, aes(x=Arm, y=y,fill=QuantileArm,Group=Arm))+
  ggtitle("by of the efficacy of the nutrition intervention")+
  scale_fill_manual(guide="none",values=c("#1dc956","#169c43","#117733", "#0d5926", "#0a431d",
                                          "#7fccbf","#5bbead","#44AA99","#3a9283","#2c6d62",
                                          "#d2ecf9","#a6d9f2","#88CCEE","#63bce9","#36a9e2",
                                          "#f0e8c2","#e6d999","#DDCC77","#d6c25c","#ccb333",
                                          "#e6b3bb","#d98c99","#CC6677","#bf4055","#993344",
                                          "#cc7fbf","#be5bad","#AA4499","#923a83","#6d2c62"))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1)+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=QuantileArm,bg=QuantileArm),position = position_dodge(width=0.9))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  scale_y_continuous(expand = c(0,0),limits=c(-10,100))+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6))
p8

plot_list= list(p1,p2,p3,p4,p5,p6,p7,p8)
plot_grid(plotlist=plot_list,labels="auto",nrow=4, ncol=2,align="v",axis="lr")
ggsave("Figure 4.pdf",height=8, width=6.0)

##########################################################################

DATA=NULL

for (i in 2:10){
  prevalences0 = readRDS("Bprevalences_factual.RDS")
  prevalences1 = readRDS(paste0("Bprevalences_coverage",i*10,"_with_var.RDS"))
  
  prevC0_0    = prevalences0[,1]
  prevH0_0    = prevalences0[,5]
  prevN0_0    = prevalences0[,9]
  prevWSHN0_0 = prevalences0[,13]
  prevS0_0    = prevalences0[,17]
  prevW0_0    = prevalences0[,21]
  prevWSH0_0  = prevalences0[,25]
  
  prevC12_0    = prevalences0[,4]
  prevH12_0    = prevalences0[,8]
  prevN12_0    = prevalences0[,12]
  prevWSHN12_0 = prevalences0[,16]
  prevS12_0    = prevalences0[,20]
  prevW12_0    = prevalences0[,24]
  prevWSH12_0  = prevalences0[,28]
  
  prevC0_1    = prevalences1[,1]
  prevH0_1    = prevalences1[,5]
  prevN0_1    = prevalences1[,9]
  prevWSHN0_1 = prevalences1[,13]
  prevS0_1    = prevalences1[,17]
  prevW0_1    = prevalences1[,21]
  prevWSH0_1  = prevalences1[,15]
  
  prevC12_1    = prevalences1[,4]
  prevH12_1    = prevalences1[,8]
  prevN12_1    = prevalences1[,12]
  prevWSHN12_1 = prevalences1[,16]
  prevS12_1    = prevalences1[,20]
  prevW12_1    = prevalences1[,24]
  prevWSH12_1  = prevalences1[,28]
  
  W0=1-prevW12_0/prevC12_0
  W1=1-prevW12_1/prevC12_1
  H0=1-prevH12_0/prevC12_0
  H1=1-prevH12_1/prevC12_1
  S0=1-prevS12_0/prevC12_0
  S1=1-prevS12_1/prevC12_1
  WSH0=1-prevWSH12_0/prevC12_0
  WSH1=1-prevWSH12_1/prevC12_1
  N0=1-prevN12_0/prevC12_0
  N1=1-prevN12_1/prevC12_1
  WSHN0=1-prevWSHN12_0/prevC12_0
  WSHN1=1-prevWSHN12_1/prevC12_1
  
  index=which(prevC0_1==1E-10)
  plot_sample=setdiff(1:nrow(prevalences0),index)
  Color = (1-resample_sample[plot_sample,2])*(1-resample_sample[plot_sample,3]) #fraction other pathway
  plot_sample2=plot_sample[which(Color<quantile(Color,0.20))]
  plot_sample3=plot_sample[which(Color>quantile(Color,0.80))]
  
  y2 = cbind(c(W1[plot_sample2]-W0[plot_sample2],
              S1[plot_sample2]-S0[plot_sample2],
              H1[plot_sample2]-H0[plot_sample2],
              WSH1[plot_sample2]-WSH0[plot_sample2],
              N1[plot_sample2]-N0[plot_sample2],
              WSHN1[plot_sample2]-WSHN0[plot_sample2]),
              rep(i*10,length(plot_sample2)),
              rep(c("W","S","H",
                "WSH","N","WSHN"),
                each=length(plot_sample2)),
              rep("Highest",times=6*length(plot_sample2)))


  y3 = cbind(c(W1[plot_sample3]-W0[plot_sample3],
                 S1[plot_sample3]-S0[plot_sample3],
                 H1[plot_sample3]-H0[plot_sample3],
                 WSH1[plot_sample3]-WSH0[plot_sample3],
                 N1[plot_sample3]-N0[plot_sample3],
                 WSHN1[plot_sample3]-WSHN0[plot_sample3]),
               rep(i*10,length(plot_sample3)),
               rep(c("W","S","H",
                     "WSH","N","WSHN"),
                   each=length(plot_sample3)),
               rep("Lowest",times=6*length(plot_sample3)))

  
  DATA=rbind(DATA,y2,y3)
}

DATA=as.data.frame(DATA)
colnames(DATA)= c("y","x","Intervention","Completeness")
DATA$x=as.numeric(DATA$x)
DATA$y=as.numeric(DATA$y)
DATA$Intervention=factor(DATA$Intervention)

DATA_W = DATA[DATA$Intervention =="W",]
write.csv(DATA_W,"Table_Fig5A.csv")
p1 = ggplot(DATA_W, aes(x=factor(x),y=y,fill=Completeness))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1,alpha=0.75)+
  stat_summary(fun="median",geom="line",linewidth=1.25,aes(group=Completeness),position = position_dodge(width=0.9))+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=Completeness,bg=Completeness),position = position_dodge(width=0.9),show.legend = FALSE)+
  scale_fill_manual(values=c("#1dc956","#0a431d"))+
    guides(fill=guide_legend(title="Intervenable fraction quintile"))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  xlab("Coverage fraction")+
  ggtitle("Water (W) arm")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(0,1))+
  theme_classic(base_size = 8)+
  theme(legend.position = "bottom")+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=6))
p1
median(DATA_W[DATA_W$x==50 & DATA_W$Completeness=="Highest","y"])
median(DATA_W[DATA_W$x==50 & DATA_W$Completeness=="Lowest","y"])


DATA_S = DATA[DATA$Intervention =="S",]
write.csv(DATA_S,"Table_Fig5B.csv")
p2 = ggplot(DATA_S, aes(x=factor(x),y=y,fill=Completeness))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1,alpha=0.75)+
  stat_summary(fun="median",geom="line",linewidth=1.25,aes(group=Completeness),position = position_dodge(width=0.9))+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=Completeness,bg=Completeness),position = position_dodge(width=0.9),show.legend = FALSE)+
  scale_fill_manual(values=c("#7fccbf","#2c6d62"))+
  guides(fill=guide_legend(title="Intervenable fraction quintile"))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  xlab("Coverage fraction")+
  ggtitle("Sanitation (S) arm")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(0,1))+
  theme_classic(base_size = 8)+
  theme(legend.position = "bottom")+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=6))
p2

DATA_H = DATA[DATA$Intervention =="H",]
write.csv(DATA_H,"Table_Fig5C.csv")
p3 = ggplot(DATA_H, aes(x=factor(x),y=y,fill=Completeness))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1,alpha=0.75)+
  stat_summary(fun="median",geom="line",linewidth=1.25,aes(group=Completeness),position = position_dodge(width=0.9))+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=Completeness,bg=Completeness),position = position_dodge(width=0.9),show.legend = FALSE)+
  scale_fill_manual(values=c("#d2ecf9","#36a9e2"))+
  guides(fill=guide_legend(title="Intervenable fraction quintile"))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  xlab("Coverage fraction")+
  ggtitle("Hygiene (H) arm")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(0,1))+
  theme_classic(base_size = 8)+
  theme(legend.position = "bottom")+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=6))
p3
median(DATA_H[DATA_H$x==50 & DATA_H$Completeness=="Highest","y"])
median(DATA_H[DATA_H$x==50 & DATA_H$Completeness=="Lowest","y"])


DATA_WSH = DATA[DATA$Intervention =="WSH",]
write.csv(DATA_WSH,"Table_Fig5D.csv")
p4 = ggplot(DATA_WSH, aes(x=factor(x),y=y,fill=Completeness))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1,alpha=0.75)+
  stat_summary(fun="median",geom="line",linewidth=1.25,aes(group=Completeness),position = position_dodge(width=0.9))+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=Completeness,bg=Completeness),position = position_dodge(width=0.9),show.legend = FALSE)+
  scale_fill_manual(values=c("#f0e8c2","#ccb333"))+
  guides(fill=guide_legend(title="Intervenable fraction quintile"))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  xlab("Coverage fraction")+
  ggtitle("WSH arm")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(0,1))+
  theme_classic(base_size = 8)+
  theme(legend.position = "bottom")+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=6))
p4

DATA_N = DATA[DATA$Intervention =="N",]
write.csv(DATA_N,"Table_Fig5E.csv")
p5 = ggplot(DATA_N, aes(x=factor(x),y=y,fill=Completeness))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1,alpha=0.75)+
  stat_summary(fun="median",geom="line",linewidth=1.25,aes(group=Completeness),position = position_dodge(width=0.9))+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=Completeness,bg=Completeness),position = position_dodge(width=0.9),show.legend = FALSE)+
  scale_fill_manual(values=c("#e6b3bb","#993344"))+
  guides(fill=guide_legend(title="Intervenable fraction quintile"))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  xlab("Coverage fraction")+
  ggtitle("Nutrition (N) arm")+
  geom_hline(yintercept=0)+
  coord_cartesian(ylim=c(0,1))+
  theme_classic(base_size = 8)+
  theme(legend.position = "bottom")+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=6))
p5

DATA_WSHN = DATA[DATA$Intervention =="WSHN",]
write.csv(DATA_WSHN,"Table_Fig5F.csv")
p6 = ggplot(DATA_WSHN, aes(x=factor(x),y=y,fill=Completeness))+
  geom_violin(adjust=2.5,scale="width",color="black",lwd=0.1,alpha=0.75)+
  stat_summary(fun="median",geom="line",linewidth=1.25,aes(group=Completeness),position = position_dodge(width=0.9))+
  stat_summary(fun="median",geom="point",size=0.5,pch=21,aes(group=Completeness,bg=Completeness),position = position_dodge(width=0.9),show.legend = FALSE)+
  scale_fill_manual(values=c("#cc7fbf","#6d2c62"))+
  guides(fill=guide_legend(title="Intervenable fraction quintile"))+
  ylab("Change in intervention\neffectiveness (percentage points)")+
  xlab("Coverage fraction")+
  coord_cartesian(ylim=c(0,1))+
  ggtitle("WSHN arm")+
  geom_hline(yintercept=0)+
  theme_classic(base_size = 8)+
  theme(legend.position = "bottom")+
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 6),
        legend.title=element_text(size=7), 
        legend.text=element_text(size=6))
p6

ggarrange(p1,p2,p3,p4,p5,p6,labels="auto",label.x=0.9,nrow=3, ncol=2)#,common.legend = TRUE,legend="bottom")
ggsave("Figure 5.pdf",height=6, width=6)
