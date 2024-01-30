library(SobolSequence)
library(tictoc)
library(deSolve)
library(parallel)
library(ggplot2)
library(ggpubr)
library(nleqslv)


data = read.csv("washmodeldata_bangladesh.csv")
data=data[!is.na(data$ageyrs) & data$ageyrs>0,] #remove children with non-positive ages
data=data[!is.na(data$diar7d),] #remove those entries without diarrhea info


#Armid
# 1: Control
# 2: Handwashing
# 3: Nutrition
# 4: Nutrition + WSH
# 5: Sanitation
# 6: Water
# 7: WSH

weights = rep(1,nrow(data))
data=cbind(data,weights)


#Some assumptions
data[data$svy==0,"freechl"]=0 #No chlorine at baseline 
data[data$armid%in%c(1,2,5,6,7),"rlnsp"]=0 #no nutrition in non-N arms
data[is.na(data$rlnsp) & !is.na(data$lnsp),"rlnsp"]=data[is.na(data$rlnsp) & !is.na(data$lnsp),"lnsp"] #Use observed rather than reported if reported not available.
data[data$armid%in%c(3,4) & grepl("C",data$individualid,fixed=T),"rlnsp"]=0 #no nutrition for siblings

#Require full adherence data
data=data[!is.na(data$freechl),]
data=data[!is.na(data$hwsws),]
data=data[!is.na(data$latseal),]
data=data[!is.na(data$rlnsp),]

#Define nutrition adherence
data$nutrition = array(NA,nrow(data))
data[data$rlnsp >= 0.5,"nutrition"] = 1
data[data$rlnsp < 0.5,"nutrition"] = 0

#Define adherence groups
data$adh_cat= array(length(data[1,]))
data[data$nutrition ==0 & data$freechl ==0 & data$hwsws == 0 & data$latseal ==0,"adh_cat"] = "C"
data[data$nutrition ==1 & data$freechl ==0 & data$hwsws == 0 & data$latseal ==0,"adh_cat"] = "N"
data[data$nutrition ==0 & data$freechl ==1 & data$hwsws == 0 & data$latseal ==0,"adh_cat"] = "W"
data[data$nutrition ==0 & data$freechl ==0 & data$hwsws == 1 & data$latseal ==0,"adh_cat"] = "H"
data[data$nutrition ==0 & data$freechl ==0 & data$hwsws == 0 & data$latseal ==1,"adh_cat"] = "S"
data[data$nutrition ==1 & data$freechl ==1 & data$hwsws == 0 & data$latseal ==0,"adh_cat"] = "NW"
data[data$nutrition ==1 & data$freechl ==0 & data$hwsws == 1 & data$latseal ==0,"adh_cat"] = "NH"
data[data$nutrition ==1 & data$freechl ==0 & data$hwsws == 0 & data$latseal ==1,"adh_cat"] = "NS"
data[data$nutrition ==0 & data$freechl ==1 & data$hwsws == 1 & data$latseal ==0,"adh_cat"] = "WH"
data[data$nutrition ==0 & data$freechl ==1 & data$hwsws == 0 & data$latseal ==1,"adh_cat"] = "WS"
data[data$nutrition ==0 & data$freechl ==0 & data$hwsws == 1 & data$latseal ==1,"adh_cat"] = "HS"
data[data$nutrition ==1 & data$freechl ==1 & data$hwsws == 1 & data$latseal ==0,"adh_cat"] = "NWH"
data[data$nutrition ==1 & data$freechl ==1 & data$hwsws == 0 & data$latseal ==1,"adh_cat"] = "NWS"
data[data$nutrition ==1 & data$freechl ==0 & data$hwsws == 1 & data$latseal ==1,"adh_cat"] = "NHS"
data[data$nutrition ==0 & data$freechl ==1 & data$hwsws == 1 & data$latseal ==1,"adh_cat"] = "WHS"
data[data$nutrition ==1 & data$freechl ==1 & data$hwsws == 1 & data$latseal ==1,"adh_cat"] = "NWHS"

#What is the distribution of states for people not in the intervention?
#Assume it is the distribution at baseline
adh_cat_groups = c("C","N","W","H","S","NW","NH","NS","WH","WS","HS","NWH","NWS","NHS","WHS","NWHS")
baseline_adherence_vec = array(0,16)
baseline_adherence_vec[match(names(table(data[data$svy==0,"adh_cat"])),adh_cat_groups)] = table(data[data$svy==0,"adh_cat"])/sum(table(data[data$svy==0,"adh_cat"]))

########################################
#Adherence distributions

W_prev = matrix(NA,7,2)
for (i in 1:7){
  W_prev[i,1] = sum(data[data$armid==i & data$svy==0,"freechl"])/length(data[data$armid==i & data$svy==0,"freechl"])
  W_prev[i,2] = sum(data[data$armid==i & data$svy%in%c(1,2),"freechl"])/length(data[data$armid==i & data$svy%in%c(1,2),"freechl"])
}

S_prev = matrix(NA,7,2)
for (i in 1:7){
  S_prev[i,1] = sum(data[data$armid==i & data$svy==0,"latseal"])/length(data[data$armid==i & data$svy==0,"latseal"])
  S_prev[i,2] = sum(data[data$armid==i & data$svy%in%c(1,2),"latseal"])/length(data[data$armid==i & data$svy%in%c(1,2),"latseal"])
}

H_prev = matrix(NA,7,2)
for (i in 1:7){
  H_prev[i,1] = sum(data[data$armid==i & data$svy==0,"hwsws"])/length(data[data$armid==i & data$svy==0,"hwsws"])
  H_prev[i,2] = sum(data[data$armid==i & data$svy%in%c(1,2),"hwsws"])/length(data[data$armid==i & data$svy%in%c(1,2),"hwsws"])
}

N_prev = matrix(NA,7,2)
for (i in 1:7){
  N_prev[i,1] = sum(data[data$armid==i & data$svy==0,"nutrition"])/length(data[data$armid==i & data$svy==0,"nutrition"])
  N_prev[i,2] = sum(data[data$armid==i & data$svy%in%c(1,2),"nutrition"])/length(data[data$armid==i & data$svy%in%c(1,2),"nutrition"])
}

prev = c(W_prev[c(1,6,5,2,7,3,4),1],W_prev[c(1,6,5,2,7,3,4),2])
Arm = rep(c("C","W","S","H","WSH","N","WSH-N"),2)
Data = rep(c("Baseline","Midline/endline"),each=7)
DATAW=data.frame(prev,Arm,Data)
DATAW$Arm=factor(DATAW$Arm,levels=c("C","W","S","H","WSH","N","WSH-N"))
DATAW$Data=factor(DATAW$Data,levels=c("Baseline","Midline/endline"))
DW = ggplot( )+
  geom_bar(data = DATAW, aes(x=Arm, weight=prev, fill=Data),position = position_dodge(width = 0.75))+
  scale_fill_manual(guide="none",values=c("firebrick2","dodgerblue1"))+
  labs(fill="")+
  coord_cartesian(ylim=c(0,1))+
  ylab("Fraction using\nchlorination")+
  theme_classic(base_size=10)

prev = c(S_prev[c(1,6,5,2,7,3,4),1],S_prev[c(1,6,5,2,7,3,4),2])
Arm = rep(c("C","W","S","H","WSH","N","WSH-N"),2)
Data = rep(c("Baseline","Midline/endline"),each=7)
DATAS=data.frame(prev,Arm,Data)
DATAS$Arm=factor(DATAS$Arm,levels=c("C","W","S","H","WSH","N","WSH-N"))
DATAS$Data=factor(DATAS$Data,levels=c("Baseline","Midline/endline"))
DS = ggplot( )+
  geom_bar(data = DATAS, aes(x=Arm, weight=prev, fill=Data),position = position_dodge(width = 0.75))+
  scale_fill_manual(guide="none",values=c("firebrick2","dodgerblue1"))+
  labs(fill="")+
  coord_cartesian(ylim=c(0,1))+
  ylab("Fraction using latrine\nwith water seal")+
  theme_classic(base_size=10)

prev = c(H_prev[c(1,6,5,2,7,3,4),1],H_prev[c(1,6,5,2,7,3,4),2])
Arm = rep(c("C","W","S","H","WSH","N","WSH-N"),2)
Data = rep(c("Baseline","Midline/endline"),each=7)
DATAH=data.frame(prev,Arm,Data)
DATAH$Arm=factor(DATAH$Arm,levels=c("C","W","S","H","WSH","N","WSH-N"))
DATAH$Data=factor(DATAH$Data,levels=c("Baseline","Midline/endline"))
DH = ggplot( )+
  geom_bar(data = DATAH, aes(x=Arm, weight=prev, fill=Data),position = position_dodge(width = 0.75))+
  scale_fill_manual(guide="none", values=c("firebrick2","dodgerblue1"))+
  labs(fill="")+
  coord_cartesian(ylim=c(0,1))+
  ylab("Fraction using handwashing\nwith soap & water")+
  theme_classic(base_size=10)

prev = c(N_prev[c(1,6,5,2,7,3,4),1],N_prev[c(1,6,5,2,7,3,4),2])
Arm = rep(c("C","W","S","H","WSH","N","WSH-N"),2)
Data = rep(c("Baseline","Midline/endline"),each=7)
DATAN=data.frame(prev,Arm,Data)
DATAN$Arm=factor(DATAN$Arm,levels=c("C","W","S","H","WSH","N","WSH-N"))
DATAN$Data=factor(DATAN$Data,levels=c("Baseline","Midline/endline"))
DN = ggplot( )+
  geom_bar(data = DATAN, aes(x=Arm, weight=prev, fill=Data),position = position_dodge(width = 0.75))+
  scale_fill_manual(values=c("firebrick2","dodgerblue1"))+
  labs(fill="")+
  coord_cartesian(ylim=c(0,1))+
  ylab("Fraction consuming\nnutrition sachets")+
  theme_classic(base_size=10)+
  theme(legend.position = c(0.25,0.4),legend.direction = "vertical",
        legend.text = element_text(size=8))

Fig4=ggarrange(DW, DS, DH, DN, ncol=2, nrow=2,hjust=c(-7.5,-7,-7.5,-7), labels="auto",heights=c(1,1,1,1))
ggsave(file="Fig4.pdf", width = 7, height = 3.5)



#Compute cluster-level averages for each intervention
clusters=unique(data$clusterid)
n_clusters = length(clusters)

for (j in 0:2){
  for (i in 1:length(clusters)){
    
    data[data$svy==j & data$clusterid == clusters[i],"cluster_nutrition"]=mean(data[data$svy==j & data$clusterid == clusters[i],"nutrition"],na.rm=T)
    data[data$svy==j & data$clusterid == clusters[i],"cluster_freechl"]=mean(data[data$svy==j & data$clusterid == clusters[i],"freechl"],na.rm=T)
    data[data$svy==j & data$clusterid == clusters[i],"cluster_hwsws"]=mean(data[data$svy==j & data$clusterid == clusters[i],"hwsws"],na.rm=T)
    data[data$svy==j & data$clusterid == clusters[i],"cluster_latseal"]=mean(data[data$svy==j & data$clusterid == clusters[i],"latseal"],na.rm=T)
    
  }
}

#List all clusters in a way for fast evaluation
cluster_list = vector(mode = "list", length = n_clusters)
for (j in 0:2){
  for (i in 1:n_clusters){
    cluster_list[[(n_clusters*j+i)]] = list(data[data$svy==j & data$clusterid == clusters[i],],j,clusters[i])
  }
}

#For faster evaluation
data$model_prev=NA
data$NLL=NA

###################################################################################
#Model

model = function(x, model_par,N_vec){
  
  x=abs(x) #transformation keeps things positive
  x[N_vec==0]=0
  
  piN = model_par[1]
  piW = model_par[2]
  piH = model_par[3]
  piS = model_par[4]
  R0W = model_par[5]
  R0H = model_par[6]
  R0O = model_par[7]
  
  phiN = 1 - piN
  phiW = 1 - piW
  phiH = 1 - piH
  phiS = 1 - piS
  
  I = x[1]
  I_N = x[2]
  I_W= x[3]
  I_H= x[4]
  I_S= x[5]
  I_NW= x[6]
  I_NH= x[7]
  I_NS= x[8]
  I_WH= x[9]
  I_WS= x[10]
  I_HS= x[11]
  I_NWH= x[12]
  I_NWS= x[13]
  I_NHS = x[14]
  I_WHS= x[15]
  I_NWHS= x[16]
  
  S = N_vec[1] - x[1]
  S_N = N_vec[2] - x[2]
  S_W= N_vec[3] - x[3]
  S_H= N_vec[4] - x[4]
  S_S= N_vec[5] - x[5]
  S_NW= N_vec[6] - x[6]
  S_NH= N_vec[7] - x[7]
  S_NS= N_vec[8] - x[8]
  S_WH= N_vec[9] - x[9]
  S_WS= N_vec[10] - x[10]
  S_HS= N_vec[11] - x[11]
  S_NWH= N_vec[12] - x[12]
  S_NWS= N_vec[13] - x[13]
  S_NHS= N_vec[14] - x[14]
  S_WHS= N_vec[15] - x[15]
  S_NWHS= N_vec[16] - x[16]
  
  EW = phiS * sum(x[c(5,8,10,11,13:16)])  + sum(x[c(1:4,6:7,9,12)])
  EH = sum(x[1:16])
  EO = sum(x[1:16])
  
  y = numeric(length(x))
  y[1]  =  S * ( R0W * EW + R0H * EH + R0O * EO) -  I
  y[2]  =  S_N * phiN * ( R0W * EW  + R0H * EH + R0O * EO) -  I_N
  y[3]  =  S_W * ( R0W * EW * phiW + R0H * EH + R0O * EO) -  I_W
  y[4]  =  S_H * ( R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_H
  y[5]  =  S_S * ( R0W * EW + R0H * EH + R0O * EO) -  I_S
  y[6]  =  S_NW * phiN * ( R0W * EW * phiW + R0H * EH + R0O * EO) -  I_NW
  y[7]  =  S_NH * phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_NH
  y[8]  =  S_NS * phiN * ( R0W * EW + R0H * EH + R0O * EO) -  I_NS
  y[9]  =  S_WH * ( R0W * phiW * EW + R0H * EH * phiH  + R0O * EO) -  I_WH
  y[10] =  S_WS * ( R0W * phiW * EW + R0H * EH + R0O * EO) -  I_WS
  y[11] =  S_HS * ( R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_HS
  y[12] =  S_NWH * phiN * ( R0W * EW * phiW + R0H * EH * phiH  + R0O * EO) -  I_NWH
  y[13] =  S_NWS * phiN * ( R0W * EW * phiW + R0H * EH + R0O * EO) -  I_NWS
  y[14] =  S_NHS * phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_NHS
  y[15] =  S_WHS * ( R0W * EW * phiW + R0H  * EH * phiH  + R0O * EO) -  I_WHS
  y[16] =  S_NWHS * phiN * ( R0W * EW * phiW + R0H * EH * phiH  + R0O * EO) -  I_NWHS
  
  return(y)
  
}


model_Jac = function(x, model_par ,N_vec){
  
  x=abs(x)
  x[N_vec==0]=0

  piN = model_par[1]
  piW = model_par[2]
  piH = model_par[3]
  piS = model_par[4]
  R0W = model_par[5]
  R0H = model_par[6]
  R0O = model_par[7]
  
  phiN = 1 - piN
  phiW = 1 - piW
  phiH = 1 - piH
  phiS = 1 - piS
  
  I = x[1]
  I_N = x[2]
  I_W= x[3]
  I_H= x[4]
  I_S= x[5]
  I_NW= x[6]
  I_NH= x[7]
  I_NS= x[8]
  I_WH= x[9]
  I_WS= x[10]
  I_HS= x[11]
  I_NWH= x[12]
  I_NWS= x[13]
  I_NHS = x[14]
  I_WHS= x[15]
  I_NWHS= x[16]
  
  S = N_vec[1] - x[1]
  S_N = N_vec[2] - x[2]
  S_W= N_vec[3] - x[3]
  S_H= N_vec[4] - x[4]
  S_S= N_vec[5] - x[5]
  S_NW= N_vec[6] - x[6]
  S_NH= N_vec[7] - x[7]
  S_NS= N_vec[8] - x[8]
  S_WH= N_vec[9] - x[9]
  S_WS= N_vec[10] - x[10]
  S_HS= N_vec[11] - x[11]
  S_NWH= N_vec[12] - x[12]
  S_NWS= N_vec[13] - x[13]
  S_NHS= N_vec[14] - x[14]
  S_WHS= N_vec[15] - x[15]
  S_NWHS= N_vec[16] - x[16]
  
  EW = phiS * sum(x[c(5,8,10,11,13:16)])  + sum(x[c(1:4,6:7,9,12)])
  EH = sum(x[1:16])
  EO = sum(x[1:16])
  
  Df = matrix(0, length(x), length(x))
  
  # no intervention 
  # S * ( R0W * EW + R0H * EH + R0O * EO) -  I
  Df[1,c(1:4,6:7,9,12)]   = S * (R0W + R0H + R0O) #All groups without S. Note that this is wrong for 1=I, but we'll overwrite it -- this way we aren't messing around with indices too much
  Df[1,c(5,8,10,11,13:16)]= S * (R0W*phiS + R0H + R0O) #All groups with S.
  Df[1,1]   =  -(R0W * EW + R0H * EH + R0O * EO) + S * (R0W + R0H + R0O) - 1 #We correct 1=I here
  
  #N 
  #S_N * phiN * (R0W * EW  + R0H * EH + R0O * EO) -  I_N
  Df[2,c(1:4,6:7,9,12)]   = S_N * phiN * (R0W + R0H + R0O)
  Df[2,c(5,8,10,11,13:16)]= S_N * phiN * (R0W*phiS + R0H + R0O)
  Df[2,2]   = -phiN * (R0W * EW  +  R0H * EH + R0O * EO) + S_N * phiN * (R0W + R0H + R0O) -  1
  
  #W
  #S_W * ( R0W * EW * phiW + R0H * EH + R0O * EO) -  I_W
  Df[3,c(1:4,6:7,9,12)]   = S_W * (R0W * phiW + R0H + R0O)
  Df[3,c(5,8,10,11,13:16)]= S_W * (R0W * phiW * phiS + R0H + R0O)
  Df[3,3]   = -(R0W * EW * phiW + R0H * EH + R0O * EO) + S_W * (R0W * phiW + R0H + R0O) -  1
  
  #H
  # S_H * (R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_H
  Df[4,c(1:4,6:7,9,12)]   = S_H * (R0W + R0H * phiH  + R0O)
  Df[4,c(5,8,10,11,13:16)]= S_H * (R0W*phiS + R0H * phiH  + R0O)
  Df[4,4]   = -(R0W * EW + R0H * EH * phiH  + R0O * EO) + S_H * (R0W + R0H * phiH  + R0O) -  1
  
  #S
  # S_S * ( R0W * EW + R0H * EH + R0O * EO) -  I_S
  Df[5,c(1:4,6:7,9,12)]   = S_S * (R0W + R0H + R0O)
  Df[5,c(5,8,10,11,13:16)]= S_S * (R0W*phiS + R0H + R0O)
  Df[5,5]   = -(R0W * EW + R0H * EH + R0O * EO) + S_S * (R0W*phiS + R0H + R0O) -  1
  
  #NW
  # S_NW * phiN * ( R0W * EW * phiW + R0H * EH + R0O * EO) -  I_NW
  Df[6,c(1:4,6:7,9,12)]   = S_NW * phiN * (R0W * phiW + R0H + R0O)
  Df[6,c(5,8,10,11,13:16)]= S_NW * phiN * (R0W * phiW * phiS + R0H + R0O)
  Df[6,6]   = -phiN * (R0W * EW * phiW + R0H * EH + R0O * EO) + S_NW * phiN * (R0W * phiW + R0H + R0O) -  1
  
  #NH
  # S_NH * phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_NH
  Df[7,c(1:4,6:7,9,12)]   = S_NH * phiN * (R0W + R0H * phiH  + R0O)
  Df[7,c(5,8,10,11,13:16)]= S_NH * phiN * (R0W*phiS + R0H * phiH  + R0O)
  Df[7,7] = -phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) + S_NH * phiN * (R0W + R0H * phiH  + R0O) -  1
  
  #NS
  # S_NS * phiN * ( R0W * EW + R0H * EH + R0O * EO) -  I_NS
  Df[8,c(1:4,6:7,9,12)]   = S_NH * phiN * (R0W + R0H * phiH  + R0O)
  Df[8,c(5,8,10,11,13:16)]= S_NH * phiN * (R0W*phiS + R0H * phiH  + R0O)
  Df[8,8]   = -phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) + S_NH * phiN * (R0W + R0H * phiH  + R0O) -  1
  
  #WH
  # S_WH * ( R0W * phiW * EW + R0H * EH * phiH  + R0O * EO) -  I_WH
  Df[9,c(1:4,6:7,9,12)]   = S_WH * ( R0W * phiW + R0H * phiH  + R0O)
  Df[9,c(5,8,10,11,13:16)]= S_WH * ( R0W * phiW * phiS + R0H * phiH  + R0O)
  Df[9,9]   = -(R0W * phiW * EW + R0H * EH * phiH  + R0O * EO) + S_WH * ( R0W * phiW + R0H * phiH  + R0O) -  1
  
  #WS
  # S_WS * ( R0W * phiW * EW + R0H * EH + R0O * EO) -  I_WS
  Df[10,c(1:4,6:7,9,12)]   = S_WS * ( R0W * phiW + R0H + R0O)
  Df[10,c(5,8,10,11,13:16)]= S_WS * ( R0W * phiW * phiS + R0H + R0O)
  Df[10,10] = -(R0W * phiW * EW + R0H * EH + R0O * EO) + S_WS * ( R0W * phiW * phiS + R0H + R0O) -  1
  
  #HS
  # S_HS * ( R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_HS
  Df[11,c(1:4,6:7,9,12)]   = S_HS * ( R0W + R0H * phiH  + R0O)
  Df[11,c(5,8,10,11,13:16)]= S_HS * ( R0W * phiS + R0H * phiH  + R0O)
  Df[11,11] = -(R0W * EW + R0H * EH * phiH  + R0O * EO) + S_HS * ( R0W * phiS + R0H * phiH  + R0O) -  1
  
  #NWH
  # S_NWH * phiN * ( R0W * EW * phiW + R0H * EH * phiH  + R0O * EO) -  I_NWH
  Df[12,c(1:4,6:7,9,12)]   = S_NWH * phiN * (R0W * phiW + R0H * phiH  + R0O)
  Df[12,c(5,8,10,11,13:16)]= S_NWH * phiN * (R0W * phiW * phiS + R0H * phiH  + R0O)
  Df[12,12] = -phiN * (R0W * EW * phiW + R0H * EH * phiH  + R0O * EO) + S_NWH * phiN * (R0W * phiW + R0H * phiH  + R0O) -  1
  
  #NWS
  # S_NWS * phiN * ( R0W * EW * phiW + R0H * EH + R0O * EO) -  I_NWS
  Df[13,c(1:4,6:7,9,12)]   = S_NWS * phiN * (R0W * phiW + R0H + R0O)
  Df[13,c(5,8,10,11,13:16)]= S_NWS * phiN * (R0W * phiW *phiS + R0H + R0O)
  Df[13,13] = -phiN * (R0W * EW * phiW + R0H * EH + R0O * EO) + S_NWS * phiN * (R0W * phiW *phiS + R0H + R0O) -  1
  
  #NHS
  # S_NHS * phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_NHS
  Df[14,c(1:4,6:7,9,12)]   = S_NHS * phiN * (R0W + R0H * phiH  + R0O)
  Df[14,c(5,8,10,11,13:16)]= S_NHS * phiN * (R0W * phiS + R0H * phiH  + R0O)
  Df[14,14] = - phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) + S_NHS * phiN * (R0W * phiS + R0H * phiH  + R0O) -  1
  
  #WHS
  # S_WHS * (R0W * EW * phiW + R0H  * EH * phiH  + R0O * EO) -  I_WHS
  Df[15,c(1:4,6:7,9,12)]   = S_WHS * (R0W * phiW + R0H * phiH  + R0O)
  Df[15,c(5,8,10,11,13:16)]= S_WHS * (R0W * phiW * phiS + R0H * phiH  + R0O)
  Df[15,15] = - phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) + S_WHS * (R0W * phiW * phiS + R0H * phiH  + R0O) -  1
  
  #NWHS
  # S_NWHS * phiN * (R0W * EW * phiW + R0H * EH * phiH  + R0O * EO) -  I_NWHS
  Df[16,c(1:4,6:7,9,12)]   = S_NWHS * phiN * (R0W * phiW + R0H * phiH  + R0O)
  Df[16,c(5,8,10,11,13:16)]= S_NWHS * phiN * (R0W * phiW * phiS + R0H * phiH  + R0O)
  Df[16,16] = -phiN * (R0W * EW * phiW + R0H * EH * phiH  + R0O * EO) + S_NWHS * phiN * (R0W * phiW * phiS + R0H * phiH  + R0O) -  1
  
  return(Df)
  
}

evaluate_cluster_NLL = function(k, par){
  par=abs(par)
  
  #R0 = par[1]
  #The pathway R0s are awkwardly parameterized, but it keeps p2 and p3 between
  #0 and 1, which is easier for sampling.
  R0W = par[1]*par[2]
  R0H = par[1]*(1-par[2])*par[3]
  R0O = par[1]*(1-par[2])*(1-par[3])
  piN = par[4]
  piW = par[5]
  piH = par[6]
  piS = par[7]
  omega = par[8]
  #Par 9 and 10 are the relative effectiveness of the H and S pre-existing conditions
  #compared to the intervention
  #Par 11 and 12 are the relative R0 of midline and endline vs baseline
  #Par 13-18 are the relative R0 of each arm relative to the control arm

  data_cluster_temp = cluster_list[[k]][[1]]
  j = cluster_list[[k]][[2]]
  
  if (nrow(data_cluster_temp)==0){return(0)} else{ #skip the rest if there is no eligible people in this cluster for this survey
    
    #Set arm and svy specific R0 adjustment
    #Armid
    # 1: Control
    # 2: Handwashing
    # 3: Nutrition
    # 4: Nutrition + WSH
    # 5: Sanitation
    # 6: Water
    # 7: WSH
    
    svy_temp = j+1
    arm_temp = data_cluster_temp$armid[1]
    R0_adj = c(1,par[11:12])[svy_temp]*c(1,par[13:18])[arm_temp]
    
    piH_adj=1
    piS_adj=1
    if (j == 0 | is.element(arm_temp,c(1,3,5,6))){
      piH_adj =par[9]
    }
    if (j == 0 | is.element(arm_temp,c(1,2,3,6))){
      piS_adj = par[10]
    }
    
    model_par= c(piN, piW, piH*piH_adj, piS*piS_adj, R0W*R0_adj, R0H*R0_adj, R0O*R0_adj)
    
    #Distribution of WSH groups
    rho_vec = c(sum(data_cluster_temp$adh_cat == "C"),
                sum(data_cluster_temp$adh_cat == "N"),
                sum(data_cluster_temp$adh_cat == "W"),
                sum(data_cluster_temp$adh_cat == "H"),
                sum(data_cluster_temp$adh_cat == "S"),
                sum(data_cluster_temp$adh_cat == "NW"),
                sum(data_cluster_temp$adh_cat == "NH"),
                sum(data_cluster_temp$adh_cat == "NS"),
                sum(data_cluster_temp$adh_cat == "WH"),
                sum(data_cluster_temp$adh_cat == "WS"),
                sum(data_cluster_temp$adh_cat == "HS"),
                sum(data_cluster_temp$adh_cat == "NWH"),
                sum(data_cluster_temp$adh_cat == "NWS"),
                sum(data_cluster_temp$adh_cat == "NHS"),
                sum(data_cluster_temp$adh_cat == "WHS"),
                sum(data_cluster_temp$adh_cat == "NWHS"))/length(data_cluster_temp$adh_cat)
    
    prev= max(1-1/sum(model_par[5:7]),1E-10) #initial condition prevalence 
    N_vec= rho_vec*omega + baseline_adherence_vec*(1-omega)
    I0 = prev*N_vec

    steady_state = abs(nleqslv(x=I0,fn=model,jac = model_Jac, 
                               model_par=model_par, N_vec = N_vec,
                               method="Broyden", control=list(allowSingular=1))$x)
    steady_state[N_vec=0]=0
    
    #Calculate prevalence within each group
    steady_state_normalized = steady_state/N_vec

    data_cluster_temp[data_cluster_temp$adh_cat == "C","model_prev"] = steady_state_normalized[1]
    data_cluster_temp[data_cluster_temp$adh_cat == "N","model_prev"] = steady_state_normalized[2]
    data_cluster_temp[data_cluster_temp$adh_cat == "W","model_prev"] = steady_state_normalized[3]
    data_cluster_temp[data_cluster_temp$adh_cat == "H","model_prev"] = steady_state_normalized[4]
    data_cluster_temp[data_cluster_temp$adh_cat == "S","model_prev"] = steady_state_normalized[5]
    data_cluster_temp[data_cluster_temp$adh_cat == "NW","model_prev"] = steady_state_normalized[6]
    data_cluster_temp[data_cluster_temp$adh_cat == "NH","model_prev"] = steady_state_normalized[7]
    data_cluster_temp[data_cluster_temp$adh_cat == "NS","model_prev"] = steady_state_normalized[8]
    data_cluster_temp[data_cluster_temp$adh_cat == "WH","model_prev"] = steady_state_normalized[9]
    data_cluster_temp[data_cluster_temp$adh_cat == "WS","model_prev"] = steady_state_normalized[10]
    data_cluster_temp[data_cluster_temp$adh_cat == "HS","model_prev"] = steady_state_normalized[11]
    data_cluster_temp[data_cluster_temp$adh_cat == "NWH","model_prev"] = steady_state_normalized[12]
    data_cluster_temp[data_cluster_temp$adh_cat == "NWS","model_prev"] = steady_state_normalized[13]
    data_cluster_temp[data_cluster_temp$adh_cat == "NHS","model_prev"] = steady_state_normalized[14]
    data_cluster_temp[data_cluster_temp$adh_cat == "WHS","model_prev"] = steady_state_normalized[15]
    data_cluster_temp[data_cluster_temp$adh_cat == "NWHS","model_prev"] = steady_state_normalized[16]
    
    data_cluster_temp[,"NLL"] = -data_cluster_temp$weights * (data_cluster_temp$diar7d * log(data_cluster_temp$model_prev) + (1-data_cluster_temp$diar7d)* log(1-data_cluster_temp$model_prev))
    
    NLL = sum(data_cluster_temp[,"NLL"],na.rm=T)
    
    return(NLL)
  }
}

fullLogLikelihood = function(par, list){
  
  NLL_vec= sapply(1:length(list), evaluate_cluster_NLL, par= par)
  
  return(sum(NLL_vec))
  
}

fullLogLikelihood_wrapper = function(logit2_par_R0, par_other,list){
  
  logit2_par_R0 =  pmin(logit2_par_R0,rep(700,9)) #To prevent NAN
  par_R0=exp(logit2_par_R0)/(1+exp(logit2_par_R0))
  #Transform R0 parameters so that R0 is between 1 and 2 and the adj parameters 
  #are between 0 and 2.
  NLL = fullLogLikelihood(c((1+par_R0[1]), par_other,2*par_R0[2:9]),list)
  
  return(NLL)
  
}

compute_NLL_function=function(x){
  
  par_other = scaled_sample[x,]
  logit2_par_R0_init = c(-1, rep(0,8))
  
  out = nlm(f=fullLogLikelihood_wrapper, p=logit2_par_R0_init, par_other=par_other, list=cluster_list,
            steptol=1e-4,gradtol=1e-4)
  
  par_R0 = exp(out$estimate)/(1+exp(out$estimate))
  
  print(x)
  
  return( c((1+par_R0[1]),par_other,2*par_R0[2:9], out$minimum))
}

###################################################################################
set.seed(0)
nsample = 1E3*50
sobol_sample = sobolSequence.points(9,count=nsample)
paramlow =  c( 0,  0,  0,    0,     0,   0, 1E-5,0,0)
paramhigh = c(1,  1,    1,  1,     1,  1,  0.2,1,1)
scaled_sample = t(paramlow + (paramhigh-paramlow)*t(sobol_sample))

n_cores = 120
print("starting cluster")
cluster = makeCluster(n_cores)
clusterExport(cluster,"cluster_list")
clusterExport(cluster,"scaled_sample")
clusterExport(cluster,"fullLogLikelihood")
clusterExport(cluster,"fullLogLikelihood_wrapper")
clusterExport(cluster,"evaluate_cluster_NLL")
clusterExport(cluster,"baseline_adherence_vec")
clusterExport(cluster,"model")
clusterExport(cluster,"model_Jac")
clusterEvalQ(cluster,library("deSolve"))
clusterEvalQ(cluster,library("nleqslv"))
print("cluster loaded")

tic()
sample_and_NLL = parSapply(cluster,1:50000,FUN=compute_NLL_function)
toc()
saveRDS(sample_and_NLL,"sample_and_NLL_iterates.RDS")

stopCluster(cluster)
print("stopping cluster")

# ###################################################################################
# #Sampling-importance resampling

sample_and_NLL = as.data.frame(t(readRDS("sample_and_NLL_iterates.RDS")))

names(sample_and_NLL)[19]="NLL"
sample_and_NLL$rNLL = sample_and_NLL$NLL- min(sample_and_NLL$NLL)
sample_and_NLL$prob = (exp(-sample_and_NLL$rNLL))
saveRDS(sample_and_NLL,"sample_and_NLL.RDS")
sample_and_NLL = readRDS("sample_and_NLL.RDS")

set.seed(0)
resample = sample(1:nrow(sample_and_NLL), 5.0E4, replace=TRUE, prob=sample_and_NLL$prob)
saveRDS(resample,"resample.RDS")
resample = readRDS("resample.RDS")
index= unique(resample)

# ###################################################################################
# #Calculate prevalence in each arm at each svy for the given parameters

prevalence_by_arm = function(par, dataset){
  par=abs(par)
  R0W = par[1]*par[2]
  R0H = par[1]*(1-par[2])*par[3]
  R0O = par[1]*(1-par[2])*(1-par[3])
  piN = par[4]
  piW = par[5]
  piH = par[6]
  piS = par[7]
  omega = par[8]

  clusters=unique(dataset$clusterid)
  for (j in 0:2){
    for (i in 1:length(clusters)){

      data_cluster_temp = dataset[dataset$svy==j & dataset$clusterid == clusters[i],]

      if (nrow(data_cluster_temp)==0){next} #skip the rest if there is no eligible people in this cluster for this survey

      #Set arm and svy specific R0 adjustment
      svy_temp = j+1
      arm_temp = data_cluster_temp$armid[1]
      R0_adj = c(1,par[11:12])[svy_temp]*c(1,par[13:18])[arm_temp]

      piH_adj=1
      piS_adj=1
      if (j == 0 | is.element(arm_temp,c(1,3,5,6))){
        piH_adj =par[9]
      }
      if (j == 0 | is.element(arm_temp,c(1,2,3,6))){
        piS_adj = par[10]
      }

      model_par= c(piN, piW, piH*piH_adj, piS*piS_adj, R0W*R0_adj, R0H*R0_adj, R0O*R0_adj)

      #Distribution of WSH groups
      rho_vec = c(sum(data_cluster_temp$adh_cat == "C"),
                  sum(data_cluster_temp$adh_cat == "N"),
                  sum(data_cluster_temp$adh_cat == "W"),
                  sum(data_cluster_temp$adh_cat == "H"),
                  sum(data_cluster_temp$adh_cat == "S"),
                  sum(data_cluster_temp$adh_cat == "NW"),
                  sum(data_cluster_temp$adh_cat == "NH"),
                  sum(data_cluster_temp$adh_cat == "NS"),
                  sum(data_cluster_temp$adh_cat == "WH"),
                  sum(data_cluster_temp$adh_cat == "WS"),
                  sum(data_cluster_temp$adh_cat == "HS"),
                  sum(data_cluster_temp$adh_cat == "NWH"),
                  sum(data_cluster_temp$adh_cat == "NWS"),
                  sum(data_cluster_temp$adh_cat == "NHS"),
                  sum(data_cluster_temp$adh_cat == "WHS"),
                  sum(data_cluster_temp$adh_cat == "NWHS"))/length(data_cluster_temp$adh_cat)

      prev= max(1-1/sum(model_par[5:7]),1E-10) #initial condition prevalence 
      N_vec= rho_vec*omega + baseline_adherence_vec*(1-omega)
      I0 = prev*N_vec
      
      steady_state = abs(nleqslv(x=I0,fn=model,jac = model_Jac, 
                                 model_par=model_par, N_vec = N_vec,
                                 method="Broyden", control=list(allowSingular=1))$x)
      steady_state[N_vec=0]=0
      
      #Calculate prevalence within each group
      steady_state_normalized = steady_state/N_vec

      data_cluster_temp[data_cluster_temp$adh_cat == "C","model_prev"] = steady_state_normalized[1]
      data_cluster_temp[data_cluster_temp$adh_cat == "N","model_prev"] = steady_state_normalized[2]
      data_cluster_temp[data_cluster_temp$adh_cat == "W","model_prev"] = steady_state_normalized[3]
      data_cluster_temp[data_cluster_temp$adh_cat == "H","model_prev"] = steady_state_normalized[4]
      data_cluster_temp[data_cluster_temp$adh_cat == "S","model_prev"] = steady_state_normalized[5]
      data_cluster_temp[data_cluster_temp$adh_cat == "NW","model_prev"] = steady_state_normalized[6]
      data_cluster_temp[data_cluster_temp$adh_cat == "NH","model_prev"] = steady_state_normalized[7]
      data_cluster_temp[data_cluster_temp$adh_cat == "NS","model_prev"] = steady_state_normalized[8]
      data_cluster_temp[data_cluster_temp$adh_cat == "WH","model_prev"] = steady_state_normalized[9]
      data_cluster_temp[data_cluster_temp$adh_cat == "WS","model_prev"] = steady_state_normalized[10]
      data_cluster_temp[data_cluster_temp$adh_cat == "HS","model_prev"] = steady_state_normalized[11]
      data_cluster_temp[data_cluster_temp$adh_cat == "NWH","model_prev"] = steady_state_normalized[12]
      data_cluster_temp[data_cluster_temp$adh_cat == "NWS","model_prev"] = steady_state_normalized[13]
      data_cluster_temp[data_cluster_temp$adh_cat == "NHS","model_prev"] = steady_state_normalized[14]
      data_cluster_temp[data_cluster_temp$adh_cat == "WHS","model_prev"] = steady_state_normalized[15]
      data_cluster_temp[data_cluster_temp$adh_cat == "NWHS","model_prev"] = steady_state_normalized[16]

      dataset[dataset$svy==j & dataset$clusterid == clusters[i],"model_prev"] = data_cluster_temp[,"model_prev"]

    }
  }

  temp_prev_model = as.data.frame(matrix(NA,7,4))
  for (i in 1:7){
    for (j in 0:2){
      temp = dataset[dataset$svy==j & dataset$armid == i,]
      prev = sum(temp$model_prev) / nrow(temp)
      temp_prev_model[i,(j+1)] = prev
    }
    temp = dataset[dataset$svy%in%c(1,2) & dataset$armid == i,]
    prev = sum(temp$model_prev) / nrow(temp)
    temp_prev_model[i,4] = prev
  }
  return(cbind(temp_prev_model))
}

counterfactual_simulation = function(i){
  print(which(unique(index)==i))
  prev_mat = prevalence_by_arm(par=unlist(sample_and_NLL[i,1:18]), dataset=data)

 #By wave, prevalence at baseline, midline, endline, midline/endline average
  prevalences = c(prev_mat[1,1],prev_mat[1,2],prev_mat[1,3],prev_mat[1,4],
                  prev_mat[2,1],prev_mat[2,2],prev_mat[2,3],prev_mat[2,4],
                  prev_mat[3,1],prev_mat[3,2],prev_mat[3,3],prev_mat[3,4],
                  prev_mat[4,1],prev_mat[4,2],prev_mat[4,3],prev_mat[4,4],
                  prev_mat[5,1],prev_mat[5,2],prev_mat[5,3],prev_mat[5,4],
                  prev_mat[6,1],prev_mat[6,2],prev_mat[6,3],prev_mat[6,4],
                  prev_mat[7,1],prev_mat[7,2],prev_mat[7,3],prev_mat[7,4])
  return(prevalences)
}

index_prevalences = sapply(index,FUN=counterfactual_simulation)
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "prevalences.RDS")
prevalences = readRDS("prevalences.RDS")

#Extract arm-svy prevalences
#Baseline
prevC0    = prevalences[,1]
prevH0    = prevalences[,5]
prevN0    = prevalences[,9]
prevWSHN0 = prevalences[,13]
prevS0    = prevalences[,17]
prevW0    = prevalences[,21]
prevWSH0  = prevalences[,25]
#Midline/Endline average
prevC1    = prevalences[,2]
prevH1    = prevalences[,6]
prevN1    = prevalences[,10]
prevWSHN1 = prevalences[,14]
prevS1    = prevalences[,18]
prevW1    = prevalences[,22]
prevWSH1  = prevalences[,26]
#Midline/Endline average
prevC2    = prevalences[,3]
prevH2    = prevalences[,7]
prevN2    = prevalences[,11]
prevWSHN2 = prevalences[,15]
prevS2    = prevalences[,19]
prevW2    = prevalences[,23]
prevWSH2  = prevalences[,27]
#Midline/Endline average
prevC12    = prevalences[,4]
prevH12    = prevalences[,8]
prevN12    = prevalences[,12]
prevWSHN12 = prevalences[,16]
prevS12    = prevalences[,20]
prevW12    = prevalences[,24]
prevWSH12  = prevalences[,28]

obs_prev = matrix(NA,7,4)
for (i in 1:7){
  obs_prev[i,1] = sum(data[data$armid==i & data$svy==0,"diar7d"])/length(data[data$armid==i & data$svy==0,"diar7d"])
  obs_prev[i,2] = sum(data[data$armid==i & data$svy==1,"diar7d"])/length(data[data$armid==i & data$svy==1,"diar7d"])
  obs_prev[i,3] = sum(data[data$armid==i & data$svy==2,"diar7d"])/length(data[data$armid==i & data$svy==2,"diar7d"])
  obs_prev[i,4] = sum(data[data$armid==i & data$svy%in%c(1,2),"diar7d"])/length(data[data$armid==i & data$svy%in%c(1,2),"diar7d"])
}
obs_prev_se = matrix(NA,7,4)
for (i in 1:7){
  p1= sum(data[data$armid==i & data$svy==0,"diar7d"])/length(data[data$armid==i & data$svy==0,"diar7d"])
  p2= sum(data[data$armid==i & data$svy==1,"diar7d"])/length(data[data$armid==i & data$svy==1,"diar7d"])
  p3= sum(data[data$armid==i & data$svy==2,"diar7d"])/length(data[data$armid==i & data$svy==2,"diar7d"])
  p4= sum(data[data$armid==i & data$svy%in%c(1,2),"diar7d"])/length(data[data$armid==i & data$svy%in%c(1,2),"diar7d"])
  obs_prev_se[i,1] = sqrt(p1*(1-p1)/nrow(data[data$armid==i & data$svy==0,]))
  obs_prev_se[i,2] = sqrt(p2*(1-p2)/nrow(data[data$armid==i & data$svy==1,]))
  obs_prev_se[i,3] = sqrt(p3*(1-p3)/nrow(data[data$armid==i & data$svy==2,]))
  obs_prev_se[i,4] = sqrt(p4*(1-p4)/nrow(data[data$armid==i & data$svy%in%c(1,2),]))
}
obs_rel_risk = matrix(NA,7,3)
for (i in 1:7){
  p1= sum(data[data$armid==i & data$svy==0,"diar7d"])/length(data[data$armid==i & data$svy==0,"diar7d"])
  p2= sum(data[data$armid==i & data$svy%in%c(1,2),"diar7d"])/length(data[data$armid==i & data$svy%in%c(1,2),"diar7d"])
  obs_rel_risk[i,1] = p2/p1
  se = sqrt((1-p1)/(p1*nrow(data[data$armid==i & data$svy==0,]))+(1-p2)/(p2*nrow(data[data$armid==i & data$svy%in%c(1,2),])))
  z= qnorm(0.975,0,1)
  obs_rel_risk[i,2] = (p2/p1) * exp(- z*se)
  obs_rel_risk[i,3] = (p2/p1) * exp(+ z*se)
}

#Data prevalence
prev14 = c(obs_prev[c(1,6,5,2,7,3,4),1],obs_prev[c(1,6,5,2,7,3,4),4])
prev_low14 = c(obs_prev[c(1,6,5,2,7,3,4),1],obs_prev[c(1,6,5,2,7,3,4),4]) -qnorm(0.975,0,1)*c(obs_prev_se[c(1,6,5,2,7,3,4),1],obs_prev_se[c(1,6,5,2,7,3,4),4])
prev_high14 = c(obs_prev[c(1,6,5,2,7,3,4),1],obs_prev[c(1,6,5,2,7,3,4),4]) +qnorm(0.975,0,1)*c(obs_prev_se[c(1,6,5,2,7,3,4),1],obs_prev_se[c(1,6,5,2,7,3,4),4])
Arm14 = rep(c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"),2)
Data14 = rep(c("Baseline","Midline/endline"),each=7)
DATA1a=data.frame(prev14,prev_low14,prev_high14,Arm14,Data14)
DATA1a$Arm=factor(DATA1a$Arm,levels=c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"))
DATA1a$Data=factor(DATA1a$Data,levels=c("Baseline","Midline/endline"))

prev123 = c(obs_prev[c(1,6,5,2,7,3,4),1],obs_prev[c(1,6,5,2,7,3,4),2],obs_prev[c(1,6,5,2,7,3,4),3])
prev_low123 = c(obs_prev[c(1,6,5,2,7,3,4),1],obs_prev[c(1,6,5,2,7,3,4),2],obs_prev[c(1,6,5,2,7,3,4),3])-
  qnorm(0.975,0,1)*c(obs_prev_se[c(1,6,5,2,7,3,4),1],obs_prev_se[c(1,6,5,2,7,3,4),2],obs_prev_se[c(1,6,5,2,7,3,4),3])
prev_high123 = c(obs_prev[c(1,6,5,2,7,3,4),1],obs_prev[c(1,6,5,2,7,3,4),2],obs_prev[c(1,6,5,2,7,3,4),3])+
  qnorm(0.975,0,1)*c(obs_prev_se[c(1,6,5,2,7,3,4),1],obs_prev_se[c(1,6,5,2,7,3,4),2],obs_prev_se[c(1,6,5,2,7,3,4),3])
Arm = rep(c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"),3)
Data = rep(c("Baseline","Midline","Endline"),each=7)
DATA1b=data.frame(prev123,prev_low123,prev_high123,Arm,Data)
DATA1b$Arm=factor(DATA1b$Arm,levels=c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"))
DATA1b$Data=factor(DATA1b$Data,levels=c("Baseline","Midline", "Endline"))

#Data relative risk
prev = obs_rel_risk[c(1,6,5,2,7,3,4),1]
prev_low = obs_rel_risk[c(1,6,5,2,7,3,4),2]
prev_high = obs_rel_risk[c(1,6,5,2,7,3,4),3]
Arm = rep(c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"),1)
Data = rep(c("Midline/endline"),each=7)
DATA2=data.frame(prev,prev_low,prev_high,Arm,Data)
DATA2$Arm=factor(DATA2$Arm,levels=c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"))
DATA2$Data=factor(DATA2$Data,levels=c("Midline/endline"))

#Model prevalence
Prevalence = c(prevC0,  prevC12,
               prevW0,  prevW12,
               prevS0,  prevS12,
               prevH0,  prevH12,
               prevWSH0,  prevWSH12,
               prevN0,  prevN12,
               prevWSHN0,  prevWSHN12)
Arm = rep(c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"),each=2*nrow(prevalences))
Simulation = rep(rep(c("Baseline","Midline/endline"),each=nrow(prevalences)),7)
DATA3a=data.frame(Prevalence,Arm,Simulation)
DATA3a$Arm=factor(DATA3a$Arm,levels=c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"))
DATA3a$Simulation=factor(DATA3a$Simulation,levels=c("Baseline","Midline/endline"))

Prevalence = c(prevC0,  prevC1, prevC2,
               prevW0,  prevW1, prevW2,
               prevS0,  prevS1, prevS2,
               prevH0,  prevH1, prevH2,
               prevWSH0,  prevWSH1, prevWSH2,
               prevN0,  prevN1, prevN2,
               prevWSHN0,  prevWSHN1,  prevWSHN2)
Arm = rep(c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"),each=3*nrow(prevalences))
Simulation = rep(rep(c("Baseline","Midline","Endline"),each=nrow(prevalences)),7)
DATA3b=data.frame(Prevalence,Arm,Simulation)
DATA3b$Arm=factor(DATA3b$Arm,levels=c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"))
DATA3b$Simulation=factor(DATA3b$Simulation,levels=c("Baseline","Midline", "Endline"))


#Model relative risk
Rel_risk = c( prevC12/prevC0,
              prevW12/prevW0,
              prevS12/prevS0,
              prevH12/prevH0,
              prevWSH12/prevWSH0,
              prevN12/prevN0,
              prevWSHN12/prevWSHN0)
Arm = rep(c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"),each=nrow(prevalences))
Simulation = rep(rep(c("Midline/endline"),each=nrow(prevalences)),7)
DATA4=data.frame(Rel_risk,Arm,Simulation)
DATA4$Arm=factor(DATA4$Arm,levels=c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"))
DATA4$Simulation=factor(DATA4$Simulation,levels=c("Baseline","Midline/endline"))

#Data and model plots
#prevalence
Pplota= ggplot()+
  geom_violin(data = DATA3a, aes(x=Arm, y=Prevalence, fill=Simulation), adjust = 4,scale="count",width=2.5,position = position_dodge(width =1), color = NA)+
  scale_fill_manual(values=c("firebrick2","dodgerblue1"))+
  geom_pointrange(data = DATA1a, aes(x=Arm, y=prev14, ymin = prev_low14, ymax = prev_high14, color=Data),position = position_dodge(width = 0.5))+
  scale_color_manual(values=c("firebrick2","dodgerblue1"))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.12))+
  ylab("Prevalence")+
  theme_classic()
Pplota

Pplotb= ggplot()+
  geom_violin(data = DATA3b, aes(x=Arm, y=Prevalence, fill=Simulation), adjust = 4,scale="count",width=2.5,position = position_dodge(width =1), color = NA,alpha=0.75)+
  scale_fill_manual(values=c("firebrick2","skyblue","blue3"))+
  geom_pointrange(data = DATA1b, aes(x=Arm, y=prev123, ymin = prev_low123, ymax = prev_high123, color=Data,shape=Data),position = position_dodge(width = 0.5))+
  scale_shape_discrete()+
  scale_color_manual(values=c("firebrick2","skyblue","blue3"))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.12))+
  ylab("Prevalence")+
  theme_classic()
Pplotb
ggsave(file="FigS2.pdf", width = 7, height = 3.5)

#Relative risk
Rplot = ggplot()+
  geom_violin(data = DATA4, aes(x=Arm, y=Rel_risk, fill=Simulation), adjust = 4,scale="count",width=1, color = NA)+
  scale_fill_manual(values=c("dodgerblue1"))+
  geom_pointrange(data = DATA2, aes(x=Arm, y=prev, ymin = prev_low, ymax = prev_high, color=Data))+
  scale_color_manual(values=c("black"))+
  scale_y_log10(expand = c(0, 0),limits=c(0.25,1.25),breaks=c(0.25, 0.5,0.75, 1,1.25))+
  ylab("Prevalence relative to baseline (log-scale)")+
  theme_classic()+
  guides(fill = guide_legend(order =1),col = guide_legend(order = 2))
Rplot

Fig_fit=ggarrange(Pplota, Rplot, ncol=1, nrow=2,hjust=-6, labels="auto",heights=c(1,1))
ggsave(file="Fig5.pdf", width = 7, height = 7)

###################################################################################
#Parameter histograms
resample_sample = sample_and_NLL[resample,]

#NLL
p0 = ggplot() +
  geom_histogram(aes(x=resample_sample$NLL,y = ..density..),binwidth = 1, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=sample_and_NLL$NLL,y = ..density..),binwidth = 1, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Negative log-likelihood")+xlim(c(3240,3340))+
  scale_y_continuous(expand = expansion(mult=c(0, 0.1)))+
  annotate("text",x=3255, y=0.28, hjust=0, label = paste("Posterior"))+
  annotate("text",x=3260, y=0.05, hjust=0, label = paste("Prior"))+
  theme_classic()
ggsave(file="FigS1.pdf", width = 7, height = 3.5)

#R0 and pathway-specific R0s
nbins= 15
p1= ggplot() +
  geom_histogram(aes(x=resample_sample[,1],y = ..density..),binwidth = 1/nbins, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  geom_vline(xintercept = mean(resample_sample[,1]),col = "grey10")+
  ylab("Density") + xlab("Overall R0")+
  xlim(c(0,1.4))+
  scale_y_continuous(limits=c(0,14),expand = expansion(mult=c(0,0.1)),labels = scales::number_format(accuracy = 0.1))+
  theme_classic()

p2=ggplot() +
  geom_histogram(aes(x=sample_and_NLL[,1]*sample_and_NLL[,2],y = ..density..),binwidth = 1/nbins, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,1]*resample_sample[,2],y = ..density..),binwidth = 1/nbins, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  geom_vline(xintercept = median(resample_sample[,1]*resample_sample[,2]),col = "grey10")+
  ylab("Density") + xlab("Water R0")+
  xlim(c(0,1.4))+
  scale_y_continuous(limits=c(0,14),expand = expansion(mult=c(0,0.1)),labels = scales::number_format(accuracy = 0.1))+
  theme_classic()

p3=ggplot() +
  geom_histogram(aes(x=sample_and_NLL[,1]*(1-sample_and_NLL[,2])*sample_and_NLL[,3],y = ..density..),binwidth = 1/nbins, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,1]*(1-resample_sample[,2])*resample_sample[,3],y = ..density..),binwidth = 1/nbins, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  geom_vline(xintercept = median(resample_sample[,1]*(1-resample_sample[,2])*resample_sample[,3]),col = "grey10")+
  ylab("Density") + xlab("Fomite R0")+
  xlim(c(0,1.4))+
  scale_y_continuous(limits=c(0,14),expand = expansion(mult=c(0,0.1)),labels = scales::number_format(accuracy = 0.1))+
  theme_classic()

p4=ggplot() +
  geom_histogram(aes(x=sample_and_NLL[,1]*(1-sample_and_NLL[,2])*(1-sample_and_NLL[,3]),y = ..density..),binwidth = 1/nbins, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,1]*(1-resample_sample[,2])*(1-resample_sample[,3]),y = ..density..),binwidth = 1/nbins, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  geom_vline(xintercept = median(resample_sample[,1]*(1-resample_sample[,2])*(1-resample_sample[,3])),col = "grey10")+
  ylab("Density") + xlab("Other R0")+
  xlim(c(0,1.4))+
  scale_y_continuous(limits=c(0,14),expand = expansion(mult=c(0,0.1)),labels = scales::number_format(accuracy = 0.1))+
  theme_classic()

FigR0=ggarrange(p1, p2, p3, p4, ncol=2, nrow=2,hjust=0, labels="auto",heights=c(1,1,1,1))
ggsave(file="Fig6.pdf", width = 7, height = 3.5)

#Efficacy
p5=ggplot() +
  geom_histogram(aes(x=sample_and_NLL[,5],y = ..density..),binwidth = 1/nbins, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,5],y = ..density..),binwidth = 1/nbins, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Efficacy of chlorination")+xlim(c(0,1))+
  scale_y_continuous(limits=c(0,3.5),expand = c(0, 0.1))+
  geom_vline(xintercept = median(resample_sample[,5]),col = "grey10")+
  theme_classic()

p6=ggplot() +
  geom_histogram(aes(x=sample_and_NLL[,7],y = ..density..),binwidth = 1/nbins, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,10]*resample_sample[,7],y = ..density..),binwidth = 1/nbins, fill = "blue", color = "black",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,7],y = ..density..),binwidth = 1/nbins, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Efficacy of latrine water seal")+xlim(c(0,1))+
  scale_y_continuous(limits=c(0,10),expand = c(0, 0.1))+
  geom_vline(xintercept = median(resample_sample[,7]),col = "grey10")+
  geom_vline(xintercept = median(resample_sample[,10]*resample_sample[,7]),col = "blue")+
  theme_classic()

p7=ggplot() +
  geom_histogram(aes(x=sample_and_NLL[,6],y = ..density..),binwidth = 1/nbins, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,9]*resample_sample[,6],y = ..density..),binwidth = 1/nbins, fill = "blue", color = "black",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,6],y = ..density..),binwidth = 1/nbins, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Efficacy of handwashing")+xlim(c(0,1))+
  scale_y_continuous(limits=c(0,3.5),expand = c(0, 0.1))+
  geom_vline(xintercept = median(resample_sample[,6]),col = "grey10")+
  geom_vline(xintercept = median(resample_sample[,9]*resample_sample[,6]),col = "blue")+
  theme_classic()

p8=ggplot() +
  geom_histogram(aes(x=sample_and_NLL[,4],y = ..density..),binwidth = 1/nbins, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,4],y = ..density..),binwidth = 1/nbins, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Efficacy of nutrition")+xlim(c(0,1))+
  scale_y_continuous(limits=c(0,6),expand = c(0, 0.1))+
  geom_vline(xintercept = median(resample_sample[,4]),col = "grey10")+
  theme_classic()

Fig_eff=ggarrange(p5, p6, p7, p8, ncol=2, nrow=2,hjust=-5, labels="auto",heights=c(1,1,1,1))
ggsave(file="Fig7.pdf", width = 7, height = 3.5)

#Relative R0s
q1=ggplot() +
  geom_histogram(aes(x=resample_sample[,11],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of midline")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = median(resample_sample[,11]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q2=ggplot() +
  geom_histogram(aes(x=resample_sample[,12],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of endline")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = median(resample_sample[,12]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q3=ggplot() +
  geom_histogram(aes(x=resample_sample[,17],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of W arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = median(resample_sample[,17]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q4=ggplot() +
  geom_histogram(aes(x=resample_sample[,16],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of S arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = median(resample_sample[,16]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q5=ggplot() +
  geom_histogram(aes(x=resample_sample[,13],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of H arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = median(resample_sample[,13]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q6=ggplot() +
  geom_histogram(aes(x=resample_sample[,18],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of WSH arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = median(resample_sample[,18]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q7=ggplot() +
  geom_histogram(aes(x=resample_sample[,14],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of N arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = median(resample_sample[,14]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q8=ggplot() +
  geom_histogram(aes(x=resample_sample[,15],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of WSHN arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = median(resample_sample[,15]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

Fig_relR0=ggarrange(q1, q2, q3, q4, q5, q6, q7, q8, ncol=2, nrow=4,hjust=c(-6, -5, -6,-5,-6, -10,-6,-5), labels="auto",heights=c(1,1,1,1))
ggsave(file="FigS3.pdf", width = 7, height = 7)

#Coverage
p9=ggplot() +
  geom_histogram(aes(x=sample_and_NLL[,8],y = ..density..),binwidth = 0.02, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,8],y = ..density..),binwidth = 0.02, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Community coverage")+xlim(c(0,0.2))+
  scale_y_continuous(expand = expansion(mult=c(0, 0.1)))+
  geom_vline(xintercept = median(resample_sample[,8]),col = "grey10")+
  theme_classic()
ggsave(file="FigS4.pdf", width = 3.5, height = 3.5/2)
