library(SobolSequence)
library(tictoc)
library(deSolve)
library(parallel)
library(ipw)
library(ggplot2)
library(ggpubr)

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

Fig_adh=ggarrange(DW, DS, DH, DN, ncol=2, nrow=2,hjust=c(-7.5,-7,-7.5,-7), labels="auto",heights=c(1,1,1,1))



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

###################################################################################
#Differential equation model

coverage_model = function(t, x, model_par){
  
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
  
  S = x[17]
  S_N = x[18]
  S_W= x[19]
  S_H= x[20]
  S_S= x[21]
  S_NW= x[22]
  S_NH= x[23]
  S_NS= x[24]
  S_WH= x[25]
  S_WS= x[26]
  S_HS= x[27]
  S_NWH= x[28]
  S_NWS= x[29]
  S_NHS= x[30]
  S_WHS= x[31]
  S_NWHS= x[32]
  
  EW = phiS * sum(x[c(5,8,10,11,13:16)])  + sum(x[c(1:4,6:7,9,12)])
  EH = sum(x[1:16])
  EO = sum(x[1:16])
  
  dxdt = numeric(length(x))
  dxdt[1]  =  S * ( R0W * EW + R0H * EH + R0O * EO) -  I
  dxdt[2]  =  S_N * phiN * ( R0W * EW  + R0H * EH + R0O * EO) -  I_N
  dxdt[3]  =  S_W * ( R0W * EW * phiW + R0H * EH + R0O * EO) -  I_W
  dxdt[4]  =  S_H * ( R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_H
  dxdt[5]  =  S_S * ( R0W * EW + R0H * EH + R0O * EO) -  I_S
  dxdt[6]  =  S_NW * phiN * ( R0W * EW * phiW + R0H * EH + R0O * EO) -  I_NW
  dxdt[7]  =  S_NH * phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_NH
  dxdt[8]  =  S_NS * phiN * ( R0W * EW + R0H * EH + R0O * EO) -  I_NS
  dxdt[9]  =  S_WH * ( R0W * phiW * EW + R0H * EH * phiH  + R0O * EO) -  I_WH
  dxdt[10] =  S_WS * ( R0W * phiW * EW + R0H * EH + R0O * EO) -  I_WS
  dxdt[11] =  S_HS * ( R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_HS
  dxdt[12] =  S_NWH * phiN * ( R0W * EW * phiW + R0H * EH * phiH  + R0O * EO) -  I_NWH
  dxdt[13] =  S_NWS * phiN * ( R0W * EW * phiW + R0H * EH + R0O * EO) -  I_NWS
  dxdt[14] =  S_NHS * phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) -  I_NHS
  dxdt[15] =  S_WHS * ( R0W * EW + R0H * phiW * EH * phiH  + R0O * EO) -  I_WHS
  dxdt[16] =  S_NWHS * phiN * ( R0W * EW * phiW + R0H * EH * phiH  + R0O * EO) -  I_NWHS
  
  dxdt[17] =   I - S * ( R0W * EW + R0H * EH + R0O * EO) 
  dxdt[18] =   I_N - S_N * phiN * ( R0W * EW + R0H * EH + R0O * EO) 
  dxdt[19] =   I_W - S_W * ( R0W * EW * phiW + R0H * EH + R0O * EO) 
  dxdt[20] =   I_H - S_H * ( R0W * EW + R0H * EH * phiH  + R0O * EO) 
  dxdt[21] =   I_S - S_S * ( R0W * EW + R0H * EH + R0O * EO) 
  dxdt[22] =   I_NW - S_NW * phiN * ( R0W * EW * phiW + R0H * EH + R0O * EO)
  dxdt[23] =   I_NH - S_NH * phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) 
  dxdt[24] =   I_NS - S_NS * phiN * ( R0W * EW + R0H * EH + R0O * EO) 
  dxdt[25] =   I_WH - S_WH * ( R0W * phiW * EW + R0H * EH * phiH  + R0O * EO) 
  dxdt[26] =   I_WS - S_WS * ( R0W * phiW * EW + R0H * EH + R0O * EO) 
  dxdt[27] =   I_HS - S_HS * ( R0W * EW + R0H * EH * phiH  + R0O * EO) 
  dxdt[28] =   I_NWH - S_NWH * phiN * ( R0W * EW * phiW + R0H * EH * phiH  + R0O * EO) 
  dxdt[29] =   I_NWS - S_NWS * phiN * ( R0W * EW * phiW + R0H * EH + R0O * EO) 
  dxdt[30] =   I_NHS - S_NHS * phiN * ( R0W * EW + R0H * EH * phiH  + R0O * EO) 
  dxdt[31] =   I_WHS - S_WHS * ( R0W * EW + R0H * phiW * EH * phiH  + R0O * EO) 
  dxdt[32] =   I_NWHS - S_NWHS * phiN * ( R0W * EW * phiW + R0H * EH * phiH  + R0O * EO)
  
  return(list(dxdt))
  
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
  #Par 9 and 10 are the relative R0 of midline and endline vs baseline
  #Par 11-16 are the relative R0 of each arm relative to the control arm
  
  
  data_cluster_temp = cluster_list[[k]][[1]]
  j = cluster_list[[k]][[2]]
  
  if (nrow(data_cluster_temp)==0){return(0)} else{ #skip the rest if there is no eligible people in this cluster for this survey
    
    #Set arm and svy specific R0 adjustment
    svy_temp = j+1
    arm_temp = data_cluster_temp$armid[1]
    
    R0_adj = c(1,par[9:10])[svy_temp]*c(1,par[11:16])[arm_temp]
    
    model_par= c(piN, piW, piH, piS, R0W*R0_adj, R0H*R0_adj, R0O*R0_adj)
    
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
    
    prev= 1E-3 #initial condition prevalence 
    x0 = c(rep(prev,16)*(rho_vec*omega + baseline_adherence_vec*(1-omega)),
           rep(1-prev, 16)*(rho_vec*omega + baseline_adherence_vec*(1-omega)))
    
    out_coverage = ode(x0, times = c(0,100), coverage_model, model_par,method="vode")
    
    steady_state = as.vector(pmax(tail(out_coverage[,2:17],1),1E-10))
    #It is possible for the  simulations go negative as they straddle 0. 
    #Capping the negative values addresses this issue 

    #Calculate prevalence within each group
    steady_state_normalized = steady_state/(omega*rho_vec + (1-omega)*baseline_adherence_vec)
    #Anything that is going to 0, set to very small to avoid dependence on the specific times in the ODE.
    steady_state_normalized[steady_state_normalized<0.005]=1E-10 
    
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
  par_R0=2*exp(logit2_par_R0)/(1+exp(logit2_par_R0))
  NLL = fullLogLikelihood(c(par_R0[1], par_other,par_R0[2:9]),list)
  
  return(NLL)
  
}

compute_NLL_function=function(x){
  
  par_other = scaled_sample[x,]
  logit2_par_R0_init = c(0.6, rep(0,8))
  
  out = nlm(f=fullLogLikelihood_wrapper, p=logit2_par_R0_init, par_other=par_other, list=cluster_list,steptol=1e-4,gradtol=1e-4)
  
  par_R0 = 2*exp(out$estimate)/(1+exp(out$estimate))
  
  print(x)
  
  return( c(par_R0[1],par_other,par_R0[2:9], out$minimum))
}


###################################################################################

set.seed(0)
nsample = 1E3*25
sobol_sample = sobolSequence.points(7,count=nsample)
paramlow =  c( 0,  0,  0,    0,     0,   0, 1E-5)
paramhigh = c(1,  1,    1,  1,     1,  1,  0.1)
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
clusterExport(cluster,"coverage_model")
clusterEvalQ(cluster,library("deSolve"))

tic()
sample_and_NLL1 = parSapply(cluster,1:25000,FUN=compute_NLL_function)
toc()
saveRDS(sample_and_NLL1,"sample_and_NLL_cluster_coverage.RDS")

stopCluster(cluster)
print("stopping cluster")

###################################################################################
#Sampling-importance resampling

sample_and_NLL = as.data.frame(t(readRDS("sample_and_NLL_cluster_coverage.RDS")))
names(sample_and_NLL)[17]="NLL"
sample_and_NLL$rNLL = sample_and_NLL$NLL- min(sample_and_NLL$NLL)
sample_and_NLL$prob = (exp(-sample_and_NLL$rNLL))

set.seed(0)
resample = sample(1:nrow(sample_and_NLL), 2.5E4, replace=TRUE, prob=sample_and_NLL$prob)
# saveRDS(resample,"resample.RDS")
resample = readRDS("resample.RDS")
index= unique(resample)

###################################################################################
#Calculate prevalence in each arm at each svy for the given parameters

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
      R0_adj = c(1,par[9:10])[svy_temp]*c(1,par[11:16])[arm_temp]
      
      model_par= c(piN, piW, piH, piS, R0W*R0_adj, R0H*R0_adj, R0O*R0_adj)
      
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
      
      prev= 0.06 #baseline prevalence 
      x0 = c(rep(prev,16)*(rho_vec*omega + baseline_adherence_vec*(1-omega)),
             rep(1-prev, 16)*(rho_vec*omega + baseline_adherence_vec*(1-omega)))
      
      out_coverage = ode(x0, times = c(0,100), coverage_model, model_par,method="vode")
      
      steady_state = as.vector(pmax(tail(out_coverage[,2:17],1),1E-10))
      #It is possible for the  simulations go negative as they straddle 0. 
      #Capping the negative values addresses this issue 
      
      #Calculate prevalence within each group
      steady_state_normalized = steady_state/(omega*rho_vec + (1-omega)*baseline_adherence_vec)
      #Anything that is going to 0, set to very small to avoid dependence on the specific times in the ODE.
      steady_state_normalized[steady_state_normalized<0.005]=1E-10
      
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
  prev_mat = prevalence_by_arm(par=unlist(sample_and_NLL[i,1:16]), dataset=data)

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
prevC12    = prevalences[,4]
prevH12    = prevalences[,8]
prevN12    = prevalences[,12]
prevWSHN12 = prevalences[,16]
prevS12    = prevalences[,20]
prevW12    = prevalences[,24]
prevWSH12  = prevalences[,28]

obs_prev = matrix(NA,7,2)
for (i in 1:7){
  obs_prev[i,1] = sum(data[data$armid==i & data$svy==0,"diar7d"])/length(data[data$armid==i & data$svy==0,"diar7d"])
  obs_prev[i,2] = sum(data[data$armid==i & data$svy%in%c(1,2),"diar7d"])/length(data[data$armid==i & data$svy%in%c(1,2),"diar7d"])
}
obs_prev_se = matrix(NA,7,2)
for (i in 1:7){
  p1= sum(data[data$armid==i & data$svy==0,"diar7d"])/length(data[data$armid==i & data$svy==0,"diar7d"])
  p2= sum(data[data$armid==i & data$svy%in%c(1,2),"diar7d"])/length(data[data$armid==i & data$svy%in%c(1,2),"diar7d"])
  obs_prev_se[i,1] = sqrt(p1*(1-p1)/nrow(data[data$armid==i & data$svy==0,]))
  obs_prev_se[i,2] = sqrt(p2*(1-p2)/nrow(data[data$armid==i & data$svy%in%c(1,2),]))
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
prev = c(obs_prev[c(1,6,5,2,7,3,4),1],obs_prev[c(1,6,5,2,7,3,4),2])
prev_low = c(obs_prev[c(1,6,5,2,7,3,4),1],obs_prev[c(1,6,5,2,7,3,4),2]) -qnorm(0.975,0,1)*c(obs_prev_se[c(1,6,5,2,7,3,4),1],obs_prev_se[c(1,6,5,2,7,3,4),2])
prev_high = c(obs_prev[c(1,6,5,2,7,3,4),1],obs_prev[c(1,6,5,2,7,3,4),2]) +qnorm(0.975,0,1)*c(obs_prev_se[c(1,6,5,2,7,3,4),1],obs_prev_se[c(1,6,5,2,7,3,4),2])
Arm = rep(c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"),2)
Data = rep(c("Baseline","Midline/endline"),each=7)
DATA1=data.frame(prev,prev_low,prev_high,Arm,Data)
DATA1$Arm=factor(DATA1$Arm,levels=c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"))
DATA1$Data=factor(DATA1$Data,levels=c("Baseline","Midline/endline"))

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
DATA3=data.frame(Prevalence,Arm,Simulation)
DATA3$Arm=factor(DATA3$Arm,levels=c("Control","Water","Sanitation","Hygiene","WSH","Nutrition","WSH-N"))
DATA3$Simulation=factor(DATA3$Simulation,levels=c("Baseline","Midline/endline"))

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
Pplot= ggplot()+
  geom_violin(data = DATA3, aes(x=Arm, y=Prevalence, fill=Simulation), adjust = 4,scale="count",width=2.5,position = position_dodge(width =1), color = NA)+
  scale_fill_manual(values=c("firebrick2","dodgerblue1"))+
  geom_pointrange(data = DATA1, aes(x=Arm, y=prev, ymin = prev_low, ymax = prev_high, color=Data),position = position_dodge(width = 0.5))+
  scale_color_manual(values=c("firebrick2","dodgerblue1"))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.12))+
  ylab("Prevalence")+
  theme_classic()
Pplot

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

###################################################################################
#Parameter histograms
resample_sample = sample_and_NLL[resample,]

#NLL
p0 = ggplot() + 
  geom_histogram(aes(x=resample_sample$NLL,y = ..density..),binwidth = 1, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=sample_and_NLL$NLL,y = ..density..),binwidth = 1, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Negative log-likelihood")+xlim(c(3245,3300))+
  scale_y_continuous(expand = expansion(mult=c(0, 0.1)))+
  annotate("text",x=3250, y=0.28, hjust=0, label = paste("Posterior"))+
  annotate("text",x=3260, y=0.05, hjust=0, label = paste("Prior"))+
  theme_classic()

#R0 and pathway-specific R0s
p1= ggplot() + 
  geom_histogram(aes(x=resample_sample[,1],y = ..density..),binwidth = 0.05, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  geom_vline(xintercept = mean(resample_sample[,1]),col = "grey10")+
  ylab("Density") + xlab("Overall R0")+
  xlim(c(0,1.4))+ #xlim(c(1,1.4))+
  scale_y_continuous(limits=c(0,13),expand = expansion(mult=c(0,0.1)),labels = scales::number_format(accuracy = 0.1))+
  theme_classic()

p2=ggplot() +
  geom_histogram(aes(x=sample_and_NLL[,1]*sample_and_NLL[,2],y = ..density..),binwidth = 0.05, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,1]*resample_sample[,2],y = ..density..),binwidth = 0.05, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  geom_vline(xintercept = mean(resample_sample[,1]*resample_sample[,2]),col = "grey10")+
  ylab("Density") + xlab("Water R0")+
  xlim(c(0,1.4))+
  scale_y_continuous(limits=c(0,13),expand = expansion(mult=c(0,0.1)),labels = scales::number_format(accuracy = 0.1))+
  theme_classic()

p3=ggplot() + 
  geom_histogram(aes(x=sample_and_NLL[,1]*(1-sample_and_NLL[,2])*sample_and_NLL[,3],y = ..density..),binwidth = 0.05, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,1]*(1-resample_sample[,2])*resample_sample[,3],y = ..density..),binwidth = 0.05, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  geom_vline(xintercept = mean(resample_sample[,1]*(1-resample_sample[,2])*resample_sample[,3]),col = "grey10")+
  ylab("Density") + xlab("Fomite R0")+
  xlim(c(0,1.4))+
  scale_y_continuous(limits=c(0,13),expand = expansion(mult=c(0,0.1)),labels = scales::number_format(accuracy = 0.1))+
  theme_classic()

p4=ggplot() + 
  geom_histogram(aes(x=sample_and_NLL[,1]*(1-sample_and_NLL[,2])*(1-sample_and_NLL[,3]),y = ..density..),binwidth = 0.05, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,1]*(1-resample_sample[,2])*(1-resample_sample[,3]),y = ..density..),binwidth = 0.05, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  geom_vline(xintercept = mean(resample_sample[,1]*(1-resample_sample[,2])*(1-resample_sample[,3])),col = "grey10")+
  ylab("Density") + xlab("Other R0")+
  xlim(c(0,1.4))+
  scale_y_continuous(limits=c(0,13),expand = expansion(mult=c(0,0.1)),labels = scales::number_format(accuracy = 0.1))+
  theme_classic()

FigR0=ggarrange(p1, p2, p3, p4, ncol=2, nrow=2,hjust=0, labels="auto",heights=c(1,1,1,1))

#Efficacy
p5=ggplot() + 
  geom_histogram(aes(x=sample_and_NLL[,5],y = ..density..),binwidth = 0.05, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,5],y = ..density..),binwidth = 0.05, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Efficacy of chlorination")+xlim(c(0,1))+
  scale_y_continuous(limits=c(0,3.5),expand = c(0, 0.1))+
  geom_vline(xintercept = mean(resample_sample[,5]),col = "grey10")+
  theme_classic()

p6=ggplot() + 
  geom_histogram(aes(x=sample_and_NLL[,7],y = ..density..),binwidth = 0.05, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,7],y = ..density..),binwidth = 0.05, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Efficacy of latrine water seal")+xlim(c(0,1))+
  scale_y_continuous(limits=c(0,3.5),expand = c(0, 0.1))+
  geom_vline(xintercept = mean(resample_sample[,7]),col = "grey10")+
  theme_classic()

p7=ggplot() + 
  geom_histogram(aes(x=sample_and_NLL[,6],y = ..density..),binwidth = 0.05, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,6],y = ..density..),binwidth = 0.05, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Efficacy of handwashing")+xlim(c(0,1))+
  scale_y_continuous(limits=c(0,3.5),expand = c(0, 0.1))+
  geom_vline(xintercept = mean(resample_sample[,6]),col = "grey10")+
  theme_classic()

p8=ggplot() + 
  geom_histogram(aes(x=sample_and_NLL[,4],y = ..density..),binwidth = 0.05, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,4],y = ..density..),binwidth = 0.05, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Efficacy of nutrition")+xlim(c(0,1))+
  scale_y_continuous(limits=c(0,5),expand = c(0, 0.1))+
  geom_vline(xintercept = mean(resample_sample[,4]),col = "grey10")+
  theme_classic()

Fig_eff=ggarrange(p5, p6, p7, p8, ncol=2, nrow=2,hjust=-5, labels="auto",heights=c(1,1,1,1))

#Coverage
p9=ggplot() + 
  geom_histogram(aes(x=sample_and_NLL[,8],y = ..density..),binwidth = 0.01, fill = "grey90", color = "grey50",alpha=0.75,boundary=0)+
  geom_histogram(aes(x=resample_sample[,8],y = ..density..),binwidth = 0.01, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Coverage fraction")+xlim(c(0,0.1))+
  scale_y_continuous(expand = expansion(mult=c(0, 0.1)))+
  geom_vline(xintercept = mean(resample_sample[,8]),col = "grey10")+
  theme_classic()

#Relative R0s
q1=ggplot() + 
  geom_histogram(aes(x=resample_sample[,9],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of midline")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = mean(resample_sample[,9]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q2=ggplot() + 
  geom_histogram(aes(x=resample_sample[,10],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of endline")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = mean(resample_sample[,10]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q3=ggplot() + 
  geom_histogram(aes(x=resample_sample[,15],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of W arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = mean(resample_sample[,15]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q4=ggplot() + 
  geom_histogram(aes(x=resample_sample[,14],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of S arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = mean(resample_sample[,14]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q5=ggplot() + 
  geom_histogram(aes(x=resample_sample[,11],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of H arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = mean(resample_sample[,11]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q6=ggplot() + 
  geom_histogram(aes(x=resample_sample[,16],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of WSH arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = mean(resample_sample[,16]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

q7=ggplot() + 
  geom_histogram(aes(x=resample_sample[,12],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of N arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = mean(resample_sample[,12]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic() 

q8=ggplot() + 
  geom_histogram(aes(x=resample_sample[,13],y = ..density..),binwidth = 0.002, fill = "grey30", color = "black",alpha=0.75,boundary=0)+
  ylab("Density") + xlab("Relative R0 of WSHN arm")+
  geom_vline(xintercept = 1,col = "grey30",linetype="dashed")+
  geom_vline(xintercept = mean(resample_sample[,13]),col = "grey10")+
  xlim(c(0.97,1.05))+
  scale_y_continuous(limits=c(0,210), expand = c(0, 0.1))+
  theme_classic()

Fig_relR0=ggarrange(q1, q2, q3, q4, q5, q6, q7, q8, ncol=2, nrow=4,hjust=c(-6, -5, -6,-5,-6, -10,-6,-5), labels="auto",heights=c(1,1,1,1))
