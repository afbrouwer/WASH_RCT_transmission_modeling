library(SobolSequence)
library(tictoc)
library(deSolve)
library(parallel)
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

#Arm-level intervention adherence in baseline/midline
intervention_adherence = matrix(0,7,16)
intervention_adherence[1,match(names(table(data[data$svy!=0 & data$armid==1,"adh_cat"])),adh_cat_groups)] = table(data[data$svy!=0 & data$armid==1,"adh_cat"])/sum(table(data[data$svy!=0 & data$armid==1,"adh_cat"]))
intervention_adherence[2,match(names(table(data[data$svy!=0 & data$armid==2,"adh_cat"])),adh_cat_groups)] = table(data[data$svy!=0 & data$armid==2,"adh_cat"])/sum(table(data[data$svy!=0 & data$armid==2,"adh_cat"]))
intervention_adherence[3,match(names(table(data[data$svy!=0 & data$armid==3,"adh_cat"])),adh_cat_groups)] = table(data[data$svy!=0 & data$armid==3,"adh_cat"])/sum(table(data[data$svy!=0 & data$armid==3,"adh_cat"]))
intervention_adherence[4,match(names(table(data[data$svy!=0 & data$armid==4,"adh_cat"])),adh_cat_groups)] = table(data[data$svy!=0 & data$armid==4,"adh_cat"])/sum(table(data[data$svy!=0 & data$armid==4,"adh_cat"]))
intervention_adherence[5,match(names(table(data[data$svy!=0 & data$armid==5,"adh_cat"])),adh_cat_groups)] = table(data[data$svy!=0 & data$armid==5,"adh_cat"])/sum(table(data[data$svy!=0 & data$armid==5,"adh_cat"]))
intervention_adherence[6,match(names(table(data[data$svy!=0 & data$armid==6,"adh_cat"])),adh_cat_groups)] = table(data[data$svy!=0 & data$armid==6,"adh_cat"])/sum(table(data[data$svy!=0 & data$armid==6,"adh_cat"]))
intervention_adherence[7,match(names(table(data[data$svy!=0 & data$armid==7,"adh_cat"])),adh_cat_groups)] = table(data[data$svy!=0 & data$armid==7,"adh_cat"])/sum(table(data[data$svy!=0 & data$armid==7,"adh_cat"]))

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



###################################################################################
#Calculate prevalence in each arm at each svy for the given parameters

prevalence_C0 = function(par, dataset){
  
  dataset=dataset[dataset$svy==0 & dataset$armid==1,]
  
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
  for (j in 0){
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
      #intervention distribution * fraction covered + non-intervention distribution * fraction not covered
      I0 = prev*N_vec
      
      steady_state = abs(nleqslv(x=I0,fn=model,jac = model_Jac, 
                                 model_par=model_par, N_vec = N_vec,
                                 method="Broyden", control=list(allowSingular=1))$x)
      steady_state[N_vec=0]=0
      
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
  
  prev = sum(dataset$model_prev) / nrow(dataset)
  
  return(prev)
}


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
      #intervention distribution * fraction covered + non-intervention distribution * fraction not covered
      I0 = prev*N_vec
      
      steady_state = abs(nleqslv(x=I0,fn=model,jac = model_Jac, 
                                 model_par=model_par, N_vec = N_vec,
                                 method="Broyden", control=list(allowSingular=1))$x)
      steady_state[N_vec=0]=0
      
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


prevalence_by_arm_new_coverage = function(par, new_omega, dataset){
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
      
      if (j == 0){
        intervention_adherence_vec = baseline_adherence_vec
      } else {
        armid = data_cluster_temp$armid[1] # should all be the same
        intervention_adherence_vec = intervention_adherence[armid,]
      }
      
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
      N_vec= rho_vec*omega + intervention_adherence_vec*(new_omega- omega)+ baseline_adherence_vec*(1-new_omega)
      #intervention distribution (in specific cluster) * original fraction covered + 
      # intervention distribution (in arm overall) * (difference between original and new fraction covered) +
      # non-intervention distribution * fraction not covered
      I0 = prev*N_vec
      
      steady_state = abs(nleqslv(x=I0,fn=model,jac = model_Jac, 
                                 model_par=model_par, N_vec = N_vec,
                                 method="Broyden", control=list(allowSingular=1))$x)
      steady_state[N_vec=0]=0
      
      #Calculate prevalence within each group
      steady_state_normalized = steady_state/(omega*rho_vec+ intervention_adherence_vec*(new_omega- omega) +  (1-new_omega)*baseline_adherence_vec)
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

prevalence_by_arm_no_conditions = function(par, dataset){
  
  #Armid
  # 1: Control
  # 2: Handwashing
  # 3: Nutrition
  # 4: Nutrition + WSH
  # 5: Sanitation
  # 6: Water
  # 7: WSH
  
  #Change adherence group to have no baseline conditions
  new_baseline_adherence_vec = c(1,rep(0,15))
  
  dataset$adh_cat_new=dataset$adh_cat
  dataset[dataset$svy==0, "adh_cat_new"] = "C"
  dataset[dataset$armid ==1 & dataset$svy %in% c(1,2), "adh_cat_new"] = "C"
  dataset[dataset$armid ==2 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("C","N","W","S","NW","NS","WS","NWS"), "adh_cat_new"] = "C"
  dataset[dataset$armid ==2 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("H","NH","WH","HS","NWH","NHS","WHS","NWHS"), "adh_cat_new"] = "H"
  dataset[dataset$armid ==3 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("C","W","H","S","WH","WS","HS","WHS"), "adh_cat_new"] = "C"
  dataset[dataset$armid ==3 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("N","NW","NH","NS","NWH","NWS","NHS","NWHS"), "adh_cat_new"] = "N"
  #Don't need to do anything with arm 4 (WSH-N)
  dataset[dataset$armid ==5 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("C","N","W","H","NW","NH","WH","NWH"), "adh_cat_new"] = "C"
  dataset[dataset$armid ==5 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("S","NS","WS","HS","NWS","NHS","WHS","NWHS"), "adh_cat_new"] = "S"
  dataset[dataset$armid ==6 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("C","N","H","S","NH","NS","HS","NHS"), "adh_cat_new"] = "C"
  dataset[dataset$armid ==6 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("W","NW","WH","WS","NWH","NWS","WHS","NWHS"), "adh_cat_new"] = "W"
  #Don't need to do anything with arm 7 (WSH) because N isn't a condition.
  
  ##########################################################
  
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
      rho_vec = c(sum(data_cluster_temp$adh_cat_new == "C"),
                  sum(data_cluster_temp$adh_cat_new == "N"),
                  sum(data_cluster_temp$adh_cat_new == "W"),
                  sum(data_cluster_temp$adh_cat_new == "H"),
                  sum(data_cluster_temp$adh_cat_new == "S"),
                  sum(data_cluster_temp$adh_cat_new == "NW"),
                  sum(data_cluster_temp$adh_cat_new == "NH"),
                  sum(data_cluster_temp$adh_cat_new == "NS"),
                  sum(data_cluster_temp$adh_cat_new == "WH"),
                  sum(data_cluster_temp$adh_cat_new == "WS"),
                  sum(data_cluster_temp$adh_cat_new == "HS"),
                  sum(data_cluster_temp$adh_cat_new == "NWH"),
                  sum(data_cluster_temp$adh_cat_new == "NWS"),
                  sum(data_cluster_temp$adh_cat_new == "NHS"),
                  sum(data_cluster_temp$adh_cat_new == "WHS"),
                  sum(data_cluster_temp$adh_cat_new == "NWHS"))/length(data_cluster_temp$adh_cat)
      
      prev= max(1-1/sum(model_par[5:7]),1E-10) #initial condition prevalence 
      N_vec= rho_vec*omega + new_baseline_adherence_vec*(1-omega)
      I0 = prev*N_vec
      
      steady_state = abs(nleqslv(x=I0,fn=model,jac = model_Jac, 
                                 model_par=model_par, N_vec = N_vec,
                                 method="Broyden", control=list(allowSingular=1))$x)
      steady_state[N_vec=0]=0
      
      #Calculate prevalence within each group
      steady_state_normalized = steady_state/(omega*rho_vec + (1-omega)*new_baseline_adherence_vec)
      #Anything that is going to 0, set to very small to avoid dependence on the specific times in the ODE.
      steady_state_normalized[steady_state_normalized<0.005]=1E-10
      
      data_cluster_temp[data_cluster_temp$adh_cat_new == "C","model_prev"] = steady_state_normalized[1]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "N","model_prev"] = steady_state_normalized[2]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "W","model_prev"] = steady_state_normalized[3]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "H","model_prev"] = steady_state_normalized[4]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "S","model_prev"] = steady_state_normalized[5]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NW","model_prev"] = steady_state_normalized[6]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NH","model_prev"] = steady_state_normalized[7]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NS","model_prev"] = steady_state_normalized[8]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "WH","model_prev"] = steady_state_normalized[9]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "WS","model_prev"] = steady_state_normalized[10]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "HS","model_prev"] = steady_state_normalized[11]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NWH","model_prev"] = steady_state_normalized[12]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NWS","model_prev"] = steady_state_normalized[13]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NHS","model_prev"] = steady_state_normalized[14]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "WHS","model_prev"] = steady_state_normalized[15]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NWHS","model_prev"] = steady_state_normalized[16]
      
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

prevalence_by_arm_full_adherence = function(par, dataset){
  
  #Armid
  # 1: Control
  # 2: Handwashing
  # 3: Nutrition
  # 4: Nutrition + WSH
  # 5: Sanitation
  # 6: Water
  # 7: WSH
  
  #Change adherence group to have full adherence
  
  dataset$adh_cat_new=dataset$adh_cat
  dataset[dataset$armid ==2 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("C"), "adh_cat_new"] = "H"
  dataset[dataset$armid ==2 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("N"), "adh_cat_new"] = "NH"
  dataset[dataset$armid ==2 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("W"), "adh_cat_new"] = "WH"
  dataset[dataset$armid ==2 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("S"), "adh_cat_new"] = "HS"
  dataset[dataset$armid ==2 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("NW"), "adh_cat_new"] = "NWH"
  dataset[dataset$armid ==2 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("NS"), "adh_cat_new"] = "NHS"
  dataset[dataset$armid ==2 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("WS"), "adh_cat_new"] = "WHS"
  dataset[dataset$armid ==2 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("NWS"), "adh_cat_new"] = "NWHS"
  
  dataset[dataset$armid ==3 & dataset$svy %in% c(1,2) & grepl("T",dataset$individualid,fixed=T) & dataset$adh_cat %in%  c("C"), "adh_cat_new"] = "N"
  dataset[dataset$armid ==3 & dataset$svy %in% c(1,2) & grepl("T",dataset$individualid,fixed=T) & dataset$adh_cat %in%  c("W"), "adh_cat_new"] = "NW"
  dataset[dataset$armid ==3 & dataset$svy %in% c(1,2) & grepl("T",dataset$individualid,fixed=T) & dataset$adh_cat %in%  c("H"), "adh_cat_new"] = "NH"
  dataset[dataset$armid ==3 & dataset$svy %in% c(1,2) & grepl("T",dataset$individualid,fixed=T) & dataset$adh_cat %in%  c("S"), "adh_cat_new"] = "NS"
  dataset[dataset$armid ==3 & dataset$svy %in% c(1,2) & grepl("T",dataset$individualid,fixed=T) & dataset$adh_cat %in%  c("WH"), "adh_cat_new"] = "NWH"
  dataset[dataset$armid ==3 & dataset$svy %in% c(1,2) & grepl("T",dataset$individualid,fixed=T) & dataset$adh_cat %in%  c("WS"), "adh_cat_new"] = "NWS"
  dataset[dataset$armid ==3 & dataset$svy %in% c(1,2) & grepl("T",dataset$individualid,fixed=T) & dataset$adh_cat %in%  c("HS"), "adh_cat_new"] = "NHS"
  dataset[dataset$armid ==3 & dataset$svy %in% c(1,2) & grepl("T",dataset$individualid,fixed=T) & dataset$adh_cat %in%  c("WHS"), "adh_cat_new"] = "NWHS"
  
  dataset[dataset$armid ==4 & dataset$svy %in% c(1,2) & grepl("T",dataset$individualid,fixed=T), "adh_cat_new"] = "NWHS"
  dataset[dataset$armid ==4 & dataset$svy %in% c(1,2) & grepl("C",dataset$individualid,fixed=T), "adh_cat_new"] = "WHS"
  
  
  dataset[dataset$armid ==5 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("C"), "adh_cat_new"] = "S"
  dataset[dataset$armid ==5 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("N"), "adh_cat_new"] = "NS"
  dataset[dataset$armid ==5 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("W"), "adh_cat_new"] = "WS"
  dataset[dataset$armid ==5 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("H"), "adh_cat_new"] = "HS"
  dataset[dataset$armid ==5 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("NW"), "adh_cat_new"] = "NWS"
  dataset[dataset$armid ==5 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("NH"), "adh_cat_new"] = "NHS"
  dataset[dataset$armid ==5 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("WH"), "adh_cat_new"] = "WHS"
  dataset[dataset$armid ==5 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("NWH"), "adh_cat_new"] = "NWHS"
  
  dataset[dataset$armid ==6 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("C"), "adh_cat_new"] = "W"
  dataset[dataset$armid ==6 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("N"), "adh_cat_new"] = "NW"
  dataset[dataset$armid ==6 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("H"), "adh_cat_new"] = "WH"
  dataset[dataset$armid ==6 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("S"), "adh_cat_new"] = "WS"
  dataset[dataset$armid ==6 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("NH"), "adh_cat_new"] = "NWH"
  dataset[dataset$armid ==6 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("NS"), "adh_cat_new"] = "NWS"
  dataset[dataset$armid ==6 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("HS"), "adh_cat_new"] = "WHS"
  dataset[dataset$armid ==6 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("NHS"), "adh_cat_new"] = "NWHS"
  
  dataset[dataset$armid ==7 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("C","W","H","S","WH","WS","HS","WHS"), "adh_cat_new"] = "WHS"
  dataset[dataset$armid ==7 & dataset$svy %in% c(1,2) & dataset$adh_cat %in%  c("N","NW","NH","NS","NWH","NWS","NHS","NWHS"), "adh_cat_new"] = "NWHS" #This one shouldn't be needed but I'm writing it anyway
  
  ##########################################################
  
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
      rho_vec = c(sum(data_cluster_temp$adh_cat_new == "C"),
                  sum(data_cluster_temp$adh_cat_new == "N"),
                  sum(data_cluster_temp$adh_cat_new == "W"),
                  sum(data_cluster_temp$adh_cat_new == "H"),
                  sum(data_cluster_temp$adh_cat_new == "S"),
                  sum(data_cluster_temp$adh_cat_new == "NW"),
                  sum(data_cluster_temp$adh_cat_new == "NH"),
                  sum(data_cluster_temp$adh_cat_new == "NS"),
                  sum(data_cluster_temp$adh_cat_new == "WH"),
                  sum(data_cluster_temp$adh_cat_new == "WS"),
                  sum(data_cluster_temp$adh_cat_new == "HS"),
                  sum(data_cluster_temp$adh_cat_new == "NWH"),
                  sum(data_cluster_temp$adh_cat_new == "NWS"),
                  sum(data_cluster_temp$adh_cat_new == "NHS"),
                  sum(data_cluster_temp$adh_cat_new == "WHS"),
                  sum(data_cluster_temp$adh_cat_new == "NWHS"))/length(data_cluster_temp$adh_cat)
      
      prev= max(1-1/sum(model_par[5:7]),1E-10) #initial condition prevalence 
      N_vec= rho_vec*omega + baseline_adherence_vec*(1-omega)
      I0 = prev*N_vec
      
      steady_state = abs(nleqslv(x=I0,fn=model,jac = model_Jac, 
                                 model_par=model_par, N_vec = N_vec,
                                 method="Broyden", control=list(allowSingular=1))$x)
      steady_state[N_vec=0]=0
      
      #Calculate prevalence within each group
      steady_state_normalized = steady_state/(omega*rho_vec + (1-omega)*baseline_adherence_vec)
      #Anything that is going to 0, set to very small to avoid dependence on the specific times in the ODE.
      steady_state_normalized[steady_state_normalized<0.005]=1E-10
      
      data_cluster_temp[data_cluster_temp$adh_cat_new == "C","model_prev"] = steady_state_normalized[1]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "N","model_prev"] = steady_state_normalized[2]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "W","model_prev"] = steady_state_normalized[3]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "H","model_prev"] = steady_state_normalized[4]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "S","model_prev"] = steady_state_normalized[5]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NW","model_prev"] = steady_state_normalized[6]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NH","model_prev"] = steady_state_normalized[7]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NS","model_prev"] = steady_state_normalized[8]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "WH","model_prev"] = steady_state_normalized[9]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "WS","model_prev"] = steady_state_normalized[10]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "HS","model_prev"] = steady_state_normalized[11]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NWH","model_prev"] = steady_state_normalized[12]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NWS","model_prev"] = steady_state_normalized[13]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NHS","model_prev"] = steady_state_normalized[14]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "WHS","model_prev"] = steady_state_normalized[15]
      data_cluster_temp[data_cluster_temp$adh_cat_new == "NWHS","model_prev"] = steady_state_normalized[16]
      
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

counterfactual_simulation = function(i, simul){
  print(which(unique(index)==i))
  #factual
  if (simul == "factual"){ 
    par=unlist(sample_and_NLL[i,1:18])
    prev_mat = prevalence_by_arm(par=par, data)}
  #counterfactual 
  if (simul == "no_arm_time_var"){ 
    par=unlist(sample_and_NLL[i,1:18])
    par[11:18] = 1
    prev_mat = prevalence_by_arm(par=par, data)}
  #counterfactual - coverage = 20%
  if (simul == "coverage20"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    prev_mat = prevalence_by_arm_new_coverage(par=par, new_omega = 0.2, data)}
  #counterfactual - coverage = 30%
  if (simul == "coverage30"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    prev_mat = prevalence_by_arm_new_coverage(par=par, new_omega = 0.3, data)}
  #counterfactual - coverage = 40%
  if (simul == "coverage40"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    prev_mat = prevalence_by_arm_new_coverage(par=par, new_omega = 0.4, data)}
  #counterfactual - coverage = 50%
  if (simul == "coverage50"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    prev_mat = prevalence_by_arm_new_coverage(par=par, new_omega = 0.5, data)}
  #counterfactual - coverage = 60%
  if (simul == "coverage60"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    prev_mat = prevalence_by_arm_new_coverage(par=par, new_omega = 0.6, data)}
  #counterfactual - coverage = 70%
  if (simul == "coverage70"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    prev_mat = prevalence_by_arm_new_coverage(par=par, new_omega = 0.7, data)}
  #counterfactual - coverage = 80%
  if (simul == "coverage80"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    prev_mat = prevalence_by_arm_new_coverage(par=par, new_omega = 0.8, data)}
  #counterfactual - coverage = 90%
  if (simul == "coverage90"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    prev_mat = prevalence_by_arm_new_coverage(par=par, new_omega = 0.9, data)}
  #counterfactual - coverage = 100%
  if (simul == "coverage100"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    prev_mat = prevalence_by_arm_new_coverage(par=par, new_omega = 1, data)}
  # counterfactual - no conditions 
  if (simul == "no_conditions"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    prev_mat = prevalence_by_arm_no_conditions(par=par, data)}
  # counterfactual - full adherence 
  if (simul == "full_adherence"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    prev_mat = prevalence_by_arm_full_adherence(par=par, data)}
  # counterfactual - 50% less "other" transmission
  if (simul == "50_less_other"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    R0O = par[1]*(1-par[2])*(1-par[3]) * 0.5 
    R0W = (par[1] - par[1]*(1-par[2])*(1-par[3]) * 0.5) * (par[2])/(par[2]+(1-par[2])*par[3])
    R0H = (par[1] - par[1]*(1-par[2])*(1-par[3]) * 0.5) * ((1-par[2])*par[3])/(par[2]+(1-par[2])*par[3])
    prev_mat = prevalence_by_arm(c(par[1],(R0W/par[1]), R0H/(par[1]-R0W), par[4:18]), data)}
  # counterfactual - 100% less "other" transmission
  if (simul == "100_less_other"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    R0W = (par[1]*par[2])/(par[2]+(1-par[2])*par[3])
    #Note how R0H is determined by setting the new par[3] to 1 in the following argument
    prev_mat = prevalence_by_arm(c(par[1],(R0W/par[1]), 1, par[4:16]), data)}
  # counterfactual - 25% more efficacy of N
  if (simul == "25efficacy_N"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    par[4]= min(1, 1.25*par[4])
    prev_mat = prevalence_by_arm(par, data)}
  # counterfactual - 50% more efficacy of N
  if (simul == "50efficacy_N"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    par[4]= min(1, 1.5*par[4])
    prev_mat = prevalence_by_arm(par, data)}
  if (simul == "100efficacy_N"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    par[4]= min(1, 2*par[4])
    prev_mat = prevalence_by_arm(par, data)}
  # counterfactual - 25% more efficacy of W
  if (simul == "25efficacy_W"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    par[5]= min(1, 1.25*par[5])
    prev_mat = prevalence_by_arm(par, data)}
  # counterfactual - 50% more efficacy of W
  if (simul == "50efficacy_W"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    par[5]= min(1, 1.5*par[5])
    prev_mat = prevalence_by_arm(par, data)}
  if (simul == "100efficacy_W"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    par[5]= min(1, 2*par[5])
    prev_mat = prevalence_by_arm(par, data)}
  # counterfactual - 25% more efficacy of H
  if (simul == "25efficacy_H"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    new_par = min(1, 1.25*par[6])
    par[9]= par[9]*par[6]/new_par
    par[6]= new_par 
    prev_mat = prevalence_by_arm(par, data)}
  # counterfactual - 50% more efficacy of H
  if (simul == "50efficacy_H"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    new_par = min(1, 1.5*par[6])
    par[9]= par[9]*par[6]/new_par
    par[6]= new_par 
    prev_mat = prevalence_by_arm(par, data)}
  if (simul == "100efficacy_H"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    new_par = min(1, 2*par[6])
    par[9]= par[9]*par[6]/new_par
    par[6]= new_par 
    prev_mat = prevalence_by_arm(par, data)}
  # counterfactual - 25% more efficacy of S
  if (simul == "25efficacy_S"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    new_par = min(1, 1.25*par[7])
    par[10]= par[10]*par[7]/new_par
    par[7]= new_par 
    prev_mat = prevalence_by_arm(par, data)}
  # counterfactual - 50% more efficacy of S
  if (simul == "50efficacy_S"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    new_par = min(1, 1.5*par[7])
    par[10]= par[10]*par[7]/new_par
    par[7]= new_par     
    prev_mat = prevalence_by_arm(par, data)}
  # counterfactual - 100% more efficacy of S
  if (simul == "100efficacy_S"){ 
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    new_par = min(1, 2*par[7])
    par[10]= par[10]*par[7]/new_par
    par[7]= new_par 
    prev_mat = prevalence_by_arm(par, data)}
  # counterfactual - double baseline prevalence
  if (simul == "2prev"){
    par=unlist(sample_and_NLL[i,1:18])
    # par[11:18] = 1 
    base_prev = prevalence_C0(par,data)

    new_R0_func = function(par1){
      new_par= par
      new_par[1]=par1
      
      new_base_prev = prevalence_C0(new_par,data)
      
      cost = (base_prev*2 - new_base_prev)^2
      return(cost)
    }
    mle = optim(par[1]/(par[1]+2-2*par[1]),new_R0_func,method="Brent",lower=1,upper=2)
    par[1]=mle$par
    
    prev_mat = prevalence_by_arm(par, data=data)}
  
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

###################################################################################
#Sampling-importance resampling

sample_and_NLL = readRDS("sample_and_NLL.RDS")
resample = readRDS("resample.RDS")
index= unique(resample)
###################################################################################


n_cores = 120
print("starting cluster")
cluster = makeCluster(n_cores)
clusterExport(cluster,"data")
clusterExport(cluster,"sample_and_NLL")
clusterExport(cluster,"index")
clusterExport(cluster,"counterfactual_simulation")
clusterExport(cluster,"model")
clusterExport(cluster,"model_Jac")
clusterExport(cluster,"prevalence_C0")
clusterExport(cluster,"prevalence_by_arm")
clusterExport(cluster,"prevalence_by_arm_new_coverage")
clusterExport(cluster,"prevalence_by_arm_no_conditions")
clusterExport(cluster,"prevalence_by_arm_full_adherence")
clusterExport(cluster,"baseline_adherence_vec")
clusterExport(cluster,"intervention_adherence")
clusterEvalQ(cluster,library("deSolve"))
clusterEvalQ(cluster,library("nleqslv"))

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="factual")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_factual.RDS")
print("factuals")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="no_arm_time_var")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_no_arm_time_var.RDS")
print("no arm and time variation")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="coverage20")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_coverage20_with_var.RDS")
print("coverage20")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="coverage30")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_coverage30_with_var.RDS")
print("coverage30")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="coverage40")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_coverage40_with_var.RDS")
print("coverage40")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="coverage50")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_coverage50_with_var.RDS")
print("coverage50")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="coverage60")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_coverage60_with_var.RDS")
print("coverage60")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="coverage70")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_coverage70_with_var.RDS")
print("coverage70")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="coverage80")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_coverage80_with_var.RDS")
print("coverage80")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="coverage90")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_coverage90_with_var.RDS")
print("coverage90")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="coverage100")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_coverage100_with_var.RDS")
print("coverage100")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="no_conditions")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_no_conditions_with_var.RDS")
print("no conditions")

tic()
index_prevalences = parSapply(cluster,index,FUN=counterfactual_simulation, simul="full_adherence")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_full_adherence_with_var.RDS")
print("full adherence")

tic()
index_prevalences =  parSapply(cluster,index,FUN=counterfactual_simulation, simul="50_less_other")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_50_less_other_with_var.RDS")
print("50% less other")

tic()
index_prevalences =  parSapply(cluster,index,FUN=counterfactual_simulation, simul="100efficacy_N")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_100efficacy_N_with_var.RDS")
print("N100")

tic()
index_prevalences =  parSapply(cluster,index,FUN=counterfactual_simulation, simul="100efficacy_W")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_100efficacy_W_with_var.RDS")
print("W100")

tic()
index_prevalences =  parSapply(cluster,index,FUN=counterfactual_simulation, simul="100efficacy_S")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_100efficacy_S_with_var.RDS")
print("S100")

tic()
index_prevalences =  parSapply(cluster,index,FUN=counterfactual_simulation, simul="100efficacy_H")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_100efficacy_H_with_var.RDS")
print("H100")

tic()
index_prevalences =  parSapply(cluster,index,FUN=counterfactual_simulation, simul="2prev")
toc()
index_prevalences = cbind(index, t(index_prevalences))
full_prevalences = index_prevalences[match(resample,index),-1]
saveRDS(full_prevalences, "Bprevalences_2prev_with_var.RDS")
print("2prev")

stopCluster(cluster)
print("stopping cluster")

