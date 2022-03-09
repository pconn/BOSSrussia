# create sensitivity tables

load('sens_out_wBS2012.RData')

#Base model extrapolation calcs
Max_obs= apply(Base$Z_s[,Data$S_i],1,"max")
Extrap_n_base = rep(0,4)
Extrap_ratio_base = rep(0,4)
for(isp in 1:4){
  Extrap_n_base[isp]=sum(Base$Z_s[isp,]>Max_obs[isp])
  Extrap_ratio_base[isp]= max(Base$Z_s[isp,])/Max_obs[isp]
}

Table_df_2012 = data.frame(
  Model=c("main","no misID","rd lair","rd Bering","-strata (bd)","-dist_land (bd)","-dist_shelf (bd)","-dist_edge (bd)",
          "-strata (rn)","-dist_shelf (rn)","+dist_edge (rn)","+dist_land (rn)",
          "-strata (rd)","-dist_land (rd)","-dist_shelf (rd)","-dist_edge (rd)",
          "-strata (sd)","-dist_shelf (sd)","+dist_edge (sd)","+dist_land (sd)"),
  LogL = rep(NA,20),k=rep(NA,20),AIC=rep(NA,20),Nbd=rep(NA,20),Nrn=rep(NA,20),Nrd=rep(NA,20),Nsd=rep(NA,20))
Table_df_2012$k[1:4]=length(Base$beta_bd)+length(Base$beta_rd)+length(Base$beta_rn)+length(Base$beta_sd)
Table_df_2012$k[5:20]=Table_df_2012$k[1]+c(-9,-2,-2,-2,-2,-1,+1,+1,-9,-2,-2,-2,-2,-1,+1,+1)
Table_df_2012$LogL[1]=-Base$jnll
Table_df_2012$LogL[3:4]=-c(Report_chukchi_lair$jnll,Report_bering$jnll)
Table_df_2012$LogL[5:20]=-c(LogL[9:12],LogL[4:7],LogL[13:16],LogL[1:4])
aic_base = 2*Base$jnll+2*Table_df_2012$k[1]
Table_df_2012$AIC[1]=0
Table_df_2012$AIC[3:20]=aic_base-(-2*Table_df_2012$Log[3:20]+2*(Table_df_2012$k[3:20]))
Table_df_2012$Nbd[1]=Base$N[3]
Table_df_2012$Nrn[1]=Base$N[2]
Table_df_2012$Nrd[1]=Base$N[4]
Table_df_2012$Nsd[1]=Base$N[1]
Table_df_2012$Nbd[3:4]=c(Report_chukchi_lair$N[3],Report_bering$N[3])
Table_df_2012$Nrn[3:4]=c(Report_chukchi_lair$N[2],Report_bering$N[2])
Table_df_2012$Nrd[3:4]=c(Report_chukchi_lair$N[4],Report_bering$N[4])
Table_df_2012$Nsd[3:4]=c(Report_chukchi_lair$N[1],Report_bering$N[1])
Table_df_2012$Nbd[5:20] = c(N[3,9:12],N[3,5:8],N[3,13:16],N[3,1:4]) #reordering since species order in analysis was sd, rn, bd, rd
Table_df_2012$Nrd[5:20] = c(N[4,9:12],N[4,5:8],N[4,13:16],N[4,1:4])
Table_df_2012$Nrn[5:20] = c(N[2,9:12],N[2,5:8],N[2,13:16],N[2,1:4])
Table_df_2012$Nsd[5:20] = c(N[1,9:12],N[1,5:8],N[1,13:16],N[1,1:4])

xtable_2012 = xtable::xtable(Table_df_2012,digits=c(NA,NA,2,0,2,0,0,0,0))
print(xtable_2012, include.rownames=FALSE)

#extrapolation metrics
Table_df_2012 = data.frame(
  Model=c("main","no misID","rd lair","rd Bering","-strata (bd)","-dist_land (bd)","-dist_shelf (bd)","-dist_edge (bd)",
          "-strata (rn)","-dist_shelf (rn)","+dist_edge (rn)","+dist_land (rn)",
          "-strata (rd)","-dist_land (rd)","-dist_shelf (rd)","-dist_edge (rd)",
          "-strata (sd)","-dist_shelf (sd)","+dist_edge (sd)","+dist_land (sd)"),
  I.bd = rep(NA,20), I.rn =rep(NA,20), I.rd =rep(NA,20), I.sd =rep(NA,20), 
  Omega.bd= rep(NA,20), Omega.rn =rep(NA,20), Omega.rd =rep(NA,20), Omega.sd =rep(NA,20))

Max_obs= apply(Base$Z_s[,Data$S_i],1,"max")
Extrap_n_base = rep(0,4)
Extrap_ratio_base = rep(0,4)
for(isp in 1:4){
  Extrap_n_base[isp]=sum(Base$Z_s[isp,]>Max_obs[isp])
  Extrap_ratio_base[isp]= max(Base$Z_s[isp,])/Max_obs[isp]
}
Table_df_2012[1,2:5]=Extrap_n_base[c(3,2,4,1)]  #reordering species in alphabetical order
Table_df_2012[1,6:9]=Extrap_ratio_base[c(3,2,4,1)]  #reordering species in alphabetical order

#Extrapolation calcs for no misID, alternative ringed seal availabilty models (don't do no misID since it didn't converge)
Max_obs= apply(Report_chukchi_lair$Z_s[,Data$S_i],1,"max")
Extrap_n_cur = rep(0,4)
Extrap_ratio_cur = rep(0,4)
for(isp in 1:4){
  Extrap_n_cur[isp]=sum(Report_chukchi_lair$Z_s[isp,]>Max_obs[isp])
  Extrap_ratio_cur[isp]= max(Report_chukchi_lair$Z_s[isp,])/Max_obs[isp]
}
Table_df_2012[3,2:5]=Extrap_n_cur[c(3,2,4,1)]  #reordering species in alphabetical order
Table_df_2012[3,6:9]=Extrap_ratio_cur[c(3,2,4,1)]  #reordering species in alphabetical order

Max_obs= apply(Report_bering$Z_s[,Data$S_i],1,"max")
Extrap_n_cur = rep(0,4)
Extrap_ratio_cur = rep(0,4)
for(isp in 1:4){
  Extrap_n_cur[isp]=sum(Report_bering$Z_s[isp,]>Max_obs[isp])
  Extrap_ratio_cur[isp]= max(Report_bering$Z_s[isp,])/Max_obs[isp]
}
Table_df_2012[4,2:5]=Extrap_n_cur[c(3,2,4,1)]  #reordering species in alphabetical order
Table_df_2012[4,6:9]=Extrap_ratio_cur[c(3,2,4,1)]  #reordering species in alphabetical order

Table_df_2012[5:20,2:5] = cbind(Extrap_num[3,],Extrap_num[2,],Extrap_num[4,],Extrap_num[1,])
Table_df_2012[5:20,6:9] = cbind(Extrap_ratio[3,],Extrap_ratio[2,],Extrap_ratio[4,],Extrap_ratio[1,])

xtable_2012 = xtable::xtable(Table_df_2012,digits=c(NA,NA,0,0,0,0,2,2,2,2))
print(xtable_2012, include.rownames=FALSE)


#now, repeat for 2013
load('sens_out_wBS2013.RData')

#Base model extrapolation calcs
Max_obs= apply(Base$Z_s[,Data$S_i],1,"max")
Extrap_n_base = rep(0,4)
Extrap_ratio_base = rep(0,4)
for(isp in 1:4){
  Extrap_n_base[isp]=sum(Base$Z_s[isp,]>Max_obs[isp])
  Extrap_ratio_base[isp]= max(Base$Z_s[isp,])/Max_obs[isp]
}

Table_df_2013 = data.frame(
  Model=c("main","no misID","rd lair","rd Bering","-strata (bd)","-dist_land (bd)","-dist_shelf (bd)","+dist_edge (bd)",
          "-strata (rn)","-dist_shelf (rn)","+dist_edge (rn)","+dist_land (rn)",
          "-strata (rd)","-dist_land (rd)","-dist_shelf (rd)","-dist_edge (rd)",
          "-strata (sd)","-dist_shelf (sd)","+dist_edge (sd)","+dist_land (sd)"),
  LogL = rep(NA,20),k=rep(NA,20),AIC=rep(NA,20),Nbd=rep(NA,20),Nrn=rep(NA,20),Nrd=rep(NA,20),Nsd=rep(NA,20))
Table_df_2013$k[1:4]=length(Base$beta_bd)+length(Base$beta_rd)+length(Base$beta_rn)+length(Base$beta_sd)
Table_df_2013$k[5:20]=Table_df_2013$k[1]+c(-9,-2,-2,-2,-2,-1,+1,+1,-9,-2,-2,-2,-2,-1,+1,+1)
Table_df_2013$LogL[1]=-Base$jnll
Table_df_2013$LogL[2:4]=-c(Report_noMisID$jnll,Report_chukchi_lair$jnll,Report_bering$jnll)
Table_df_2013$LogL[5:20]=-c(LogL[9:12],LogL[4:7],LogL[13:16],LogL[1:4])
aic_base = 2*Base$jnll+2*Table_df_2013$k[1]
Table_df_2013$AIC[1]=0
Table_df_2013$AIC[2:20]=aic_base-(-2*Table_df_2013$Log[2:20]+2*(Table_df_2013$k[2:20]))
Table_df_2013$Nbd[1]=Base$N[3]
Table_df_2013$Nrn[1]=Base$N[2]
Table_df_2013$Nrd[1]=Base$N[4]
Table_df_2013$Nsd[1]=Base$N[1]
Table_df_2013$Nbd[2:4]=c(Report_noMisID$N[3],Report_chukchi_lair$N[3],Report_bering$N[3])
Table_df_2013$Nrn[2:4]=c(Report_noMisID$N[2],Report_chukchi_lair$N[2],Report_bering$N[2])
Table_df_2013$Nrd[2:4]=c(Report_noMisID$N[4],Report_chukchi_lair$N[4],Report_bering$N[4])
Table_df_2013$Nsd[2:4]=c(Report_noMisID$N[1],Report_chukchi_lair$N[1],Report_bering$N[1])
Table_df_2013$Nbd[5:20] = c(N[3,9:12],N[3,5:8],N[3,13:16],N[3,1:4]) #reordering since species order in analysis was sd, rn, bd, rd
Table_df_2013$Nrd[5:20] = c(N[4,9:12],N[4,5:8],N[4,13:16],N[4,1:4])
Table_df_2013$Nrn[5:20] = c(N[2,9:12],N[2,5:8],N[2,13:16],N[2,1:4])
Table_df_2013$Nsd[5:20] = c(N[1,9:12],N[1,5:8],N[1,13:16],N[1,1:4])

xtable_2013 = xtable::xtable(Table_df_2013,digits=c(NA,NA,2,0,2,0,0,0,0))
print(xtable_2013, include.rownames=FALSE)

#extrapolation metrics
Table_df_2013 = data.frame(
  Model=c("main","no misID","rd lair","rd Bering","-strata (bd)","-dist_land (bd)","-dist_shelf (bd)","+dist_edge (bd)",
          "-strata (rn)","-dist_shelf (rn)","+dist_edge (rn)","+dist_land (rn)",
          "-strata (rd)","-dist_land (rd)","-dist_shelf (rd)","-dist_edge (rd)",
          "-strata (sd)","-dist_shelf (sd)","+dist_edge (sd)","+dist_land (sd)"),
  I.bd = rep(NA,20), I.rn =rep(NA,20), I.rd =rep(NA,20), I.sd =rep(NA,20), 
  Omega.bd= rep(NA,20), Omega.rn =rep(NA,20), Omega.rd =rep(NA,20), Omega.sd =rep(NA,20))

Max_obs= apply(Base$Z_s[,Data$S_i],1,"max")
Extrap_n_base = rep(0,4)
Extrap_ratio_base = rep(0,4)
for(isp in 1:4){
  Extrap_n_base[isp]=sum(Base$Z_s[isp,]>Max_obs[isp])
  Extrap_ratio_base[isp]= max(Base$Z_s[isp,])/Max_obs[isp]
}
Table_df_2013[1,2:5]=Extrap_n_base[c(3,2,4,1)]  #reordering species in alphabetical order
Table_df_2013[1,6:9]=Extrap_ratio_base[c(3,2,4,1)]  #reordering species in alphabetical order

#Extrapolation calcs for no misID, alternative ringed seal availabilty models 
Max_obs= apply(Report_noMisID$Z_s[,Data$S_i],1,"max")
Extrap_n_cur = rep(0,4)
Extrap_ratio_cur = rep(0,4)
for(isp in 1:4){
  Extrap_n_cur[isp]=sum(Report_noMisID$Z_s[isp,]>Max_obs[isp])
  Extrap_ratio_cur[isp]= max(Report_noMisID$Z_s[isp,])/Max_obs[isp]
}
Table_df_2013[2,2:5]=Extrap_n_cur[c(3,2,4,1)]  #reordering species in alphabetical order
Table_df_2013[2,6:9]=Extrap_ratio_cur[c(3,2,4,1)]  #reordering species in alphabetical order

Max_obs= apply(Report_chukchi_lair$Z_s[,Data$S_i],1,"max")
Extrap_n_cur = rep(0,4)
Extrap_ratio_cur = rep(0,4)
for(isp in 1:4){
  Extrap_n_cur[isp]=sum(Report_chukchi_lair$Z_s[isp,]>Max_obs[isp])
  Extrap_ratio_cur[isp]= max(Report_chukchi_lair$Z_s[isp,])/Max_obs[isp]
}
Table_df_2013[3,2:5]=Extrap_n_cur[c(3,2,4,1)]  #reordering species in alphabetical order
Table_df_2013[3,6:9]=Extrap_ratio_cur[c(3,2,4,1)]  #reordering species in alphabetical order

Max_obs= apply(Report_bering$Z_s[,Data$S_i],1,"max")
Extrap_n_cur = rep(0,4)
Extrap_ratio_cur = rep(0,4)
for(isp in 1:4){
  Extrap_n_cur[isp]=sum(Report_bering$Z_s[isp,]>Max_obs[isp])
  Extrap_ratio_cur[isp]= max(Report_bering$Z_s[isp,])/Max_obs[isp]
}
Table_df_2013[4,2:5]=Extrap_n_cur[c(3,2,4,1)]  #reordering species in alphabetical order
Table_df_2013[4,6:9]=Extrap_ratio_cur[c(3,2,4,1)]  #reordering species in alphabetical order

Table_df_2013[5:20,2:5] = cbind(Extrap_num[3,],Extrap_num[2,],Extrap_num[4,],Extrap_num[1,])
Table_df_2013[5:20,6:9] = cbind(Extrap_ratio[3,],Extrap_ratio[2,],Extrap_ratio[4,],Extrap_ratio[1,])

xtable_2013 = xtable::xtable(Table_df_2013,digits=c(NA,NA,0,0,0,0,2,2,2,2))
print(xtable_2013, include.rownames=FALSE)
