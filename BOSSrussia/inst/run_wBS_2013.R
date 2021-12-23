


#library( RandomFields )
library( TMB )
#library(TMBdebug)


# Compile
setwd("C:/Users/paul.conn/git/BOSSrussia/")
TmbFile = "C:/Users/paul.conn/git/BOSSrussia/BOSSrussia/src/wBeringTMB"
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 


#require(STabundance)
source('./BOSSrussia/R/util_funcs.R')
load("Grid_wBS.RData")
Grid=Grid_wBS  #for plotting
load('wBS_2013_Data_TMB.RData')  #created with format_Okhotsk_TMB.R script in OkhotskST project
Data = Data_2013
rm(Data_2013)


# Init Parameter values
n_species=4

X_s = Data$X_s
X_s$fi2=X_s$fi_nic^2

### Use Irina's models from Chernook et al. 2018
X_s$Strata=factor(X_s$Strata)
X_dm = model.matrix(~Strata*(ice_conc+ice2+dist_mainland+dist_contour),X_s)
X_dm = X_dm[,-1]  #take out intercept
X_rd = X_bd = X_dm

X_dm = model.matrix(~Strata+ice_conc+ice2+dist_contour,X_s)

#X_sd = X_sd[,-c(6,8,10,12,14,15,16)]  #remove interactions in "north" strata since no seals observed there; also remove dist_shelf
X_rn = X_sd = X_dm[,-1]



Data$X_sd = as.matrix(X_sd)
Data$X_rn = as.matrix(X_rn)
Data$X_bd = as.matrix(X_bd)
Data$X_rd = as.matrix(X_rd)
#note not currently employing priors on regression params (issues with strata covariate)
Data$Beta_VCsd = matrix(100 * solve(crossprod(Data$X_sd)),ncol(Data$X_sd),ncol(Data$X_sd))
Data$Beta_VCrn = matrix(100 * solve(crossprod(Data$X_rn)),ncol(Data$X_rn),ncol(Data$X_rn))
Data$Beta_VCbd = matrix(100 * solve(crossprod(Data$X_bd)),ncol(Data$X_bd),ncol(Data$X_bd))
Data$Beta_VCrd = matrix(100 * solve(crossprod(Data$X_rd)),ncol(Data$X_rd),ncol(Data$X_rd))
Data$phi_max = 10

Data$Melt_i = Data$S_i %/% Data$n_s -15 #just set = to day of survey relative to May 5

#Data$Beta_VC = as.matrix(Matrix::bdiag(Beta_VCsd,Beta_VCrn,Beta_VCbd,Beta_VCrd))

Data$C_i= as.matrix(Data$C_i)
Data$Ice_s = Data$X_s[,"ice_conc"]
ice_thresh = 0.01  
Data$Which_zero = which(Data$X_s[,"ice_conc"]<0.01)-1  #


Params = list("log_N" = rep(log(50000),n_species),"beta_sd"=rep(0,ncol(Data$X_sd)),"beta_rn"=rep(0,ncol(Data$X_rn)),
              "beta_bd"=rep(0,ncol(Data$X_bd)),"beta_rd"=rep(0,ncol(Data$X_rd)),
              "thin_beta_day" = rep(0,2*n_species),  
              "phi_logit"=rep(0,n_species+3),"p_logit"=rep(0,n_species+3),
                "thin_logit_i"=Data$thin_mu_logit,"MisID_pars"=Data$MisID_mu,
                "logit_Pup_prop"=Data$Pup_prop_mu,"beta0_ringed"=0,"beta1_ringed"=0)


  # Random
  #Random = c( "Etainput_s","thin_logit_i","MisID_pars","logit_Pup_prop" )
  #Random = c( "Etainput_s","thin_logit_i","logit_Pup_prop" )
  #Random = c("Etainput_s")
  #Random = c("beta","thin_logit_i","logit_Pup_prop")
  Random = c("logit_Pup_prop") 
  #Random = c("logit_Pup_prop","beta_sd","beta_rn","beta_bd","beta_rd")
  #Random = NULL
   #Random = list()
  
  # Fix parameters
  Map = list()
  #Map[["phi_log"]]=factor(rep(1,n_species+3))
  #Map[["phi_log"]]=factor(rep(2,n_species+3))
  #Map[["p_logit"]]=factor(rep(1,n_species+3))
  Map[["MisID_pars"]]=factor(rep(NA,length(Data$MisID_mu)))
  Map[["thin_beta_day"]]=factor(rep(NA,2*n_species))
  Map[["thin_logit_i"]]=factor(rep(NA,length(Data$thin_mu_logit)))
  #Map[["Etainput_s"]]=factor(rep(NA,length(Params$Etainput_s)))
  #Map[["logit_Pup_prop"]]=factor(rep(NA,length(Params$logit_Pup_prop)))
  
  
  # Make object
  #compile( paste0(Version,".cpp") )
  dyn.load( dynlib(TmbFile) )
  Start_time = Sys.time()
  #setwd( "C:/Users/paul.conn/git/OkhotskST/OkhotskSeal/src/")
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
  #Obj = MakeADFun( data=Data, parameters=Params, map=Map, silent=FALSE)
  Obj$fn( Obj$par )
  
  # Run
  #Lower = -Inf
  #Upper = Inf
  Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
  Upper = 50
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max = 10000,maxit=10000))         #
  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  Report = Obj$report()
  
  
  #XY = coordinates(Grid)
  #rownames(XY) = c(1:nrow(XY))
  #plot_N_map_xy(N=Report$Z_s[1,],XY=XY) 
  
  
  source('c:/users/paul.conn/git/ChukchiPower/ChukchiPower/R/util_funcs.R')
  library(rgeos)
  library(RColorBrewer)
  library(ggplot2)
  it=10
  plot_N_map(N=Report$Z_s[1,(it-1)*Data$n_s+c(1:Data$n_s)],Grid=Grid) 
  plot_N_map(N=Report$Z_s[2,(it-1)*Data$n_s+c(1:Data$n_s)],Grid=Grid) 
  plot_N_map(N=Report$Z_s[3,(it-1)*Data$n_s+c(1:Data$n_s)],Grid=Grid) 
  plot_N_map(N=Report$Z_s[4,(it-1)*Data$n_s+c(1:Data$n_s)],Grid=Grid) 
  
  #plot counts
  isp=3
  Count = rep(NA,nrow(Grid))
  Tmp = Data$S_i[1:nrow(Data$C_i)]+1  #recall Data$S_i defined for C++
  #Count[Tmp %% Data$n_s]=Data$C_i[,isp]/Data$P_i
  Count[Tmp %% Data$n_s]=Data$C_i[,isp]/(Data$P_i/(Data$A_s[Data$S_i+1]*Data$Ice_s[Data$S_i+1]))
  plot_N_map(N=Count,Grid=Grid) 

  #plot counts as function of ice
  S_i_tmp = (Data$S_i+1)[1:258] #take out pseudozeroes
  plot(Data$X_s[,"ice_conc"][S_i_tmp],Data$C_i[,1],main="Spotted")
  plot(Data$X_s[,"ice_conc"][S_i_tmp],Data$C_i[,2],main="Ribbon")
  plot(Data$X_s[,"ice_conc"][S_i_tmp],Data$C_i[,3],main="Bearded")
  plot(Data$X_s[,"ice_conc"][S_i_tmp],Data$C_i[,4],main="Ringed")
  
  #plot ice
  it=16
  plot_N_map(N=Data$X_s[(it-1)*Data$n_s+c(1:Data$n_s),"ice_conc"],Grid=Grid) 
  
  it=16
  plot_N_map(N=X_s[(it-1)*Data$n_s+c(1:Data$n_s),"fi_nic"],Grid=Grid) 
  
  Ice = c(0:100)/100
  plot(Ice*Report$beta_bd[1]+Ice^2*Report$beta_bd[2],ylim=c(-4,4))
  plot(Ice*Report$beta_sd[1]+Ice^2*Report$beta_sd[2],ylim=c(-4,4))
  
  
  lines(Ice*Report$Beta[3,4]+Ice^2*Report$Beta[3,7])
  
  plot(Ice*Report$Beta[1,4]+Ice^2*Report$Beta[1,7]+Ice^3*Report$Beta[1,16],ylim=c(-4,4))
  lines(Ice*Report$Beta[3,4]+Ice^2*Report$Beta[3,7]+Ice^2*Report$Beta[3,16])
  
 
  plot(Ice*Report$Beta[2,4]+Ice^2*Report$Beta[2,7])
  lines(Ice*Report$Beta[4,4]+Ice^2*Report$Beta[4,7])
  
  
 
  
  Converge=Opt$convergence
  # SD
  if(Converge==0){
    Report = Obj$report()
    #if( all(c("Etainput_s")%in%names(Map)) ){
    #SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )
    #SD$unbiased$value = c("total_abundance"=Report$total_abundance)
    #}else{
    SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=TRUE )
    #}
    Opt[["run_time"]] = Sys.time()-Start_time
  }
  
  #tabulate counts of each species by strata
  Strata_ID = as.numeric(as.character(X_s[Data$S_i,"Strata"]))[1:316]
  Counts_sum = matrix(0,3,5)
  for(i in 1:3){
    Cs = Data$C_i[Strata_ID==i,]
    Counts_sum[i,]=colSums(Cs)
  }
  Counts_sum
  
  
  #randomized quantile residuals
  library(tweedie)
  Resids = 0*cbind(Report$E_count_obs,Report$E_count_pup,Report$E_count_nophoto)
  SMALL=0.000001
  for(icol in 1:5){
    Resids[,icol] = ptweedie(Data$C_i[,icol]-1,mu=Report$E_count_obs[,icol]+SMALL,phi=Report$phi[icol],power=Report$power[icol])+runif(nrow(Data$C_i))*dtweedie(Data$C_i[,icol],mu=Report$E_count_obs[,icol]+SMALL,phi=Report$phi[icol],power=Report$power[icol])
  }
  Resids[,6] = ptweedie(Data$Pups_i-1,mu=Report$E_count_pup+SMALL,phi=Report$phi[6],power=Report$power[6])+runif(length(Data$Pups_i))*dtweedie(Data$Pups_i,mu=Report$E_count_pup+SMALL,phi=Report$phi[6],power=Report$power[6])
  Resids[,7] = ptweedie(Data$Nophoto_i-1,mu=Report$E_count_nophoto+SMALL,phi=Report$phi[7],power=Report$power[7])+runif(length(Data$Nophoto_i))*dtweedie(Data$Nophoto_i,mu=Report$E_count_nophoto+SMALL,phi=Report$phi[7],power=Report$power[7])
  
  
  Resids[Resids>1]=0.999
  
  Resid_binned = matrix(0,20,7)
  for(irow in 1:nrow(Resids)){
    for(icol in 1:7){
      Resid_binned[ceiling(Resids[irow,icol]/0.05),icol]=Resid_binned[ceiling(Resids[irow,icol]/0.05),icol]+1
    }
  }
  Xsq = rep(0,7)
  for(i in 1:7){
    Xsq[i]=20/nrow(Resids)*sum((Resid_binned[,i]-nrow(Resids)/20)^2)
  }
  Pval = 1-pchisq(Xsq,df=19)   #chi square p-value for each bin
  
  Resids.df = data.frame(Residual = as.vector(Resids))
  Labels1 = rep('',7)
  for(i in 1:7){
    Labels1[i] = paste0('Obs = ',i,', p=',format(Pval[i],digits=2))
  }
  Resids.df$Labels=rep(Labels1,each=nrow(Resids))
  
  ggplot(Resids.df)+geom_histogram(aes(x=Residual),breaks=c(0:20)/20)+facet_wrap(~Labels)
  

