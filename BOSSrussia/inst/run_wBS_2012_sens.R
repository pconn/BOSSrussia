


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
load('wBS_2012_Data_TMB.RData')  #created with format_Okhotsk_TMB.R script in OkhotskST project
Data = Data_2012
rm(Data_2012)


# Init Parameter values
n_species=4

X_s = Data$X_s

### Use Irina's models from Chernook et al. 2018
#  Note, Strata 1 = Middle, Strata 2 = North, Strata 3 = South
X_s$Strata=factor(X_s$Strata)
X_bd = model.matrix(~Strata*(ice_conc+ice2+dist_mainland+dist_shelf+dist_contour),X_s)
X_rd = model.matrix(~Strata*(ice_conc+ice2+dist_mainland+dist_shelf+dist_contour),X_s)
X_bd = X_bd[,-1]  #get rid of intercept
X_rd = X_rd[,-1]  

X_dm = model.matrix(~Strata+ice_conc+ice2+dist_shelf,X_s)
X_rn = X_sd = X_dm[,-1]


Data$X_sd = as.matrix(X_sd)
Data$X_rn = as.matrix(X_rn)
Data$X_bd = as.matrix(X_bd)
Data$X_rd = as.matrix(X_rd)


### sensitivity now by species;  ribbon/bearded -Strata, -dist_mainland, -dist_shelf, -dist_contour
Sens_mods = vector("list",4) #holds design matrix for each species & sensitivity run (dim n_species x n_sens_runs)
for(isp in 1:4){
  Sens_mods[[isp]] = vector("list",16)
}

#spotted, ribbon seals
#start by inserting full models for all [to be overwritten where appropriate]
for(imod in 1:16)Sens_mods[[1]][[imod]]=Sens_mods[[2]][[imod]]=X_sd
# minus strata
X_tmp = model.matrix(~ice_conc+ice2+dist_shelf,X_s)
Sens_mods[[1]][[1]] = as.matrix(X_tmp[,-1])
Sens_mods[[2]][[5]] = as.matrix(X_tmp[,-1])
# minus shelf
X_tmp = model.matrix(~Strata+ice_conc+ice2,X_s)
Sens_mods[[1]][[2]] = as.matrix(X_tmp[,-1])
Sens_mods[[2]][[6]] = as.matrix(X_tmp[,-1])
# plus contour
X_tmp = model.matrix(~Strata+ice_conc+ice2+dist_shelf+dist_contour,X_s)
Sens_mods[[1]][[3]] = as.matrix(X_tmp[,-1])
Sens_mods[[2]][[7]] = as.matrix(X_tmp[,-1])
# plus mainland
X_tmp = model.matrix(~Strata+ice_conc+ice2+dist_shelf+dist_mainland,X_s)
Sens_mods[[1]][[4]] = as.matrix(X_tmp[,-1])
Sens_mods[[2]][[8]] = as.matrix(X_tmp[,-1])

## bearded, ringed seals
#start by inserting full models for all [to be overwritten where appropriate]
for(imod in 1:16)Sens_mods[[3]][[imod]]=Sens_mods[[4]][[imod]]=X_bd
#minus strata
X_tmp = model.matrix(~ice_conc+ice2+dist_mainland+dist_shelf+dist_contour,X_s)
Sens_mods[[3]][[9]] = as.matrix(X_tmp[,-1])
Sens_mods[[4]][[13]] = as.matrix(X_tmp[,-1])
#minus mainland
X_tmp = model.matrix(~Strata*(ice_conc+ice2+dist_shelf+dist_contour),X_s)
Sens_mods[[3]][[10]] = as.matrix(X_tmp[,-1])
Sens_mods[[4]][[14]] = as.matrix(X_tmp[,-1])
#minus shelf
X_tmp = model.matrix(~Strata*(ice_conc+ice2+dist_mainland+dist_contour),X_s)
Sens_mods[[3]][[11]] = as.matrix(X_tmp[,-1])
Sens_mods[[4]][[15]] = as.matrix(X_tmp[,-1])
#minus contour
X_tmp = model.matrix(~Strata*(ice_conc+ice2+dist_shelf+dist_mainland),X_s)
Sens_mods[[3]][[12]] = as.matrix(X_tmp[,-1])
Sens_mods[[4]][[16]] = as.matrix(X_tmp[,-1])

#NOTE : I'm not currently putting priors on regression coefficients.  They don't seem to work well
# with stratum covariates
Data$Beta_VCsd = matrix(100 * solve(crossprod(Data$X_sd)),ncol(Data$X_sd),ncol(Data$X_sd)) #
Data$Beta_VCrn = matrix(100 * solve(crossprod(Data$X_rn)),ncol(Data$X_rn),ncol(Data$X_rn))
Data$Beta_VCbd = matrix(100 * solve(crossprod(Data$X_bd)),ncol(Data$X_bd),ncol(Data$X_bd))
Data$Beta_VCrd = matrix(100 * solve(crossprod(Data$X_rd)),ncol(Data$X_rd),ncol(Data$X_rd))
Data$phi_max = 10

Data$Melt_i = Data$S_i %/% Data$n_s -15 #just set = to day of survey relative to May 5

Data$C_i= as.matrix(Data$C_i)

ice_thresh = 0.01  
Data$Ice_s = Data$X_s[,"ice_conc"]
Data$Which_zero = which(Data$X_s[,"ice_conc"]<0.01)-1  


Params = list("log_N" = rep(log(50000),n_species),"beta_sd"=rep(0,ncol(Data$X_sd)),"beta_rn"=rep(0,ncol(Data$X_rn)),
              "beta_bd"=rep(0,ncol(Data$X_bd)),"beta_rd"=rep(0,ncol(Data$X_rd)),
              "thin_beta_day" = rep(0,2*n_species),  
              "phi_logit"=rep(0,n_species+3),"p_logit"=rep(0,n_species+3),
                "thin_logit_i"=Data$thin_mu_logit,"MisID_pars"=Data$MisID_mu,
                "logit_Pup_prop"=Data$Pup_prop_mu,"beta0_ringed"=10,"beta1_ringed"=0)
Params$thin_logit_i = Data$thin_mu_logit_nolair_chukchi
Data$thin_mu_logit = Data$thin_mu_logit_nolair_chukchi
Data$Sigma_logit_thin = Data$Sigma_logit_thin_nolair_chukchi
Data$h_mean = 0 #not used currently

Data_orig=Data
Params_orig=Params

  Random = c("logit_Pup_prop") 

  # Fix parameters
  Map = list()
  Map[["MisID_pars"]]=factor(rep(NA,length(Data$MisID_mu)))
  Map[["thin_beta_day"]]=factor(rep(NA,2*n_species))
  Map[["thin_logit_i"]]=factor(rep(NA,length(Data$thin_mu_logit)))
  Map[["beta0_ringed"]]=Map[["beta1_ringed"]]=factor(rep(NA,1))  #turn off ringed seal lair model
  
  # Make object
  dyn.load( dynlib(TmbFile) )
  Start_time = Sys.time()
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
  Obj$fn( Obj$par )
  
  # Run
  Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
  Upper = 50
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max = 10000,maxit=10000))         #
  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  Report = Obj$report()
  
  Base = Report
  SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=TRUE )

  
  MisID_orig = Data$MisID_mu
  logit<-function(x)log(x/(1-x))
  Params$phi_logit = logit(Report$phi/Data$phi_max)
  Params$p_logit = logit(Report$power-1)
  Params$MisID_pars = MisID_orig
  Params$logit_Pup_prop = logit(Report$Pup_prop)
  
  
  n_sens=16
  Extrap_num = Extrap_ratio = Max_obs = N = SE = matrix(NA,4,n_sens) 
  Converge= LogL=rep(0,n_sens)
  
  
  #now sensitivity
  for(isens in 1:n_sens){  #note haven't modified Beta_VCs because not using priors in analysis

    Data$X_sd = Sens_mods[[1]][[isens]] 
    Data$X_rn = Sens_mods[[2]][[isens]] 
    Data$X_bd = Sens_mods[[3]][[isens]] 
    Data$X_rd = Sens_mods[[4]][[isens]] 
    

    Params$beta_sd=rep(0,ncol(Data$X_sd))
    Params$beta_rn=rep(0,ncol(Data$X_rn))
    Params$beta_bd=rep(0,ncol(Data$X_bd))
    Params$beta_rd=rep(0,ncol(Data$X_rd))
 
    dyn.load( dynlib(TmbFile) )
    Start_time = Sys.time()
    Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
    Obj$fn( Obj$par )
    
    # Run
    Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
    Upper = 50
    Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max = 10000,maxit=10000))         #
    Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
    Report = Obj$report()  
    
    if(Opt$convergence==0 & Opt$message!="X_convergence (3)"){
      SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )
      Converge[isens]=1
      LogL[isens]=Report$jnll
      N[,isens]=Report$N
      SE[,isens]=SD$sd[which(names(SD$value)=="N")]
      for(isp in 1:4){
        Max_obs[isp,isens] = max(Report$Z_s[isp,Data$S_i])
        Extrap_num[isp,isens] = sum(Report$Z_s[isp,] > Max_obs[isp,isens])
        Extrap_ratio[isp,isens]=max(Report$Z_s[isp,])/Max_obs[isp,isens]
      }
      save.image('sens_out_wBS2012.RData')
    }
  }
  
 
  ### try no misID, alternative ringed seal availability analyses
  
  #no misID
  Data=Data_orig
  Params=Params_orig
  Params$MisID_pars = rep(-10,length(Params$MisID_pars))
  Params$MisID_pars[c(4,8,12,16)]=logit(c(.08,.01,.06,.01))  #maintain p(unknown) for each species at prior values
  dyn.load( dynlib(TmbFile) )
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
  Obj$fn( Obj$par )
  Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
  Upper = 50
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max = 10000,maxit=10000))         #
  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  Report = Obj$report() 
  ### Does not converge!
  
  #run with chukchi sea availability model that includes lair use by ringed seals
  Params$MisID_pars = Data$MisID_mu
  Params$thin_logit_i = Data$thin_mu_logit_lair_chukchi
  Data$thin_mu_logit = Data$thin_mu_logit_lair_chukchi
  Data$Sigma_logit_thin = Data$Sigma_logit_thin_lair_chukchi
  dyn.load( dynlib(TmbFile) )
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
  Obj$fn( Obj$par )
  Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
  Upper = 50
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max = 10000,maxit=10000))         #
  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  Report = Report_chukchi_lair = Obj$report() 

  #run w/ bering sea availability model
  Params$thin_logit_i = Data$thin_mu_logit_nolair_bering
  Data$thin_mu_logit = Data$thin_mu_logit_nolair_bering
  Data$Sigma_logit_thin = Data$Sigma_logit_thin_nolair_bering
  dyn.load( dynlib(TmbFile) )
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
  Obj$fn( Obj$par )
  Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
  Upper = 50
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max = 10000,maxit=10000))         #
  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  Report = Report_bering = Obj$report() 
  
  save.image('sens_out_wBS2012.RData')
  
  
  