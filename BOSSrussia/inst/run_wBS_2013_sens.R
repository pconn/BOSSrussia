


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

X_sens = vector("list",3)  #strata only, mainland, contour
X_strata = model.matrix(~Strata,X_s)
X_strata = X_strata[,-1]
X_sens[[1]]=list("bd" = X_strata,
                 "rn" = X_strata,
                 "rd" = X_strata,
                 "sd" = X_strata)
X_sens[[2]]=list("bd" = Data$X_bd[,-grep("dist_mainland",colnames(X_bd))],
                 "rn" = Data$X_rn,
                 "rd" = Data$X_rd[,-grep("dist_mainland",colnames(X_rd))],
                 "sd" = Data$X_sd)
X_sens[[3]]=list("bd" = Data$X_bd[,-grep("dist_contour",colnames(X_bd))],
                 "rn" = Data$X_rn,
                 "rd" = Data$X_rd[,-grep("dist_contour",colnames(X_rd))],
                 "sd" = Data$X_sd)


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
  Base = Report
  Converge=Opt$convergence
  # SD
  if(Converge==0){
    SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=TRUE )
    #}
    Opt[["run_time"]] = Sys.time()-Start_time
  }
  
  logit<-function(x)log(x)/log(1-x)
  
  Params$phi_logit = logit(Report$phi/Data$phi_max)
  Params$p_logit = logit(Report$power-1)
  #Params$MisID_pars = MisID_orig
  Params$logit_Pup_prop = logit(Report$Pup_prop)
  
  
  n_sens=length(X_sens)
  Extrap_num = Extrap_ratio = Max_obs = N = SE = matrix(NA,4,n_sens) 
  Converge= LogL=rep(0,n_sens)
  
  
  #now sensitivity
  for(isens in 1:n_sens){  #note haven't modified Beta_VCs because not using priors in analysis
    Data$X_bd = X_sens[[isens]]$bd
    Data$X_rd = X_sens[[isens]]$rd
    Data$X_sd = X_sens[[isens]]$sd
    Data$X_rn = X_sens[[isens]]$rn
    
    Params$beta_sd=rep(0,ncol(Data$X_sd))
    Params$beta_rn=rep(0,ncol(Data$X_rn))
    Params$beta_bd=rep(0,ncol(Data$X_bd))
    Params$beta_rd=rep(0,ncol(Data$X_rd))
    
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
    
    if(Opt$convergence==0 & Opt$message!="X_convergence (3)"){
      SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )
      #if(SD$pdHess==TRUE){
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
      #}
    }
  }
  
  #note: still extrapolation in 'strata' model since ice changes and Z_s is a function of ice AND total avail habitat
  
  #"base model' extrap values
  Extrap_num_base = Extrap_ratio_base =rep(0,4)
  for(isp in 1:4){
    Max_obs = max(Report$Z_s[isp,Data$S_i])
    Extrap_num_base[isp] = sum(Base$Z_s[isp,] > Max_obs)
    Extrap_ratio_base[isp]=max(Base$Z_s[isp,])/Max_obs
  }
  Extrap_num_base
  Extrap_ratio_base
  
  
  