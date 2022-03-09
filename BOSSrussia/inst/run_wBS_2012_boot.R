


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
X_s$fi2=X_s$fi_nic^2

### Use Irina's models from Chernook et al. 2018
X_s$Strata=factor(X_s$Strata)
X_bd = model.matrix(~Strata*(ice_conc+ice2+dist_mainland+dist_shelf+dist_contour),X_s)
X_rd = model.matrix(~Strata*(ice_conc+ice2+dist_mainland+dist_shelf+dist_contour),X_s)
X_bd = X_bd[,-1]  #get rid of intercept
X_rd = X_rd[,-1]  

X_dm = model.matrix(~Strata+ice_conc+ice2+dist_shelf,X_s)  #ribbon, spotted
X_rn = X_sd = X_dm[,-1]


Data$X_sd = as.matrix(X_sd)
Data$X_rn = as.matrix(X_rn)
Data$X_bd = as.matrix(X_bd)
Data$X_rd = as.matrix(X_rd)
#NOTE : I'm not currently putting priors on regression coefficients.  They don't seem to work well
# with stratum covariates
Data$Beta_VCsd = matrix(100 * solve(crossprod(Data$X_sd)),ncol(Data$X_sd),ncol(Data$X_sd)) #
Data$Beta_VCrn = matrix(100 * solve(crossprod(Data$X_rn)),ncol(Data$X_rn),ncol(Data$X_rn))
Data$Beta_VCbd = matrix(100 * solve(crossprod(Data$X_bd)),ncol(Data$X_bd),ncol(Data$X_bd))
Data$Beta_VCrd = matrix(100 * solve(crossprod(Data$X_rd)),ncol(Data$X_rd),ncol(Data$X_rd))
Data$phi_max = 10

Data$Melt_i = Data$S_i %/% Data$n_s -15 #just set = to day of survey relative to May 5

#Data$Beta_VC = as.matrix(Matrix::bdiag(Beta_VCsd,Beta_VCrn,Beta_VCbd,Beta_VCrd))

Data$C_i= as.matrix(Data$C_i)

#see how sensitive estimates are to cells 216,217
#Data$C_i[216,1]=Data$C_i[216,2]=Data$C_i[217,1]=Data$C_i[217,2]=50
#Data$Nophoto_i[c(216,217)]=50

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


  Random = c("logit_Pup_prop") 

  # Fix parameters
  Map = list()
  Map[["MisID_pars"]]=factor(rep(NA,length(Data$MisID_mu)))
  Map[["thin_beta_day"]]=factor(rep(NA,2*n_species))
  Map[["thin_logit_i"]]=factor(rep(NA,length(Data$thin_mu_logit)))
  Map[["beta0_ringed"]]=Map[["beta1_ringed"]]=factor(rep(NA,1))  #turn off ringed seal lair model
  

  
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
  
  Base = Report
  SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=TRUE )

  
  MisID_orig = Data$MisID_mu
  #now bootstrap thinning, misID params to get SD for N
  logit<-function(x)log(x/(1-x))
  Params$beta_bd = Report$beta_bd
  Params$beta_rd = Report$beta_rd
  Params$beta_rn = Report$beta_rn
  Params$beta_sd = Report$beta_sd
  Params$phi_logit = logit(Report$phi/Data$phi_max)
  Params$p_logit = logit(Report$power-1)
  Params$MisID_pars = MisID_orig
  Params$logit_Pup_prop = logit(Report$Pup_prop)
  
  
  Ests = Params
  Ests$N = Report$N
  Ests$Z = Report$Z
  
  
  set.seed(12345)
  n_boot = 500
  N_boot = SD_boot = matrix(0,4,n_boot)
  Sigma_thin = as.matrix(Data$Sigma_logit_thin)
  Sigma_MisID = as.matrix(Data$MisID_Sigma)
  Converge = rep(0,n_boot)
  
  for(iboot in 1:n_boot){
    cat(paste0("iboot = ",iboot,"\n"))
    #Data$p = rbeta(1,17,8)
    Params$thin_logit_i=mgcv::rmvn(1,Data$thin_mu_logit,Sigma_thin)
    Params$MisID_pars=mgcv::rmvn(1,Data$MisID_mu,Sigma_MisID)
    
    Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
    Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
    Upper = 50
    Opt = nlminb( start=Obj$par, objective=Obj$fn, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))         #
    Report = Obj$report()
    

    if(Opt$convergence==0 & Opt$message!="X_convergence (3)"){
      Converge[iboot]=1
      N_boot[,iboot]=Report$N
      save.image('boot_out_wBS2012.RData')
    }
  }
  
  SE.sd = sqrt(SD$sd[1]^2+var(N_boot[1,1:n_boot]))
  SE.rn = sqrt(SD$sd[2]^2+var(N_boot[2,1:n_boot]))
  SE.bd = sqrt(SD$sd[3]^2+var(N_boot[3,1:n_boot]))
  SE.rd = sqrt(SD$sd[4]^2+var(N_boot[4,1:n_boot]))
  
  save.image('boot_out_wBS2012.RData')