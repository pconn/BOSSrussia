


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
X_sd = cbind(X_s$ice_conc,X_s$ice2,X_s$depth,X_s$depth2,X_s$sqrt_edge,X_s$dist_edge,X_s$dist_shelf,X_s$dist_contour)   #fixed effects (non-spline) 
X_rn = cbind(X_s$ice_conc,X_s$ice2,X_s$depth,X_s$depth2,X_s$sqrt_edge,X_s$dist_edge,X_s$dist_shelf,X_s$dist_contour)   #fixed effects (non-spline) 
X_rd = cbind(X_s$ice_conc,X_s$ice2,X_s$fi_nic,X_s$dist_contour,X_s$dist_mainland,X_s$depth,X_s$depth2,X_s$northing,X_s$easting,X_s$northing2,X_s$easting2)
X_bd = cbind(X_s$fi_nic,X_s$dist_contour,X_s$sqrt_edge,X_s$depth,X_s$depth2,X_s$northing,X_s$easting,X_s$northing2,X_s$easting2)

### Use Irina's models from Chernook et al. 2018
X_s$Strata=factor(X_s$Strata)
X_sd = model.matrix(~Strata*(ice_conc+ice2+dist_mainland+dist_shelf+dist_contour),X_s)
X_sd = X_sd[,-c(1,8,10,12,14,16)]  #remove interactions in "north" strata since no seals observed there
X_rn = model.matrix(~Strata*(ice_conc+ice2+dist_mainland+dist_shelf+dist_contour),X_s)
X_rn = X_rn[,-c(1,8,10,12,14,16)]  #remove interactions in "north" strata since no seals observed there
X_bd = model.matrix(~Strata*(ice_conc+ice2+dist_mainland+dist_shelf+dist_contour),X_s)
X_rd = model.matrix(~Strata*(ice_conc+ice2+dist_mainland+dist_shelf+dist_contour),X_s)

Which_1 = which(X_s[,"Strata"]==1)
Which_2 = which(X_s[,"Strata"]==2)
Which_3 = which(X_s[,"Strata"]==3)
X_s1 = X_s[Which_1,]
X_s2 = X_s[Which_2,]
X_s3 = X_s[Which_3,]
X_sd3 = X_s3[,colnames(X_s)%in% c("ice_conc","ice2","dist_shelf")]
#X_sd3 = X_s3[,colnames(X_s)%in% c("dist_mainland","ice_conc","dist_contour","ice2","dist_shelf")]
X_sd1 = cbind(1,X_s1[,colnames(X_s)%in% c("ice_conc","ice2","dist_shelf")]) #don't include distance from mainland here as it can lead to some clear extrapolations for ribbon seals

X_sd2 = matrix(1,nrow=nrow(X_s2),1) #just a column of ones here since no observations
#X_sd1 = matrix(1,nrow=nrow(X_s1),1) #just a column of ones here since no observations
#X_sd3 = matrix(1,nrow=nrow(X_s3),1) #just a column of ones here since no observations

#now assemble into larger design matrix
X_sd = matrix(0,nrow(X_s),ncol(X_sd1)+ncol(X_sd2)+ncol(X_sd3))
X_sd[Which_1,1:ncol(X_sd1)]=as.matrix(X_sd1)
X_sd[Which_2,ncol(X_sd1)+1]=as.matrix(X_sd2)
X_sd[Which_3,(ncol(X_sd1)+2):(ncol(X_sd1)+1+ncol(X_sd3))]=as.matrix(X_sd3)

#take off ice concentration for southern strata for ribbon, spotted (decreasing function)
#X_sd = X_sd[,-c(3,5)]
X_rn = X_sd
X_bd = X_bd[,-1] #remove intercept
X_rd = X_rd[,-1]


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
  #Map[["beta_sd"]]=factor(c(NA,1,2))
  #Map[["beta_rn"]]=factor(c(NA,1,2))
  
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
  n_boot = 10
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
    
    # if(Opt$convergence==0 & Opt$message!="X_convergence (3)"){
    #   
    #   #now turn ringed seal melt on
    #   Params$log_N = Opt$par[1:2]
    #   Params$phi_log_rus =Opt$par[3:5]
    #   Params$p_logit_rus = Opt$par[6:8]
    #   Params$phi_log_us = Opt$par[9:11]
    #   Params$p_logit_us = Opt$par[12:14]
    #   Params$p_sp_logit = Opt$par[15]
    #   Params$log_lambda = Report$log_lambda
    #   Params$Beta = Report$Beta
    #   Map[["beta0_ringed"]]=NULL
    #   Map[["beta1_ringed"]]=NULL
    #   
    #   Obj = MakeADFun( data=CHESS_data, parameters=Params, random=Random, map=Map, silent=FALSE)
    #   Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
    #   Upper = 50
    #   Opt = try(nlminb( start=Obj$par, objective=Obj$fn, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000)))         #
    #   Report = Obj$report()
    
    if(Opt$convergence==0 & Opt$message!="X_convergence (3)"){
      Converge[iboot]=1
      #SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )
      N_boot[,iboot]=Report$N
      #length_sd = length(SD$sd)
      #SD_boot[,iboot]=SD$sd[(length_sd-1):length_sd]
      save.image('boot_out_wBS2013.RData')
    }
    #}
  }
  
  #SE.bd = sqrt(mean((SD_boot[1,1:1000])^2,na.rm=TRUE)+var(N_boot[1,1:1000]))
  #SE.rd = sqrt(mean((SD_boot[2,1:1000])^2,na.rm=TRUE)+var(N_boot[2,1:1000]))
  
  SE.sd = sqrt(SD$sd[1]^2+var(N_boot[1,1:n_boot]))
  SE.rn = sqrt(SD$sd[2]^2+var(N_boot[2,1:n_boot]))
  SE.bd = sqrt(SD$sd[3]^2+var(N_boot[3,1:n_boot]))
  SE.rd = sqrt(SD$sd[4]^2+var(N_boot[4,1:n_boot]))
  
  save.image('boot_out_wBS2012.RData')