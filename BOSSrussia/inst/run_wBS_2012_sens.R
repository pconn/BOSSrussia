


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
X_rd = X_bd
X_bd = X_bd[,-1]
X_rd = X_rd[,-1]

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
colnames(X_sd)=c("inter1","dist_shelf1","ice1","ice12","inter2","dist_shelf3","ice3","ice32")
X_rn = X_sd

Data$X_sd = as.matrix(X_sd)
Data$X_rn = as.matrix(X_rn)
Data$X_bd = as.matrix(X_bd)
Data$X_rd = as.matrix(X_rd)

X_sens = vector("list",4)  #strata only, mainland, shelf, contour
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
X_sens[[3]]=list("bd" = Data$X_bd[,-grep("dist_shelf",colnames(X_bd))],
                 "rn" = Data$X_rn[,-grep("dist_shelf",colnames(X_rn))],
                 "rd" = Data$X_rd[,-grep("dist_shelf",colnames(X_rd))],
                 "sd" = Data$X_sd[,-grep("dist_shelf",colnames(X_sd))])
X_sens[[4]]=list("bd" = Data$X_bd[,-grep("dist_contour",colnames(X_bd))],
                 "rn" = Data$X_rn,
                 "rd" = Data$X_rd[,-grep("dist_contour",colnames(X_rd))],
                 "sd" = Data$X_sd)

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
  #Params$beta_bd = Report$beta_bd
  #Params$beta_rd = Report$beta_rd
  #Params$beta_rn = Report$beta_rn
  #Params$beta_sd = Report$beta_sd
  Params$phi_logit = logit(Report$phi/Data$phi_max)
  Params$p_logit = logit(Report$power-1)
  Params$MisID_pars = MisID_orig
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
  
  
  