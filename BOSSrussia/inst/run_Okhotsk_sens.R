


library( RandomFields )
library( TMB )
library( INLA )
#library(TMBdebug)

logit<-function(x) log(x/(1-x))

# Compile
setwd("C:/Users/paul.conn/git/BOSSrussia/")
TmbFile = "C:/Users/paul.conn/git/BOSSrussia/BOSSrussia/src/OkhotskTMB"
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 


#require(STabundance)
source('./BOSSrussia/R/util_funcs.R')
load("./Okhotsk_Grid_2013.rdat")
Grid=Data$Grid  #for plotting
load('Okhotsk_TMB_data.RData')  #created with format_Okhotsk_TMB.R script in OkhotskST project



# Init Parameter values
n_species=4
Data$X_s=Data$X_s[,-which(colnames(Data$X_s)%in%c("fi"))]
Data$X_s = Data$X_s[,-which(colnames(Data$X_s)%in%c("easting","northing"))]  #take out easting, northing since we don't want to extrapolate too much into unsurveyed areas
#Data$X_s$fi2 = Data$X_s$fi^2  #add quadratic for landfast ice
n_beta = ncol(Data$X_s)
X_s = Data$X_s
Data$X_s = as.matrix(Data$X_s)
Beta_VC1 = matrix(100 * solve(crossprod(Data$X_s)),n_beta,n_beta)
Data$Beta_VC = kronecker(diag(n_species),Beta_VC1)
Data$C_i= as.matrix(Data$C_i)
Data$phi_max = 10
Data$Ice_s = Data$X_s[,"ice"]

X_sens = vector("list",5)
X_sens[[1]]=Data$X_s[,-c(1,2)]
X_sens[[2]]=Data$X_s[,-c(4,5)]
X_sens[[3]]=Data$X_s[,-6]
X_sens[[4]]=Data$X_s[,-7]
X_sens[[5]]=Data$X_s[,-8]

Etainput_s = matrix(0,n_species,nrow(Data$spde$M0))  #spatial random effects

  Params = list("beta"=matrix(0*rnorm(n_species*n_beta),n_species,n_beta), 
                "logtau_z"=rep(0,n_species), "logkappa_z"=rep(-.9,n_species),
                "phi_logit"=rep(0.0,n_species+3),"p_logit"=rep(0,n_species+3),"thin_logit_i"=Data$thin_mu_logit,"MisID_pars"=Data$MisID_mu,
                "logit_Pup_prop"=Data$Pup_prop_mu, "Etainput_s"=Etainput_s)
  #Params$beta[,4] = -5  #set intercept to low number to start out with
  Params$beta = as.numeric(t(Params$beta))
  
  MisID_orig = Params$MisID_pars
  Params$MisID_pars = rep(-10,length(Params$MisID_pars))

  # Random
  #Random = c( "Etainput_s","thin_logit_i","MisID_pars","logit_Pup_prop" )
  #Random = c( "Etainput_s","thin_logit_i","logit_Pup_prop" )
  #Random = c("Etainput_s")
  #Random = c("beta","thin_logit_i","logit_Pup_prop")
  Random = c("beta","logit_Pup_prop")
  #Random = NULL
   #Random = list()
  
  # Fix parameters
  Map = list()
  #Map[["phi_log"]]=factor(rep(1,n_species+3))
  Map[["phi_logit"]]=factor(c(1,2,3,4,1,5,6))
  Map[["p_logit"]]=factor(c(7,8,9,10,7,11,12))
  #Map[["phi_log"]]=factor(rep(2,n_species+3))
 # Map[["p_logit"]]=factor(rep(1,n_species+3))
  #Map[["logkappa_z"]] = factor(rep(3,n_species))
  #Map[["logtau_z"]] = factor(rep(4,n_species))
  Map[["logkappa_z"]] = factor(rep(NA,n_species))
  Map[["logtau_z"]] = factor(rep(NA,n_species))
  Map[["MisID_pars"]]=factor(rep(NA,length(Data$MisID_mu)))
  Map[["thin_logit_i"]]=factor(rep(NA,length(Data$thin_mu_logit)))
  Map[["Etainput_s"]]=factor(rep(NA,length(Params$Etainput_s)))
  #Map[["logit_Pup_prop"]]=factor(rep(NA,length(Params$logit_Pup_prop)))
  
  #first run: no misID
  
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
  
  
  Params$phi_logit = logit(Report$phi/Data$phi_max)
  Params$p_logit = logit(Report$power-1)
  Params$MisID_pars = MisID_orig
  Params$logit_Pup_prop = logit(Report$Pup_prop)
  Params$beta = as.numeric(t(Report$Beta))
  

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
  
  
  n_sens=length(X_sens)
  Extrap_num = Extrap_ratio = Max_obs = N = SE = matrix(NA,4,n_sens) 
  Converge= LogL=rep(0,n_sens)
  
  
  #now sensitivity
  for(isens in 1:n_sens){
    Data$X_s = X_sens[[isens]]
    n_beta= ncol(Data$X_s)
    Params$beta = rep(0,4*ncol(Data$X_s))
    Beta_VC1 = matrix(100 * solve(crossprod(Data$X_s)),n_beta,n_beta)
    Data$Beta_VC = kronecker(diag(n_species),Beta_VC1)
    
    Params$MisID_pars = rep(-10,length(Params$MisID_pars))
    
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
    
    Params$phi_logit = logit(Report$phi/Data$phi_max)
    Params$p_logit = logit(Report$power-1)
    Params$MisID_pars = MisID_orig
    Params$logit_Pup_prop = logit(Report$Pup_prop)
    Params$beta = as.numeric(t(Report$Beta))
    
    if(Opt$convergence==0 & Opt$message!="X_convergence (3)"){
      SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )
      #if(SD$pdHess==TRUE){
        Converge[isens]=1
        LogL[isens]=Report$jnll
        N[,isens]=Report$total_abundance
        SE[,isens]=SD$sd[which(names(SD$value)=="total_abundance")]
        for(isp in 1:4){
          Max_obs[isp,isens] = max(Report$Z_adj[isp,Data$S_i])
          Extrap_num[isp,isens] = sum(Report$Z_adj[isp,] > Max_obs[isp,isens])
          Extrap_ratio[isp,isens]=max(Report$Z_adj[isp,])/Max_obs[isp,isens]
        }
        save.image('sens_out_Okhotsk.RData')
      #}
    }
  }
  
  
