


library( RandomFields )
library( TMB )
library( INLA )
library(sf)
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
load('c:/users/paul.conn/git/okhotskST/Okhotsk_TMB_data_Jan2022.RData')


# Init Parameter values
n_species=4
Data$X_s=Data$X_s[,-which(colnames(Data$X_s)%in%c("fi"))]
Data$X_s = Data$X_s[,-which(colnames(Data$X_s)%in%c("easting","northing"))]  #take out easting, northing since we don't want to extrapolate too much into unsurveyed areas
#Data$X_s$fi2 = Data$X_s$fi^2  #add quadratic for landfast ice
#Data$X_s[,"ln_depth"]=log(Data$X_s[,"depth"]^3+1)
#Data$X_s=Data$X_s[,-which(colnames(Data$X_s)%in%c("depth","depth2"))]
n_beta = ncol(Data$X_s)
X_s = Data$X_s
Data$X_s = as.matrix(Data$X_s)
Beta_VC1 = matrix(100 * solve(crossprod(Data$X_s)),n_beta,n_beta)
Data$Beta_VC = kronecker(diag(n_species),Beta_VC1)
Data$C_i= as.matrix(Data$C_i)
Data$phi_max = 10
Data$Ice_s = Data$X_s[,"ice"]
Data$beta_prior=1.0

Etainput_s = matrix(0,n_species,nrow(Data$spde$M0))  #spatial random effects

#cells to includ in bearded estimate that omit northeast corner where extrapolations are evident
Cur.df=data.frame((st_coordinates(st_centroid(st_as_sf(Grid),byid=TRUE))))  
Cur.df2 = Cur.df[Cur.df$X<(-1030000),]
Cells_reduced = which(Cur.df$X<(-1030000))


#### "chukchi" ringed haul-out model - MDSDA only (no julian day)
Data$thin_mu_logit=Data$thin_mu_logit_chukchi
Data$Sigma_logit_thin = Data$Sigma_logit_thin_chukchi
Data$Cells_reduced_s = Cells_reduced-1
Params = list("beta"=matrix(0*rnorm(n_species*n_beta),n_species,n_beta), 
              "logtau_z"=rep(0,n_species), "logkappa_z"=rep(-.9,n_species),
              "phi_logit"=rep(0.0,n_species+3),"p_logit"=rep(0,n_species+3),"thin_logit_i"=Data$thin_mu_logit_chukchi,"MisID_pars"=Data$MisID_mu,
              "logit_Pup_prop"=Data$Pup_prop_mu, "Etainput_s"=Etainput_s)
Params$beta[,4] = -5  #set intercept to low number to start out with
Params$beta = as.numeric(t(Params$beta))

MisID_orig = Params$MisID_pars
Params$MisID_pars = rep(-10,length(Params$MisID_pars))

# Random
#Random = c( "Etainput_s","thin_logit_i","MisID_pars","logit_Pup_prop" )
#Random = c( "Etainput_s","thin_logit_i","logit_Pup_prop" )
#Random = c("Etainput_s")
#Random = c("beta","thin_logit_i","logit_Pup_prop")
Random = c("beta","logit_Pup_prop")
#Random = "logit_Pup_prop"
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
Map[["thin_logit_i"]]=factor(rep(NA,length(Data$thin_mu_logit_chukchi)))
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


#now rerun with misID 
logit <- function(x) log(x/(1-x))
Params$beta = as.numeric(t(Report$Beta))
Params$phi_logit = logit(Report$phi/Data$phi_max)
Params$p_logit = logit(Report$power-1)
Params$MisID_pars = MisID_orig
Params$logit_Pup_prop = logit(Report$Pup_prop)

# Make object
#compile( paste0(Version,".cpp") )
#dyn.load( dynlib(TmbFile) )
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

Converge=Opt$convergence
# SD
if(Converge==0){
  Report = Obj$report()
  SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=TRUE )
  Opt[["run_time"]] = Sys.time()-Start_time
}

  


#now bootstrap thinning, misID params to get SD for N
Params$beta = as.numeric(t(Report$Beta))
Params$phi_logit = logit(Report$phi/Data$phi_max)
Params$p_logit = logit(Report$power-1)
Params$MisID_pars = MisID_orig
Params$logit_Pup_prop = logit(Report$Pup_prop)


Ests = Params
Ests$N = Report$total_abundance
Ests$Z = Report$Z_adj


seed.base=10000
n_boot = 150
N_boot = SD_boot = matrix(0,4,n_boot)
Sigma_thin = as.matrix(Data$Sigma_logit_thin)
Sigma_MisID = as.matrix(Data$MisID_Sigma)
Converge = rep(0,n_boot)

for(iboot in 1:n_boot){
  set.seed(seed.base+iboot)
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
    N_boot[,iboot]=Report$total_abundance
    #length_sd = length(SD$sd)
    #SD_boot[,iboot]=SD$sd[(length_sd-1):length_sd]
    save.image('boot_out_Okhotsk_1.RData')
  }
  #}
}

#SE.bd = sqrt(mean((SD_boot[1,1:1000])^2,na.rm=TRUE)+var(N_boot[1,1:1000]))
#SE.rd = sqrt(mean((SD_boot[2,1:1000])^2,na.rm=TRUE)+var(N_boot[2,1:1000]))
Which_converge = which(Converge==1)
SE.sd = sqrt(SD$sd[1]^2+var(N_boot[1,Which_converge]))
SE.rn = sqrt(SD$sd[2]^2+var(N_boot[2,Which_converge]))
SE.bd = sqrt(SD$sd[3]^2+var(N_boot[3,Which_converge]))
SE.rd = sqrt(SD$sd[4]^2+var(N_boot[4,Which_converge]))

Boot_out = list(N_boot=N_boot,Converge=Converge)
save(Boot_out,'boot_out_Okhotsk_Apr2021_1.RData')
