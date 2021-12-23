


library( RandomFields )
library( TMB )
library( INLA )
library(sf)

#library(TMBdebug)


# Compile
setwd("C:/Users/paul.conn/git/BOSSrussia/")
TmbFile = "C:/Users/paul.conn/git/BOSSrussia/BOSSrussia/src/OkhotskTMB"
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 


#require(STabundance)
source('./BOSSrussia/R/util_funcs.R')
load("./Okhotsk_Grid_2013.rdat")
Grid=Data$Grid  #for plotting
load('c:/users/paul.conn/git/okhotskST/Okhotsk_TMB_data_Apr2021.RData')  #created with format_Okhotsk_TMB.R script in OkhotskST project



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
  mod_no_misID = Report
  
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
  mod_base = Report
  
  #XY = coordinates(Grid)
  #rownames(XY) = c(1:nrow(XY))
  #plot_N_map_xy(N=Report$Z_s[1,],XY=XY) 
  
  
  library(ggplot2)
  library(gridExtra)
  plot_N_map<-function(N,Grid,No.plot=NA){
    #require(rgeos)
    Tmp=st_as_sf(Grid)
    Abundance=N
    Cur.df=cbind(data.frame((st_coordinates(st_centroid(Tmp,byid=TRUE)))),Abundance)
    Cur.df[,1]=as.integer(Cur.df[,1])
    Cur.df[,2]=as.integer(Cur.df[,2])
    if(! is.na(No.plot)){
      Cur.df=Cur.df[-No.plot,]
    }
    new.colnames=colnames(Cur.df)
    new.colnames[1:2]=c("Easting","Northing")
    colnames(Cur.df)=new.colnames
    tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
    p1=ggplot(Cur.df)+aes(x=Easting,y=Northing,fill=Abundance)+geom_raster()+tmp.theme
    p1 = p1 + scale_fill_viridis_c()
    p1
  }
  sd_plot = plot_N_map(N=Report$Z_adj[1,],Grid) + ggtitle("D. Spotted")
  rn_plot = plot_N_map(N=Report$Z_adj[2,],Grid) + ggtitle("B. Ribbon")
  bd_plot = plot_N_map(N=Report$Z_adj[3,],Grid) + ggtitle("A. Bearded") 
  rd_plot = plot_N_map(N=Report$Z_adj[4,],Grid) + ggtitle("C. Ringed")
  
  pdf("Okhotsk_N_plot.pdf",height=6,width=8)
   grid.arrange(bd_plot,rn_plot,rd_plot,sd_plot,ncol=2)
  dev.off()
  

  
  
  #plot counts
  isp=3
  Count = rep(NA,nrow(Grid))
  Tmp = Data$S_i[1:nrow(Data$C_i)]+1  #recall Data$S_i defined for C++
  Count[Tmp %% Data$n_s]=Data$C_i[,isp]/Data$P_i
  plot_N_map(N=Count,Grid=Grid) 
  
  #plot fast ice
  #plot_N_map(N=X_s[,'fi'],Grid=Grid) 
  
  plot(Report$E_count_obs[,isp],Data$C_i[,isp])
  
  
  
  
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
  
  #randomized quantile residuals
  library(tweedie)
  set.seed(11111)  
  Resids = 0*cbind(Report$E_count_obs,Report$E_count_pup,Report$E_count_nophoto)
  SMALL=0.000001
  for(icol in 1:5){
    Resids[,icol] = ptweedie(Data$C_i[,icol]-1,mu=Report$E_count_obs[,icol]+SMALL,phi=Report$phi[icol],power=Report$power[icol])+runif(nrow(Data$C_i))*dtweedie(Data$C_i[,icol],mu=Report$E_count_obs[,icol]+SMALL,phi=Report$phi[icol],power=Report$power[icol])
  }
  Resids[,6] = ptweedie(Data$Pups_i-1,mu=Report$E_count_pup+SMALL,phi=Report$phi[6],power=Report$power[6])+runif(length(Data$Pups_i))*dtweedie(Data$Pups_i,mu=Report$E_count_pup+SMALL,phi=Report$phi[6],power=Report$power[6])
  Resids[,7] = ptweedie(Data$Nophoto_i-1,mu=Report$E_count_nophoto+SMALL,phi=Report$phi[7],power=Report$power[7])+runif(length(Data$Nophoto_i))*dtweedie(Data$Nophoto_i,mu=Report$E_count_nophoto+SMALL,phi=Report$phi[7],power=Report$power[7])
  
  
  Resids[Resids>1]=0.999
  
  Resid_binned = matrix(0,10,7)
  for(irow in 1:nrow(Resids)){
    for(icol in 1:7){
      Resid_binned[ceiling(Resids[irow,icol]/0.1),icol]=Resid_binned[ceiling(Resids[irow,icol]/0.1),icol]+1
    }
  }
  Xsq = rep(0,7)
  for(i in 1:7){
    Xsq[i]=10/nrow(Resids)*sum((Resid_binned[,i]-nrow(Resids)/10)^2)
  }
  Pval = 1-pchisq(Xsq,df=9)   #chi square p-value for each bin
  
  Resids.df = data.frame(Residual = as.vector(Resids))
  Labels1 = rep('',7)
  for(i in 1:7){
    Labels1[i] = paste0('Obs = ',i,', p=',format(Pval[i],digits=2))
  }
  Resids.df$Labels=rep(Labels1,each=nrow(Resids))
  
pdf('GOF_Okhotsk.pdf')
  ggplot(Resids.df)+geom_histogram(aes(x=Residual),breaks=c(0:10)/10)+facet_wrap(~Labels)
dev.off()



# scale of model
Depth = c(0:16)
Ice = c(0:100)/100
Land = c(0:45)*10000
Break = c(0:35)
Edge = c(0:32)

Depth_plot = Depth^3 
Break_plot = Break / mean(Break)
Edge_plot = Edge / mean(Edge)
Land_plot = Land/mean(Land)

Effect = c(Depth*Report$Beta[1,1]+Depth^2*Report$Beta[1,2],
           Ice * Report$Beta[1,4] + Ice^2*Report$Beta[1,5],
           Land * Report$Beta[1,6],
           Break * Report$Beta[1,7],
           Edge * Report$Beta[1,8]
)
for(isp in 2:4)Effect = c(Effect,
           c(Depth*Report$Beta[isp,1]+Depth^2*Report$Beta[isp,2],
              Ice * Report$Beta[isp,4] + Ice^2*Report$Beta[isp,5],
              Land * Report$Beta[isp,6],
              Break * Report$Beta[isp,7],
              Edge * Report$Beta[isp,8])                          
            )
n_entries = length(Effect)/4
Plot_df = data.frame(Effect = Effect, 
                     Species = rep(c("Spotted","Ribbon","Bearded","Ringed"),each=n_entries),
                     Covariate = rep(c(rep("Depth",length(Depth)),rep("Ice",length(Ice)),
                                       rep("Dist_land",length(Land)),rep("Dist_shelf",length(Break)),
                                       rep("Dist_edge",length(Edge))),4),
                     Cov_value = rep(c(Depth_plot,Ice,Land_plot,Break_plot,Edge_plot),4)
)
cov_plot = ggplot(Plot_df)+geom_line(aes(x=Cov_value,y=Effect,color=Species))+facet_wrap(~Covariate,scales="free")+
  xlab("Covariate value")  
### Depth is in meters, other distance values scaled to their mean (they are in projected space anyway)
pdf("Cov_effects_Okhotsk.pdf")
  cov_plot
dev.off()

#correlation among effects
cor(Data$X_s[,"dist_land"],Data$X_s[,"dist_edge"])
cor(Data$X_s[,"dist_break"],Data$X_s[,"dist_land"])
cor(Data$X_s[,"dist_break"],Data$X_s[,"dist_edge"])
cor(Data$X_s[,"dist_land"],Data$X_s[,"depth"])
cor(Data$X_s[,"dist_break"],Data$X_s[,"depth"])
cor(Data$X_s[,"dist_edge"],Data$X_s[,"depth"])
cor(Data$X_s[,"dist_land"],Data$X_s[,"ice"])
cor(Data$X_s[,"dist_break"],Data$X_s[,"ice"])
cor(Data$X_s[,"dist_edge"],Data$X_s[,"ice"])
cor(Data$X_s[,"depth"],Data$X_s[,"ice"])
#maximum absolute correlation about 0.54

#log-based CIs
load("c:/users/paul.conn/git/BOSSrussia/boot_out_Okhotsk_1.RData")
N_boot_tot = N_boot[,which(Converge==1)]
load("c:/users/paul.conn/git/BOSSrussia/boot_out_Okhotsk_2.RData")
N_boot_tot=cbind(N_boot_tot,N_boot[,which(Converge==1)])
load("c:/users/paul.conn/git/BOSSrussia/boot_out_Okhotsk_3.RData")
N_boot_tot=cbind(N_boot_tot,N_boot[,which(Converge==1)])
load("c:/users/paul.conn/git/BOSSrussia/boot_out_Okhotsk_4.RData")
N_boot_tot=cbind(N_boot_tot,N_boot[,which(Converge==1)])

#SEs = c(37628,42847,91829,221963) #from run_Okhotsk_boot
SE_fun <- function(x) sqrt(var(x))
SEs = apply(N_boot_tot,1,'SE_fun')
C = exp(1.96 * sqrt(log(1+(SEs/Report$total_abundance)^2)))
CI_lo<- Report$total_abundance/C
CI_hi<- Report$total_abundance*C
CI_lo
CI_hi

### gIVH
SD_Z = matrix(SD$sd[9:length(SD$sd)],4,Data$n_s)
Max_sd_obs = rep(0,4)
Which_out = vector("list",4)
for(i in 1:4){
  Max_sd_obs[i]=max(SD_Z[i,Data$S_i])
  Which_out[[i]]=which(SD_Z[i,]>Max_sd_obs[i])
}
#just 3 predictions deemed "extrapolations"; but these predictions not much bigger than for surveyed cells 


#################################
#      SENSITIVITY RUNS
##################################

#drop habitat terms one-by-one
X_sens = vector("list",5)
X_sens[[1]]=Data$X_s[,-c(1,2)]
X_sens[[2]]=Data$X_s[,-c(4,5)]
X_sens[[3]]=Data$X_s[,-6]
X_sens[[4]]=Data$X_s[,-7]
X_sens[[5]]=Data$X_s[,-8]

n_sens=length(X_sens)
Extrap_num = Extrap_ratio = Max_obs = N = SE = matrix(NA,4,n_sens) 
Converge= LogL=rep(0,n_sens)

Beta_VC_main = Data$Beta_VC
X_s_main = Data$X_s

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
  
  #rerun w misID
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
  
  if(isens==2)Report_no_ice = Report
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


# bering ringed haul-out model (jday + MDSDA)
Data$Beta_VC = Beta_VC_main
Data$X_s = X_s_main 
n_beta=ncol(Data$X_s)

Data$thin_mu_logit=Data$thin_mu_logit_jday_bering
Data$Sigma_logit_thin = Data$Sigma_logit_thin_jday_bering
Data$Cells_reduced_s = Cells_reduced-1
Params = list("beta"=matrix(0*rnorm(n_species*n_beta),n_species,n_beta), 
              "logtau_z"=rep(0,n_species), "logkappa_z"=rep(-.9,n_species),
              "phi_logit"=rep(0.0,n_species+3),"p_logit"=rep(0,n_species+3),"thin_logit_i"=Data$thin_mu_logit_jday_bering,"MisID_pars"=Data$MisID_mu,
              "logit_Pup_prop"=Data$Pup_prop_mu, "Etainput_s"=Etainput_s)
#run without misID to start with
Params$MisID_pars = rep(-10,length(Params$MisID_pars))
Start_time = Sys.time()
Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
Obj$fn( Obj$par )
# Run
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max = 10000,maxit=10000))         #
Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
Report = Obj$report() 

#now run again w/ misID
logit <- function(x) log(x/(1-x))
Params$beta = as.numeric(t(Report$Beta))
Params$phi_logit = logit(Report$phi/Data$phi_max)
Params$p_logit = logit(Report$power-1)
Params$MisID_pars = MisID_orig
Params$logit_Pup_prop = logit(Report$Pup_prop)
Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
#Obj = MakeADFun( data=Data, parameters=Params, map=Map, silent=FALSE)
Obj$fn( Obj$par )
# Run
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max = 10000,maxit=10000))         #
Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
Report = Obj$report() 

mod_bering = Report

Converge = c(Converge,c(0,0))
if(Opt$convergence==0)Converge[6]=1
LogL[6]=Report$jnll
N=cbind(N,Report$total_abundance)
Extrap_num = cbind(Extrap_num,0)
Extrap_ratio = cbind(Extrap_ratio,0)
for(isp in 1:4){
  max_obs = max(Report$Z_adj[isp,Data$S_i])
  Extrap_num[isp,6] = sum(Report$Z_adj[isp,] > max_obs)
  Extrap_ratio[isp,6]=max(Report$Z_adj[isp,])/max_obs
}


# chuckhi ringed haul-out model (jday + MDSDA)
Data$thin_mu_logit=Data$thin_mu_logit_jday_chukchi
Data$Sigma_logit_thin = Data$Sigma_logit_thin_jday_chukchi
Data$Cells_reduced_s = Cells_reduced-1
Params = list("beta"=matrix(0*rnorm(n_species*n_beta),n_species,n_beta), 
              "logtau_z"=rep(0,n_species), "logkappa_z"=rep(-.9,n_species),
              "phi_logit"=rep(0.0,n_species+3),"p_logit"=rep(0,n_species+3),"thin_logit_i"=Data$thin_mu_logit_jday_chukchi,"MisID_pars"=Data$MisID_mu,
              "logit_Pup_prop"=Data$Pup_prop_mu, "Etainput_s"=Etainput_s)
#run without misID to start with
Params$MisID_pars = rep(-10,length(Params$MisID_pars))
Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
Obj$fn( Obj$par )
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max = 10000,maxit=10000))         #
Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
Report = Obj$report() 

#now run again w/ misID
logit <- function(x) log(x/(1-x))
Params$beta = as.numeric(t(Report$Beta))
Params$phi_logit = logit(Report$phi/Data$phi_max)
Params$p_logit = logit(Report$power-1)
Params$MisID_pars = MisID_orig
Params$logit_Pup_prop = logit(Report$Pup_prop)
Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=FALSE)
Obj$fn( Obj$par )
Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
Upper = 50
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max = 10000,maxit=10000))         #
Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
Report = Obj$report() 
mod_chukchi = Report

if(Opt$convergence==0)Converge[7]=1
LogL[7]=Report$jnll
N=cbind(N,Report$total_abundance)
Extrap_num = cbind(Extrap_num,0)
Extrap_ratio = cbind(Extrap_ratio,0)
for(isp in 1:4){
  max_obs = max(Report$Z_adj[isp,Data$S_i])
  Extrap_num[isp,7] = sum(Report$Z_adj[isp,] > max_obs)
  Extrap_ratio[isp,7]=max(Report$Z_adj[isp,])/max_obs
}

Extrap_num=cbind(matrix(0,4,2),Extrap_num[,-2])
Extrap_ratio=cbind(matrix(0,4,2),Extrap_ratio[,-2])
for(isp in 1:4){
  max_obs = max(mod_base$Z_adj[isp,Data$S_i])
  Extrap_num[isp,1] = sum(mod_base$Z_adj[isp,] > max_obs)
  Extrap_ratio[isp,1]=max(mod_base$Z_adj[isp,])/max_obs
  max_obs = max(mod_no_misID$Z_adj[isp,Data$S_i])
  Extrap_num[isp,2] = sum(mod_no_misID$Z_adj[isp,] > max_obs)
  Extrap_ratio[isp,2]=max(mod_no_misID$Z_adj[isp,])/max_obs
}


Sens_table = data.frame("model"=c("base","no misID","-depth","-dist_land","-dist_break","-dist_edge",
                                  "avail_jday","avail_Bering"),
                        "LogL" = c(-mod_base$jnll,-mod_no_misID$jnll,-LogL[-2]),
                        "N_bd" = c(mod_base$total_abundance[3],mod_no_misID$total_abundance[3],N[3,-2]),
                        "N_rn" = c(mod_base$total_abundance[2],mod_no_misID$total_abundance[2],N[2,-2]),
                        "N_rd" = c(mod_base$total_abundance[4],mod_no_misID$total_abundance[4],N[4,-2]),
                        "N_sd" = c(mod_base$total_abundance[1],mod_no_misID$total_abundance[1],N[1,-2]),
                        "sum(Ik)"=colSums(Extrap_num),
                        "max(Lambda)"=apply(Extrap_ratio,2,'max')
                        )
library(xtable)
print(xtable(Sens_table,digits=c(1,0,0,0,0,0,0,0,2)),include.rownames=FALSE)