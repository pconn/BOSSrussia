


library( RandomFields )
library( TMB )
library( INLA )
#library(TMBdebug)


# Compile
setwd("C:/Users/paul.conn/git/BOSSrussia/")
TmbFile = "C:/Users/paul.conn/git/BOSSrussia/BOSSrussia/src/OkhotskTMB"
#dyn.unload( dynlib(TmbFile) )  #unload file if previously loaded to get proper write permission
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="") 


#require(STabundance)
source('./BOSSrussia/R/util_funcs.R')
#source('./OkhotskSeal/R/sim_funcs.R')
#source('./OkhotskSeal/R/sim_data_generic.R')






# Settings
grid_dim = c("x"=25, "y"=25)
n_samp = 100
n_cells = grid_dim[1]*grid_dim[2]
prop_sampled=0.2
Prop_sampled=rep(prop_sampled,n_samp)
SpatialScale = sqrt(prod(grid_dim))/5  # Range ~ 2*Scale
SD_eta = SD_x = SD_delta = 1
beta0 = 3
Use_REML = FALSE   # 
Spatial_sim_model = "GP_gaussian"
Spatial_model = "SPDE_GMRF"
Alpha = 1  # Smoothness for GMRF, 1 or 2 (1 is faster)
RandomSeed = 123456
n_sim = 100
n_species = 4
n_obs_types = 5 #doesn't include pups
Prop_photo = runif(n_samp,-.2,1.2)
Prop_photo = runif(n_samp,0.8,0.95)
Prop_photo[Prop_photo<0]=0
Prop_photo[Prop_photo>1]=1 

Thin = rep(0.6,n_species*n_samp)  #keeping this as a vector for use w MVNORM_t in TMB
thin_logit = log(Thin/(1-Thin))
Sigma_logit_thin = vector("list",n_species)
for(isp in 1:n_species){
  Sigma_logit_thin[[isp]]=0.05*diag(n_samp)  #recall this is variance, not SD 
}
#convert into block diagonal 
Sigma_logit_thin = as.matrix(bdiag(Sigma_logit_thin))

#misID  - hardwired for 4 species
load('./misID_priors_Russia.RData')
MisID_pos_rows = c(rep(1,4),rep(2,4),rep(3,4),rep(4,4))
MisID_pos_cols = c(2:5,1,3:5,1,2,4,5,1:3,5)

#proportion pups by species
Pup_prop = rep(0.1,n_species)

# Results
Results = array(NA,dim=c(n_sim,n_species,6))
# columns are "sim","N_true","N_est","N_SE","Convergence","LogLik"

# Loop through replicates
counter=1
i=1
#for(i in 1:n_sim){
    cat(paste("Simulation ",i,"\n"))
    set.seed( RandomSeed + i )

    betax = betax_prob = runif(n_species,-.5,.5)       # Impact of covariate (x) on density

    # Spatial model
    loc_s = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y'])
    model_x <- RMgauss(var=SD_x^2, scale=SpatialScale)
    model_eta <- RMgauss(var=SD_eta^2, scale=SpatialScale)
      
    # Realization from GRF
    # x_s -- spatial variation in covariate
    Eta_s = matrix(0,n_species,n_cells)
    if( Spatial_sim_model=="GP_gaussian"){
       x_s = RFsimulate(model=model_x, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
       for(isp in 1:n_species){
         Eta_s[isp,] = RFsimulate(model=model_eta, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
       }
    }
    #if( Spatial_sim_model=="ICAR"){
      #x_s = rrw( Q )[,1]
      #eta_s = rrw( Q )[,1]
    #}
      
    # Total abundance
    Ztrue_s = Eta_s
    for(isp in 1:n_species)Ztrue_s[isp,] = exp( beta0 + betax[isp]*x_s + Eta_s[isp,] )
    #Ztrue_s = rpois(n_cells,Ztrue_s)     
    #plot_N_map_xy(N=Ztrue_s[1,],XY=loc_s,leg.title="Abundance")
      
    # Sampling intensity
    R_s = rep(1/n_cells,n_cells)
      
    # Process for locating samples
    s_i = sample(1:prod(grid_dim), size=n_samp, replace=FALSE, prob=R_s)
    y_s = ifelse(1:prod(grid_dim) %in% s_i, 1, 0)
      
    # Counting process
    Ctrue_i = matrix(0,n_species,n_samp)
    Pups_i = Nophoto_i=rep(0,n_samp)
    Sum_lambda = rep(0,n_samp)
    for(isp in 1:n_species){
      Ctrue_i[isp,] = rpois( n=n_samp, lambda=Prop_photo*Ztrue_s[isp,s_i]*prop_sampled*(1-Pup_prop[isp])*Thin[(isp-1)*n_samp+c(1:n_samp)])
      Pups_i = Pups_i + rpois(n=n_samp,lambda=Prop_photo*Ztrue_s[isp,s_i]*prop_sampled*Pup_prop[isp]*Thin[(isp-1)*n_samp+c(1:n_samp)])
      Sum_lambda = Sum_lambda + Ztrue_s[isp,s_i]*Thin[(isp-1)*n_samp+c(1:n_samp)]*prop_sampled
    }
    
    Nophoto_i = rpois(n=n_samp,(1-Prop_photo)*Sum_lambda)
    
    Tot_obs = colSums(Ctrue_i)+Pups_i+Nophoto_i
    
    
    #now confusion matrix on C_i
    C_i = matrix(0,n_samp,n_obs_types)
    C_i = round(t(Ctrue_i)%*%misID_list$Confuse)  #not strictly necessary to round since using Tweedie
    #C_i =t(Ctrue_i)%*%misID_list$Confuse
    
    # Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
    mesh = inla.mesh.create( loc_s )
      
       
    # Options
    #Options_vec = c("SE"=0)  #bias correction for beta, cell-level intensity?
    
    # Data
    spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
    #Data = list( "Options_vec"=Options_vec, "c_i"=c_i, "P_i"=Prop_sampled,"A_s"=rep(1,n_cells),"s_i"=s_i-1, "X_sj"=cbind(1,x_s), "y_s"=y_s, "X_sk"=cbind(1,x_s),"X_sb"=matrix(1,n_cells,1), "spde"=spde)
    Data = list( "C_i"=C_i, "Pups_i"=Pups_i,"Nophoto_i"=Nophoto_i,"Prop_photo_i" = Prop_photo,"P_i"=Prop_sampled,"A_s"=rep(1,n_cells),"S_i"=s_i-1,"Y_s"=y_s,"X_s"=matrix(1,n_cells,2), 
                 "spde"=spde,"thin_mu_logit"=thin_logit,"Sigma_logit_thin"=as(Sigma_logit_thin,'dgTMatrix'),
                 "MisID_mu" = misID_list$MisID_Mu, "MisID_Sigma"=misID_list$MisID_Sigma, "MisID_pos_rows"=MisID_pos_rows-1,
                 "MisID_pos_cols" = MisID_pos_cols-1,"n_s"=n_cells,"n_sp"=n_species,
                 "Pup_prop_mu" = rep(log(.1/.9),n_species),"Pup_prop_sigma" = rep(0.1,n_species))
    Data$X_s[,2]=x_s  
    
    # Parameters
    Etainput_s = matrix(0,n_species,mesh$n)  #spatial random effects
    Etainput_s[,1:625]=Eta_s
    
    Params = list("beta"=matrix(0,n_species,ncol(Data$X_s)), 
                  "logtau_z"=rep(0,n_species), "logkappa_z"=rep(-.9,n_species),
                  "phi_log"=rep(0.01,n_species+3),"p_logit"=rep(0,n_species+3),"thin_logit_i"=thin_logit,"MisID_pars"=misID_list$MisID_Mu,
                  "logit_Pup_prop"=rep(-2,n_species))
    Params$beta[,1]=beta0 + runif(n_species,-0.1,0.1) #set mean expected abundance close to truth for faster optimization
    if(ncol(Data$X_s)>1)Params$beta[,2:ncol(Data$X_s)]=betax + runif(n_species,-0.1,0.1)
    Params$Etainput_s = Etainput_s

    # Random
    Random = c( "Etainput_s" )
    #if(Use_REML==TRUE) Random = c(Random,"Beta")
    #Random = list()
    
    # Fix parameters
    Map = list()
    #Map[["phi_log"]]=factor(rep(1,n_species+3))
    Map[["phi_log"]]=factor(rep(NA,n_species+3))
    Map[["p_logit"]]=factor(rep(2,n_species+3))
    Map[["logkappa_z"]] = factor(rep(3,n_species))
    Map[["logtau_z"]] = factor(rep(4,n_species))
    
      
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
    Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, maxit=1000))         #
    Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
    Report = Obj$report()

    # Potentially fix random fields with zero sample or population variance
    #if( any(Report$MargSD_z<0.001) ){
      #Which = which(Report$MargSD_z<0.001)
     # Map[["logtau_z"]] = factor( ifelse(1:2==Which,NA,1:2) )
    #  if(length(Which)==2){
    #    Map[["logkappa_z"]] = factor( c(NA,NA) )
     # }
     # Map[["beta_b"]] = factor(NA)
    #  Params[["beta_b"]] = 0
    #  Map[["etainput_s"]] = factor( rep(NA,length(Params[["etainput_s"]])) )
    #  Params[["etainput_s"]][] = 0
    #  # Re-run
    #  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE)
    #  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, maxit=1000))         #
    #  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
    #}
    
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

    # Record results
    for(isp in 1:4)Results[counter,isp,1:5]=c(i,sum(Ztrue_s[isp,]), SD$unbiased$value[which(names(SD$value)=="total_abundance")][isp],SD$sd[which(names(SD$value)=="total_abundance")][isp],Converge)
    Results[counter,,6]=-Opt$objective

    #dyn.unload( dynlib(TmbFile) )

    # output
    save( Results, file="Results.RData")

    counter=counter+1
      # Plot stuff
      #par( mfrow=c(1,3) )
      #plot( x=Ztrue_s[s_i], y=c_i)
      #plot( x=Ztrue_s, y=R_s)
      #plot( x=Report$Z_s, y=Report$R_s)
}
  
# 
# ####################
# # Read results
# ####################
# TmbFile = "C:/Users/paul.conn/git/pref_sampling/" #TMB_version/inst/executables"
# setwd(TmbFile)
# load( "Results.RData")
# WhichDone = which(Results[,"Convergence"]==0)
# Results=Results[WhichDone,]
# Results[,"RelBias_N"]=(Results[,"N_est"]-Results[,"N_true"])/Results[,"N_true"]
# True.b = Results[,"b"]
# True.b[which(True.b==1)]=0
# True.b[which(True.b==2)]=1
# True.b[which(True.b==3)]=3
# Results[,"b"]=True.b
# Which.est = which(Results[,"EM"]>1)
# Results[,"Bias_b"]=NA
# Results[Which.est,"Bias_b"]=(Results[Which.est,"b_est"]-Results[Which.est,"b"])
# 
# #
# pdf("B_est.pdf")
#   Which_plot = which(Results[,"EM"]>1)
#   plot(Results[Which_plot,"b_est"],xlab="Simulation",ylab=expression(hat(b)))
# dev.off()
# 
# 
# #For figures, limit to results for which |b estimate| < 10 
# Results[,"OnBoundary"]=rep(0,nrow(Results))
# Results[which(abs(Results[,"b_est"])>10),"OnBoundary"]=1
# 
# Results_plot=Results[Results[,"OnBoundary"]==0,]
# Results_plot[,"Est.model"]=rep("Independent",nrow(Results_plot))
# Results_plot[Results_plot[,"EM"]==2,"Est.model"]="Joint"
# 
# library(doBy)
# sd <- function(x)sqrt(var(x))
# Bias_table = summaryBy(Bias_b+RelBias_N~b+EM,data=Results_plot,FUN=c(median,mean,sd))
# 
# #plot bias_N by b_est for EstMod=2
# Cur_data=Results_plot[Results_plot[,"EM"]==2,]
# pdf("Bias_bN.pdf")
# plot(Cur_data[,"Bias_b"],Cur_data[,"RelBias_N"],xlab="Bias (b)",ylab="Relative bias (N)")
# dev.off()
# 
# 
# 
# #plot bias as function of b, estimation method
# #first, rearrange relative bias in form for ggplot
# Bias.df=Results_plot
# Bias.df[,"B"]=Bias.df[,"b"]
# Bias.df[,"Bias"]=Bias.df[,"RelBias_N"]
# 
# 
# library(ggplot2)
# #plot proportion relative bias
# bias.plot = ggplot(Bias.df,aes(factor(Est.model),Bias))+geom_boxplot()+facet_grid(~b) #,scales="free")
# bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
# bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "Estimation model", y="Proportion relative bias")
# bias.plot
# pdf("bias.pdf")
# bias.plot
# dev.off()
# 
# #bias of b parameter
# BiasB.df = Bias.df[which(!is.na(Bias.df[,"Bias_b"])),]
# bias.plot = ggplot(BiasB.df,aes(factor(b),Bias_b))+geom_boxplot() #+facet_grid(~b) #,scales="free")
# bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
# bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "True b parameter", y="Absolute bias")
# bias.plot
# pdf("biasB.pdf")
# bias.plot
# dev.off()
# 
# #produce plot of bias of species-habitat relationship parameters
# bias.plot = ggplot(Bias.df,aes(factor(Est.model),BiasBeta11))+geom_boxplot()+facet_grid(~b) #,scales="free")
# bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
# bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "Estimation model", y="Absolute bias")
# bias.plot
# pdf("Bias_beta0.pdf")
# bias.plot
# dev.off()
# 
# #produce plot of bias of species-habitat relationship parameters
# bias.plot = ggplot(Bias.df,aes(factor(Est.model),BiasBeta12))+geom_boxplot()+facet_grid(~b) #,scales="free")
# bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
# bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "Estimation model", y="Absolute bias")
# bias.plot
# pdf("Bias_beta1.pdf")
# bias.plot
# dev.off()
# 
# 
# Results[,"Converge"]=Results[,"OnBoundary"]+Results[,"Convergence"]
# Results[which(Results[,"Converge"]==2),"Converge"]=1
# Converge_table = summaryBy(Converge~b+EM,data=Results,FUN=sum)
# 
# 
