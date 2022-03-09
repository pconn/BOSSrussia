


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
                "MisID_pars"=Data$MisID_mu,
                "logit_Pup_prop"=Data$Pup_prop_mu,"beta0_ringed"=10,"beta1_ringed"=0)
Params$thin_logit_i = Data$thin_mu_logit_nolair_chukchi
Data$thin_mu_logit = Data$thin_mu_logit_nolair_chukchi
Data$Sigma_logit_thin = Data$Sigma_logit_thin_nolair_chukchi
Data$h_mean = 0 #not used currently

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
  Map[["beta0_ringed"]]=Map[["beta1_ringed"]]=factor(rep(NA,1))  #turn off ringed seal lair model
  
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
  
  
  #XY = coordinates(Grid)
  #rownames(XY) = c(1:nrow(XY))
  #plot_N_map_xy(N=Report$Z_s[1,],XY=XY) 
  
  
  #source('c:/users/paul.conn/git/ChukchiPower/ChukchiPower/R/util_funcs.R')
  #library(rgeos)
  library(RColorBrewer)
  library(ggplot2)
  plot_N_map<-function(N,Grid,No.plot=NA){
    #require(rgeos)
    require(sf)
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

  It = c(1,8,16)
  #Plot time-specific estimates; rows: bearded, ribbon, ringed, spotted, cols: Apr 19, April 26, May 4
  library(RColorBrewer)
  Bearded_plots=Ribbon_plots = Ringed_plots=Spotted_plots=vector("list",3)
  Titles = c("19 Apr","26 Apr","4 May")
  Max = matrix(0,4,3)
  for(it in 1:3){
    N_s= Report$Z_s[3,(It[it]-1)*Data$n_s+c(1:Data$n_s)]
    Bearded_plots[[it]]=plot_N_map(N_s,Grid=Grid)
    Bearded_plots[[it]]=Bearded_plots[[it]]+labs(fill='Bd N')+
      theme(axis.ticks = element_blank(), axis.text = element_blank(),axis.title = element_blank())+
      scale_fill_viridis_c(limits = c(0,7500), breaks = c(0, 2500, 5000, 7500),values=c(0,.1,.2,1))+ggtitle(Titles[it])
    Max[3,it] = max(N_s)
        
    N_s= Report$Z_s[2,(It[it]-1)*Data$n_s+c(1:Data$n_s)]
    Ribbon_plots[[it]]=plot_N_map(N_s,Grid=Grid)
    Ribbon_plots[[it]]=Ribbon_plots[[it]]+labs(fill='Rn N')+
      theme(axis.ticks = element_blank(), axis.text = element_blank(),axis.title = element_blank())+
      scale_fill_viridis_c(limits = c(0,1000), breaks = c(0, 300, 600, 900),values=c(0,.1,.2,1))
    Max[2,it] = max(N_s)
    
    N_s= Report$Z_s[4,(It[it]-1)*Data$n_s+c(1:Data$n_s)]
    Ringed_plots[[it]]=plot_N_map(N_s,Grid=Grid)
    Ringed_plots[[it]]=Ringed_plots[[it]]+labs(fill='Rd N')+
      theme(axis.ticks = element_blank(), axis.text = element_blank(),axis.title = element_blank())+
      scale_fill_viridis_c(limits = c(0,8000), breaks = c(0, 2500, 5000, 7500),values=c(0,.1,.2,1))
    Max[4,it] = max(N_s)
    
    N_s= Report$Z_s[1,(It[it]-1)*Data$n_s+c(1:Data$n_s)]
    Spotted_plots[[it]]=plot_N_map(N_s,Grid=Grid)
    Spotted_plots[[it]]=Spotted_plots[[it]]+labs(fill='Sd N')+
      theme(axis.ticks = element_blank(), axis.text = element_blank(),axis.title = element_blank())+
      scale_fill_viridis_c(limits = c(0,5000), breaks = c(0, 1500, 3000, 4500),values=c(0,.1,.2,1))
    Max[1,it] = max(N_s)
    
  }
  #pdf("Seal_N_st_plot.pdf",height=6,width=8)
  require(gridExtra)
  grid.arrange(Bearded_plots[[1]],Bearded_plots[[2]],Bearded_plots[[3]], 
               Ribbon_plots[[1]],Ribbon_plots[[2]],Ribbon_plots[[3]],
               Ringed_plots[[1]],Ringed_plots[[2]],Ringed_plots[[3]],
               Spotted_plots[[1]], Spotted_plots[[2]],Spotted_plots[[3]],ncol=3)
  
  pdf("wBS_N_plot_2012.pdf",height=10,width=8)
    grid.arrange(Bearded_plots[[1]],Bearded_plots[[2]],Bearded_plots[[3]], 
               Ribbon_plots[[1]],Ribbon_plots[[2]],Ribbon_plots[[3]],
               Ringed_plots[[1]],Ringed_plots[[2]],Ringed_plots[[3]],
               Spotted_plots[[1]], Spotted_plots[[2]],Spotted_plots[[3]],ncol=3)
  dev.off()
  
  #cowplot::plot_grid(Bearded_plot,Ringed_plot,ncol=1,plot.margin=margin(0,10,0,0))
  #dev.off
  
  
  ####### Covariate plots
  
  # scale of model
  Ice = c(0:100)/100
  Land = c(0:26)/10
  Break = c(0:23)/10
  Edge = c(0:27)/10
  Ice_df = data.frame(ice_conc=Ice,ice2=Ice^2,dist_mainland=0,dist_edge=0,dist_shelf=0)
  Break_df = data.frame(ice_conc=0,ice2=0,dist_mainland=0,dist_edge=0,dist_shelf=Break)
  Edge_df = data.frame(ice_conc=0,ice2=0,dist_mainland=0,dist_edge=Edge,dist_shelf=0)
  Land_df = data.frame(ice_conc=0,ice2=0,dist_mainland=Land,dist_edge=0,dist_shelf=0)
  DF_list = vector("list",4)
  DF_list[[1]] = Ice_df
  DF_list[[2]] = Break_df
  DF_list[[3]] = Edge_df
  DF_list[[4]] = Land_df
  Species = c("Spotted","Ribbon","Bearded","Ringed")
  Covariate = c("ice_conc","dist_shelf","dist_edge","dist_land")
  Strata = c("Middle","North","South")
  
  Effect_list = vector("list",4*4*3)
  counter=0
  for(isp in 1:4){
    for(icov in 1:4){
      for(istrata in 1:3){
        if(isp<3 & icov>2)crap=1 #don't do anything for these; not modeled
        else{
          counter=counter+1
          X_cur = DF_list[[icov]]
          X_cur$Strata = factor(istrata,levels=c("1","2","3"))
          if(isp==3 | isp==4)DM = model.matrix(~Strata*(ice_conc+ice2+dist_mainland+dist_shelf+dist_edge),X_cur)
          else DM = model.matrix(~Strata+ice_conc+ice2+dist_shelf,X_cur)
          DM = DM[,-1] #get rid of intercept
          if(isp==1) Pred = DM %*% Report$beta_sd
          if(isp==2) Pred = DM %*% Report$beta_rn
          if(isp==3) Pred = DM %*% Report$beta_bd
          if(isp==4) Pred = DM %*% Report$beta_rd
          #Note, Strata 1 = Middle, Strata 2 = North, Strata 3 = South
          
          Effect_list[[counter]]=data.frame("Effect"=Pred,"Species"=Species[isp],Strata=Strata[istrata],Series=counter,
                                           "Covariate"=Covariate[icov])
          if(icov==1)Effect_list[[counter]]$Value = Ice
          if(icov==2)Effect_list[[counter]]$Value = Break
          if(icov==3)Effect_list[[counter]]$Value = Edge
          if(icov==4)Effect_list[[counter]]$Value = Land
        }
      }
    }
  }
  Plot_df = Effect_list[[1]]
  for(i in 2:length(Effect_list))Plot_df = rbind(Plot_df,Effect_list[[i]])

  breaks_fun <- function(x) {
    if (max(x) > 2) {
      seq(0,2,1)
    } else {
      seq(0, 1, 0.5)
    }
  }
  
  cov_plot = ggplot(Plot_df)+geom_line(aes(x=Value,y=Effect,color=Strata,group=Series))+facet_grid(Species~Covariate,scales="free")+
    xlab("Covariate value") + scale_x_continuous(breaks = breaks_fun) 
  cov_plot
  
  pdf("wBS2012_cov_plot.pdf")
    cov_plot
  dev.off()
  
  
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
  
  
  
  
  

  #plot counts
  isp=3
  Count = rep(NA,nrow(Grid))
  Tmp = Data$S_i[1:nrow(Data$C_i)]+1  #recall Data$S_i defined for C++
  Count[Tmp %% Data$n_s]=Data$C_i[,isp]/(Data$P_i/Data$Ice_s[Data$S_i[1:252]])
  plot_N_map(N=Count,Grid=Grid) 
  
  #sum(Data$C_i[,isp]/Data$P_i)  #so even just sampled cells would appear to have more ribbon seals in 2012 than the whole grid in 2013.
  
  #plot adjusted counts as function of P_i (look for pref sampling)
  Tot_seals = (rowSums(Data$C_i)+Data$Nophoto_i)/Data$P_i
  plot(Data$P_i,Tot_seals)
  summary(lm(Tot_seals~Data$P_i))  #looks ok
  #take out the 2 obvious outliers
  #summary(lm(Tot_seals[-c(216,217)]~Data$P_i[-c(216,217)]))
  
  #plot Strata
  plot_N_map(N=as.numeric(as.character(X_s$Strata[1:Data$n_s])),Grid=Grid) 
  

  #plot counts as function of ice
  S_i_tmp = (Data$S_i+1)[1:252] #take out pseudozeroes
  plot(Data$X_s[,"ice_conc"][S_i_tmp],Data$C_i[,1],main="Spotted")
  plot(Data$X_s[,"ice_conc"][S_i_tmp],Data$C_i[,2],main="Ribbon")
  plot(Data$X_s[,"ice_conc"][S_i_tmp],Data$C_i[,3],main="Bearded")
  plot(Data$X_s[,"ice_conc"][S_i_tmp],Data$C_i[,4],main="Ringed")
  
  #plot ice
  it=10
  plot_N_map(N=Data$X_s[(it-1)*Data$n_s+c(1:Data$n_s),"ice_conc"],Grid=Grid) 

  plot_N_map(N=Data$X_s[(it-1)*Data$n_s+c(1:Data$n_s),"dist_shelf"],Grid=Grid) 
  
  it=16
  plot_N_map(N=Data$X_s[(it-1)*Data$n_s+c(1:Data$n_s),"fi_nic"],Grid=Grid) 
  
  Ice = c(0:100)/100
  plot(Ice*Report$beta_rn[3]+Ice^2*Report$beta_rn[5],ylim=c(-4,4))
  plot(Ice*Report$beta_sd[3]+Ice^2*Report$beta_sd[5],ylim=c(-4,4))
  
  
  lines(Ice*Report$Beta[3,4]+Ice^2*Report$Beta[3,7])
  
  plot(Ice*Report$Beta[1,4]+Ice^2*Report$Beta[1,7]+Ice^3*Report$Beta[1,16],ylim=c(-4,4))
  lines(Ice*Report$Beta[3,4]+Ice^2*Report$Beta[3,7]+Ice^2*Report$Beta[3,16])
  
 
  plot(Ice*Report$Beta[2,4]+Ice^2*Report$Beta[2,7])
  lines(Ice*Report$Beta[4,4]+Ice^2*Report$Beta[4,7])
  
  
  #tabulate counts of each species by strata
  Strata_ID = as.numeric(as.character(X_s[Data$S_i,"Strata"]))[1:252]
  Counts_sum = matrix(0,3,5)
  for(i in 1:3){
    Cs = Data$C_i[Strata_ID==i,]
    Counts_sum[i,]=colSums(Cs)
  }
  Counts_sum
  
  
  SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=TRUE )

  
  #randomized quantile residuals
  library(tweedie)
  set.seed(12345)
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
  
  pdf('wBS_GOF_2012.pdf')
  ggplot(Resids.df)+geom_histogram(aes(x=Residual),breaks=c(0:20)/20)+facet_wrap(~Labels)+ylab("Count")
  dev.off()
  

  # Extrapolation metrics
  
  #"base model' extrap values
  Extrap_num_base = Extrap_ratio_base =rep(0,4)
  for(isp in 1:4){
    Max_obs = max(Report$Z_s[isp,Data$S_i])
    Extrap_num_base[isp] = sum(Report$Z_s[isp,] > Max_obs)
    Extrap_ratio_base[isp]=max(Report$Z_s[isp,])/Max_obs
  }
  Extrap_num_base
  Extrap_ratio_base

