# format BOSS data for TMB, including haulout distributions for Russian Bering 

library('RPostgreSQL')
library('sf')
library('rgeos')
library('sp')
library('Matrix')
set.seed(12345)  #just so pseudo-zeroes will always be in the same place

remotes::install_github('jmlondon/berchukFits') #London et al. 2021 availability predictions
library(berchukFits)
data(bearded_fit)
fit_bearded=bearded_fit
data(ribbon_fit)
fit_ribbon=ribbon_fit
data(spotted_fit)
fit_spotted=spotted_fit

load("misID_priors_Russia.RData")  #created in format_okhotsk_TMB.R


load("c:/users/paul.conn/git/bass/BeringData2012_2013.rdat")
BOSS_data = Data$Grid
BOSS_grid = BOSS_data$y2012[[1]]

load("c:/users/paul.conn/git/bass/AlaskaBeringData2012_2013_14Dec2015.RDat")
BOSS_US = Data$Grid$y2012[[1]]



#remove US grid cells from combined Bering surface
Centroids = gCentroid(BOSS_grid,byid=TRUE)
I.Alaska=as.vector(gIntersects(Centroids,gUnionCascaded(BOSS_US),byid=TRUE))
BOSS_grid = BOSS_grid[I.Alaska==0,]
n_cells = length(BOSS_grid)

#add in "strata" - used in Chernook et al. 2018 Russian analysis
Centroids = Centroids[I.Alaska==0,]
LatLong = coordinates(spTransform(Centroids,CRS("+proj=longlat +datum=WGS84")))
BOSS_grid$Strata='M'
BOSS_grid$Strata[LatLong[,2]>62.4]="N"
BOSS_grid$Strata[LatLong[,1]<167.8 & LatLong[,1]>0]="S"
BOSS_grid$Strata=factor(BOSS_grid$Strata)

#plot(st_as_sf(BOSS_grid)[,"Strata"])


#limit to dates Russian surveys were flown and to grid cells in Russian waters
# 2012 BOSS_grid goes 4/5-5/22, surveys go 4/20 - 5/5 Russian time, which are 4/19-5/4 U.S. time
# 2013 BOSS grid goes 4/5-5/9, surveys go 4/11 - 4/30 Russian time, 4/10-4/29 US
t_steps_2012 = 16
t_steps_2013 = 20
BOSS_data2 = list(y2012 = vector("list",16),y2013 = vector("list",20))
for(idate in 15:30){
  BOSS_data2$y2012[[idate-14]]=BOSS_data$y2012[[idate]][I.Alaska==0,]
  BOSS_data2$y2012[[idate-14]]$Strata=BOSS_grid$Strata
}
for(idate in 6:25){
  BOSS_data2$y2013[[idate-5]]=BOSS_data$y2013[[idate]][I.Alaska==0,]
  BOSS_data2$y2013[[idate-5]]$Strata=BOSS_grid$Strata
} 


BOSS_effort = read.csv('Russia_BOSS_effort.csv',header=TRUE)
#remove a number of cells that effort wasn't randomly allocated
BOSS_effort = BOSS_effort[-which(BOSS_effort$Id %in% c(903:905,990,992,993) & BOSS_effort$Year==2012),]

laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
spdf <- SpatialPoints(coords = BOSS_effort[,5:6],
                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
BOSS_effort_sp = spTransform(spdf,CRS(laea_180_proj))
plot(BOSS_grid)
plot(BOSS_effort_sp,col='blue',add=T)
#remove effort that intersects w/ U.S. BOSS Grid (won't want to double count)
I.intersect = as.vector(gIntersects(BOSS_effort_sp,gUnionCascaded(BOSS_grid),byid=TRUE))

BOSS_effort_sp = SpatialPointsDataFrame(BOSS_effort_sp,BOSS_effort[,1:14])
BOSS_effort_sp = BOSS_effort_sp[I.intersect,]

Grid_wBS = BOSS_grid
save(Grid_wBS,file="Grid_wBS.RData")

#remove effort where number of seals counted > number detected by IR (glare issues)
BOSS_effort[which(BOSS_effort[,13]=="-"),13]="0"
BOSS_effort[,13]=as.numeric(as.character(BOSS_effort[,13]))  #for some reason this was a factor
N_seals = rowSums(BOSS_effort[,c(8:13)])
Which_greater = which(BOSS_effort[,"IR"]<N_seals)

BOSS_effort_sp = BOSS_effort_sp[-Which_greater,]  #removes 41 cells!

#rename
new_IDs = c(1:nrow(BOSS_effort_sp))
rownames(BOSS_effort_sp@data)=new_IDs


#check to see if any cells repeated
Unique.dates = unique(BOSS_effort_sp$Date)
for(i in 1:length(Unique.dates)){
  Which = which(BOSS_effort_sp$Date==Unique.dates[1])
  Unique.IDs = unique(BOSS_effort_sp@data[Which,"Id"])
  if(length(Unique.IDs)!=length(Which))cat(paste("Duplicated cell, date ",i,"\n"))
}
#none repeated


n_surveyed=nrow(BOSS_effort_sp)

#which cell surveyed
rownames(BOSS_grid@data)=c(1:length(BOSS_grid))  #rename IDs
Grid_sampled = rep(0,n_surveyed)
Intersects = gIntersects(BOSS_effort_sp,BOSS_grid,byid=TRUE)
for(i in 1:n_surveyed)Grid_sampled[i]=which(Intersects[,i]==1)

#Produce Day, Hour in solar hour for haulout predictions
Cell_IDs=as.numeric(rownames(BOSS_grid@data))
Longitude.sampled = coordinates(spTransform(BOSS_effort_sp,CRS("+proj=longlat +datum=WGS84")))[,1]
ymd = paste(BOSS_effort_sp$Year,BOSS_effort_sp$Date,sep='/') #convert date time to POSIXct
ymdt = paste(ymd,BOSS_effort_sp$Time)
BOSS_effort_sp$dt = as.POSIXct(ymdt,format="%Y/%d/%m %H:%M",tz="Asia/Magadan")
BOSS_effort_sp$SolarT = solaR::local2Solar(BOSS_effort_sp$dt,lon=Longitude.sampled)

#split into 2012, 2013
BE_2012 = BOSS_effort_sp[BOSS_effort_sp$Year==2012,]
BE_2013 = BOSS_effort_sp[BOSS_effort_sp$Year==2013,]
rownames(BE_2013@data)=c(1:nrow(BE_2013))

#Day and Hour, 2012
date.start=as.Date("2012-04-20") #note this is in Russian time
DayHour_2012=matrix(0,nrow(BE_2012),2)
DayHour_2012[,1] = as.numeric(as.Date(BE_2012$SolarT)-date.start+1)  
DayHour_2012[,2] = as.numeric(strftime(BE_2012$SolarT, format="%H"))

#Day and Hour, 2013
date.start=as.Date("2013-04-11")
DayHour_2013=matrix(0,nrow(BE_2013),2)
DayHour_2013[,1] = as.numeric(as.Date(BE_2013$SolarT)-date.start+1)  
DayHour_2013[,2] = as.numeric(strftime(BE_2013$SolarT, format="%H"))

#haulout predictions at date and time 
# 2012 BOSS_grid goes 4/5-5/22, surveys go 4/20 - 5/5 Russian time, which are 4/19-5/4 U.S. time
# 2013 BOSS grid goes 4/5-5/9, surveys go 4/11 - 4/30 Russian time, 4/10-4/29 US
con <- RPostgreSQL::dbConnect(PostgreSQL(), 
                              dbname = Sys.getenv("pep_db"), 
                              host = Sys.getenv("pep_ip"), 
                              user = Sys.getenv("pep_user"), 
                              rstudioapi::askForPassword(paste("Enter your DB password for user account: ", Sys.getenv("pep_user"), sep = "")))

grid.sf <- sf::st_read(con, 
                          query = "SELECT * FROM base.geo_analysis_grid", 
                          geometry_column = "geom")  
fast_ice_nic_sf = sf::st_read(con, 
                                 query = "SELECT * FROM surv_boss.geo_fastice_2013_nic", 
                                 geometry_column = "geom")
fast_ice_nws_sf = sf::st_read(con, 
                                 query = "SELECT * FROM surv_boss.geo_fastice_2013_nws", 
                                 geometry_column = "geom")
grid.wx.2012 <- RPostgreSQL::dbGetQuery(con, "SELECT * FROM base.tbl_analysis_grid_cov_wx 
                                   WHERE fdatetime_range_start BETWEEN '2012-04-19 00:00:00' AND '2012-05-05 00:00:00'")

grid.wx.2013 <- RPostgreSQL::dbGetQuery(con, "SELECT * FROM base.tbl_analysis_grid_cov_wx 
                                   WHERE fdatetime_range_start BETWEEN '2013-04-10 00:00:00' AND '2013-04-30 00:00:00'")
n = nrow(BOSS_effort_sp)
Covs = data.frame(matrix(0,n,6))
Hour = c(DayHour_2012[,2],DayHour_2013[,2])
colnames(Covs)=c("day","hour","precip","temp","pressure","wind")
#determine closest grid cell # to each location surveyed
Distances = st_distance(st_as_sf(BOSS_effort_sp),st_centroid(grid.sf))
Which_closest = rep(0,n)
for(i in 1:n)Which_closest[i]=which(Distances[i,]==min(Distances[i,]))
Cell_num = grid.sf[Which_closest,]$cell
for(i in 1:length(BE_2012)){
  Cur_dat = grid.wx.2012[which(grid.wx.2012$cell==Cell_num[i]),]
  DT_diff = abs(Cur_dat$fdatetime_range_start-BE_2012[i,]$dt)
  cur_row = which(DT_diff == min(DT_diff))[1]
  Covs[i,3:5]=Cur_dat[cur_row,c("rast_acpcp","rast_air2m","rast_prmsl")]
  Covs[i,"wind"] = sqrt(Cur_dat[cur_row,"rast_uwnd"]^2+Cur_dat[cur_row,"rast_vwnd"]^2)
}
for(i in 1:length(BE_2013)){
  Cur_dat = grid.wx.2012[which(grid.wx.2012$cell==Cell_num[i+length(BE_2012)]),]
  DT_diff = abs(Cur_dat$fdatetime_range_start-BE_2013[i,]$dt)
  cur_row = which(DT_diff == min(DT_diff))[1]
  Covs[i+length(BE_2012),3:5]=Cur_dat[cur_row,c("rast_acpcp","rast_air2m","rast_prmsl")]
  Covs[i+length(BE_2012),"wind"] = sqrt(Cur_dat[cur_row,"rast_uwnd"]^2+Cur_dat[cur_row,"rast_vwnd"]^2)
}
#transform covariates to scale used in GLMPMs
Covs[,"temp"]=(Covs[,"temp"]-270)/27
Covs[,"pressure"]=(Covs[,"pressure"]-100000)/10000
Covs[,"wind"]=Covs[,"wind"]/10
#now convert day, hour into format used in GLMPMs
Covs$day= (lubridate::yday(BOSS_effort_sp$dt)-120)/10
Covs$day2=Covs$day^2
Covs$day3=Covs$day^3
Covs$hour = factor(lubridate::hour(BOSS_effort_sp$SolarT),levels=c(0:23))
Hour = lubridate::hour(BOSS_effort_sp$SolarT)
Sin1 = sin(pi*Hour/12)
Cos1 = cos(pi*Hour/12)
Sin2 = sin(pi*Hour/6)
Cos2 = cos(pi*Hour/6)
Sin3 = sin(pi*Hour/4)
Cos3 = cos(pi*Hour/4)
AS.vec = c("ADULT.F","ADULT.M","SUBADULT","YOUNG OF YEAR")
L_list = vector("list",4)  
Covs$temp2 = Covs$temp
Covs$dry = rep(0,n)  #needed for model.matrix
Covs$sin1 = Sin1
Covs$sin2 = Sin2
Covs$sin3 = Sin3
Covs$cos1 = Cos1
Covs$cos2 = Cos2
Covs$cos3 = Cos3
Covs$northing = 1.0 #value 'typical' of tagged bearded seals in the Bering   #   -coordinates(BOSS_effort_sp)[,2]/3152522  #denominator is mean of spatial points used in analysis; see process_haulout_data4.R in haulout code/paper

# ------------------------------------------------------------------------------
#                   RIBBON SEALS
# ------------------------------------------------------------------------------

FE = fit_ribbon$fixed.effects
npar=nrow(fit_ribbon$covb)
for(iage in 1:4){
  #design matrix
  Covs$age_sex=factor(AS.vec[iage],levels=AS.vec)
  L_list[[iage]]=model.matrix(fit_ribbon$fixed.formula,data=Covs)
}
L = matrix(0,n*4,npar)
counter=1
for(i in 1:n){
  for(iage in 1:4){
    L[counter,]= L_list[[iage]][i,]
    counter=counter+1
  }
}

Ell = L%*%FE$estimate[-which(FE$estimate==0)]
Cell = L%*%fit_ribbon$covb%*%t(L)
Pred = plogis(Ell)
Deriv = diag(as.numeric(exp(Ell)/(1+exp(Ell))^2),n*4,n*4)
Pred_var = Deriv%*%crossprod(Cell,Deriv)

#stable stage distribution combo prediction
Pi = matrix(c(0.36,0.37,0.13,0.14),1,4)
M = kronecker(diag(n),Pi)
Pred_ho_ribbon = M %*% Pred
Var_ho_ribbon = M %*% tcrossprod(Pred_var,M)

n.12 = length(BE_2012)
n.13 = length(BE_2013)
Pred_ho_ribbon_2012 = matrix(Pred_ho_ribbon[1:n.12],ncol=1)
Pred_ho_ribbon_2013 = matrix(Pred_ho_ribbon[(n.12+1):n],ncol=1)
Var_ho_ribbon_2012 = Var_ho_ribbon[1:n.12,1:n.12]
Var_ho_ribbon_2013 = Var_ho_ribbon[(n.12+1):n,(n.12+1):n]

# ------------------------------------------------------------------------------
#                   Bearded SEALS
# ------------------------------------------------------------------------------

npar=nrow(fit_bearded$covb)
FE = fit_bearded$fixed.effects
#for(iage in 1:4){
  #design matrix
 # Covs$age.sex=factor(AS.vec[iage],levels=AS.vec)
 # L_list[[iage]]=model.matrix(fit_bearded$fixed.formula,data=Covs)
#}
#L = matrix(0,n*4,npar)
#counter=1
#for(i in 1:n){
#  for(iage in 1:4){
#    L[counter,]= L_list[[iage]][i,]
#    counter=counter+1
#  }
#}
L = model.matrix(fit_bearded$fixed.formula,data=Covs)

Ell = L%*%FE$estimate  #[-which(FE$estimate==0)]
Cell = L%*%fit_bearded$covb%*%t(L)
Pred = plogis(Ell)
Deriv = diag(as.numeric(exp(Ell)/(1+exp(Ell))^2),n,n)
Pred_var = Deriv%*%crossprod(Cell,Deriv)

#stable stage distribution combo prediction  #bearded seal model no longer uses age_sex
# Pi = matrix(c(0.27,0.23,0.38,0.12),1,4)
# M = kronecker(diag(n),Pi)
# Pred_ho_bearded = M %*% Pred
# Var_ho_bearded = M %*% tcrossprod(Pred_var,M)

Pred_ho_bearded_2012 = matrix(Pred[1:n.12],ncol=1)
Pred_ho_bearded_2013 = matrix(Pred[(n.12+1):n],ncol=1)
Var_ho_bearded_2012 = Pred_var[1:n.12,1:n.12]
Var_ho_bearded_2013 = Pred_var[(n.12+1):n,(n.12+1):n]

# ------------------------------------------------------------------------------
#                   Spotted SEALS
# ------------------------------------------------------------------------------

npar=nrow(fit_spotted$covb)
FE = fit_spotted$fixed.effects
for(iage in 1:4){
  #design matrix
  Covs$age_sex=factor(AS.vec[iage],levels=AS.vec)
  L_list[[iage]]=model.matrix(fit_spotted$fixed.formula,data=Covs)
}
L = matrix(0,n*4,npar)
counter=1
for(i in 1:n){
  for(iage in 1:4){
    L[counter,]= L_list[[iage]][i,]
    counter=counter+1
  }
}

Ell = L%*%FE$estimate[-which(FE$estimate==0)]
Cell = L%*%fit_spotted$covb%*%t(L)
Pred = plogis(Ell)
Deriv = diag(as.numeric(exp(Ell)/(1+exp(Ell))^2),n*4,n*4)
Pred_var = Deriv%*%crossprod(Cell,Deriv)

#stable stage distribution combo prediction
Pi = matrix(c(0.31,0.24,0.33,0.12),1,4)
M = kronecker(diag(n),Pi)
Pred_ho_spotted = M %*% Pred
Var_ho_spotted = M %*% tcrossprod(Pred_var,M)
Pred_ho_spotted_2012 = matrix(Pred_ho_spotted[1:n.12],ncol=1)
Pred_ho_spotted_2013 = matrix(Pred_ho_spotted[(n.12+1):n],ncol=1)
Var_ho_spotted_2012 = Var_ho_spotted[1:n.12,1:n.12]
Var_ho_spotted_2013 = Var_ho_spotted[(n.12+1):n,(n.12+1):n]



# ------------------------------------------------------------------------------
#                   Ringed SEALS
# ------------------------------------------------------------------------------

#NB: make 2 predictions: one using "Chukchi" seals, and one using "Bering" seals 
#The Chukchi model better reflects adult behavior than the Bering
#which is biased towards subadults with lower haul-out probabilities.  However, the Bering probably
#has more similar environmental conditions

# I reverted to a haul-out model with just jday*Sea and MDSDA*Sea; models with Markus or spring_temp
# gives really weird predictions (like low haulout with high temps)

#NB2: assuming subnivean lair use proportion from model estimated in CHESS paper

load("c:/users/paul.conn/git/haulout/ringed/ringed_gam_MDSDA.RData")  
load('c:/users/paul.conn/git/haulout/ringed/snow_melt_grid_Feb2021.Rda')


DT = c(BE_2012$SolarT,BE_2013$SolarT)
Covs_ringed = data.frame(sea = factor(rep("Chukchi",n),levels=c("Bering","Chukchi")),
                         jday= lubridate::yday(DT))
Grid_rows = rep(0,n)
for(i in 1:n)Grid_rows[i]=which(grid_sf$cell==Cell_num[i])
Covs_ringed$dmelt_MDSDA = Covs_ringed$jday - c(grid_sf$MDSDA_2012[Grid_rows[1:n.12]],grid_sf$MDSDA_2013[Grid_rows[(n.12+1):n]])
Covs_ringed$dmelt_Markus = Covs_ringed$jday - c(grid_sf$Markus_2012[Grid_rows[1:n.12]],grid_sf$Markus_2013[Grid_rows[(n.12+1):n]])

load("c:/users/paul.conn/git/haulout/ringed/Temps_cell_year.RDa")  #waiting to do this until temp data avail on larger grid
Covs_ringed$spring_temp = c(Cell_temp_yr[Grid_rows[1:n.12],"2012"],Cell_temp_yr[Grid_rows[(n.12+1):n],"2013"])

Pred_ho_ringed_chukchi = mgcv::predict.gam(m_MDSDA,newdata=Covs_ringed,type="response",se.fit=TRUE)
summary(Pred_ho_ringed_chukchi$fit)

Xp = predict(m_MDSDA,Covs_ringed,type="lpmatrix")
Pred = Xp%*%coef(m_MDSDA)
Varp = Xp %*% m_MDSDA$Vp %*% t(Xp)
Deriv = diag(as.numeric(exp(Pred)/(1+exp(Pred))^2))
Pred_var_chukchi = Deriv%*%crossprod(Varp,Deriv)

Covs_ringed$sea = factor(rep("Bering",n),levels=c("Bering","Chukchi"))
Pred_ho_ringed_bering = mgcv::predict.gam(m_MDSDA,newdata=Covs_ringed,type="response",se.fit=TRUE)
summary(Pred_ho_ringed_bering$fit)

Xp = predict(m_MDSDA,Covs_ringed,type="lpmatrix")
Pred = Xp%*%coef(m_MDSDA)
Varp = Xp %*% m_MDSDA$Vp %*% t(Xp)
Deriv = diag(as.numeric(exp(Pred)/(1+exp(Pred))^2))
Pred_var_bering = Deriv%*%crossprod(Varp,Deriv)

#now predict snow lair use 
set.seed(12345)
Melt_jittered = Covs_ringed$dmelt_MDSDA+rnorm(n,0,0.001) 
Covs_ringed$dmelt_MDSDA = Melt_jittered  #so we don't have repeated entries / VC matrix singularity
save(Covs_ringed,file="Covs_ringed_wBS.RData")  

#******  now need to run CHESS model with these covariates to get predictions and VC*****
#since the usual trick of doing X Beta X' 
#doesn't work because of constraints (penalty) on melt values in the TMB estimation procedure for availability as a
#function of snow melt date.  Instead, get predictions, VC straight from TMB CHESS model

load("c:/users/paul.conn/git/CHESS/melt_avail_wBS.RData") #produced in run_CHESS_ho_pspline_Mar2021.R
Avail_pred = melt_avail$Pred_melt
Avail_VC = melt_avail$VC_melt
Avail_chukchi = as.numeric(Avail_pred) * Pred_ho_ringed_chukchi$fit
Avail_bering = as.numeric(Avail_pred) * Pred_ho_ringed_bering$fit

#extension of Goodman's (1960) exact formula for variance of product; see Ver Hoef et al. 2014 appendix b
# note there is some ambiguity to do with the sign of the last term (known variance vs. unbiased estimation) - I've chosen the latter
# overall the effect will be small because the third term is much smaller than the other two
Var_melt_chukchi = as.numeric(Avail_pred*Avail_pred)*Pred_var_chukchi + as.numeric(Pred_ho_ringed_chukchi$fit*Pred_ho_ringed_chukchi$fit)*Avail_VC -  Avail_VC * Pred_var_chukchi
Var_melt_bering = as.numeric(Avail_pred*Avail_pred)*Pred_var_bering + as.numeric(Pred_ho_ringed_bering$fit*Pred_ho_ringed_bering$fit)*Avail_VC - Avail_VC * Pred_var_bering

#output means, var-cov matrices for HO predictions
HO_out = list(Mu_bd_12 = Pred_ho_bearded_2012,Mu_sd_12=Pred_ho_spotted_2012,Mu_rn_12=Pred_ho_ribbon_2012,
              Var_bd_12 = Var_ho_bearded_2012, Var_sd_12 = Var_ho_spotted_2012,Var_rn_12 = Var_ho_ribbon_2012,
              Mu_bd_13 = Pred_ho_bearded_2013,Mu_sd_13=Pred_ho_spotted_2013,Mu_rn_13=Pred_ho_ribbon_2013,
              Var_bd_13 = Var_ho_bearded_2013, Var_sd_13 = Var_ho_spotted_2013,Var_rn_13 = Var_ho_ribbon_2013,
              Mu_rd_12_chukchi = Avail_chukchi[1:n.12],
              Mu_rd_12_bering =  Avail_bering[1:n.12],
              Mu_rd_13_chukchi = Avail_chukchi[(n.12+1):n],
              Mu_rd_13_bering = Avail_bering[(n.12+1):n],
              Var_rd_12_chukchi = Var_melt_chukchi[1:n.12,1:n.12],
              Var_rd_12_bering = Var_melt_bering[1:n.12,1:n.12],
              Var_rd_13_chukchi = Var_melt_chukchi[(n.12+1):n,(n.12+1):n],
              Var_rd_13_bering = Var_melt_bering[(n.12+1):n,(n.12+1):n]
              )
save(HO_out,file="Avail_dists_BOSSrussia.RData")

# composite thinning expectation and variance
Sigma_thin = vector("list",4)
# put in disturbance
n_surveyed = n.12
mu_d = 23/24
var_d = mu_d*(1-mu_d)/24
Mu_dp = rep(mu_d,n_surveyed)
Var_dp = diag(var_d,n_surveyed)
Mu_spotted = as.numeric(HO_out$Mu_sd_12)*Mu_dp
X = Y = matrix(0,n_surveyed,n_surveyed)
for(i in 1:n_surveyed){
  for(j in 1:n_surveyed){
    X[i,j]= HO_out$Mu_sd_12[i]*HO_out$Mu_sd_12[j]
    Y[i,j]= Mu_dp[i]*Mu_dp[j]
  }
}
Sigma_thin[[1]] = Y*HO_out$Var_sd_12 +
  X*diag(Var_dp)-
  HO_out$Var_sd_12*diag(Var_dp)

# ribbon seals - no extra uncertainty other than availabilty
Mu_ribbon = as.numeric(HO_out$Mu_rn_12)
Sigma_thin[[2]] = HO_out$Var_rn_12

# bearded seals
mu_d = 48/51
var_d = mu_d*(1-mu_d)/51
Mu_dp = rep(mu_d,n_surveyed)
Var_dp = diag(var_d,n_surveyed)
Mu_bearded = as.numeric(HO_out$Mu_bd_12)*Mu_dp
X = Y = matrix(0,n_surveyed,n_surveyed)
for(i in 1:n_surveyed){
  for(j in 1:n_surveyed){
    X[i,j]= HO_out$Mu_bd_12[i]*HO_out$Mu_bd_12[j]
    Y[i,j]= Mu_dp[i]*Mu_dp[j]
  }
}
Sigma_thin[[3]] = Y*HO_out$Var_bd_12 +
  X*diag(Var_dp)-
  HO_out$Var_bd_12*diag(Var_dp)

#ringed
mu_d = 1-57/189
var_d = mu_d*(1-mu_d)/189
Mu_d = rep(mu_d,n_surveyed)
Var_d = diag(var_d,nrow=n_surveyed)
Mu_ringed_12_bering = Mu_d * HO_out$Mu_rd_12_bering
Mu_ringed_12_chukchi = Mu_d * HO_out$Mu_rd_12_chukchi
X12b = X12c = Y = matrix(0,n_surveyed,n_surveyed)
for(i in 1:n_surveyed){
  for(j in 1:n_surveyed){
    X12b[i,j]= HO_out$Mu_rd_12_bering[i]*HO_out$Mu_rd_12_bering[j]
    X12c[i,j]= HO_out$Mu_rd_12_chukchi[i]*HO_out$Mu_rd_12_chukchi[j]
    Y[i,j]= Mu_d[i]*Mu_d[j]
  }
}
Sigma_thin[[4]] = Y*HO_out$Var_rd_12_chukchi +
  X12c*diag(Var_d)-
  HO_out$Var_rd_12_chukchi*diag(Var_d)
Sigma_thin_12c = Sigma_thin

Sigma_thin[[4]] = Y*HO_out$Var_rd_12_bering +
  X12b*diag(Var_d)-
  HO_out$Var_rd_12_bering*diag(Var_d)
Sigma_thin_12b = Sigma_thin

Thin_12c = c(Mu_spotted, Mu_ribbon, Mu_bearded,Mu_ringed_12_chukchi)
Thin_12b = c(Mu_spotted, Mu_ribbon, Mu_bearded,Mu_ringed_12_bering)

Sigma_thin_12b = as.matrix(bdiag(Sigma_thin_12b))
Sigma_thin_12b = as(Sigma_thin_12b,"dgTMatrix")
Sigma_thin_12c = as.matrix(bdiag(Sigma_thin_12c))
Sigma_thin_12c = as(Sigma_thin_12c,"dgTMatrix")

#compute logit scale distribution using delta method
diff_logit <- function(x) -1/(x*(x-1))
Thin_logit_12b = log(Thin_12b/(1-Thin_12b))
Diff_logit = diag(diff_logit(Thin_12b))
Sigma_logit_thin_12b = Diff_logit %*% Sigma_thin_12b %*% t(Diff_logit)
diag(Sigma_logit_thin_12b)=diag(Sigma_logit_thin_12b)+0.000001  #adding a  small nugget to  prevent numerical problems

diff_logit <- function(x) -1/(x*(x-1))
Thin_logit_12c = log(Thin_12c/(1-Thin_12c))
Diff_logit = diag(diff_logit(Thin_12c))
Sigma_logit_thin_12c = Diff_logit %*% Sigma_thin_12c %*% t(Diff_logit)
diag(Sigma_logit_thin_12c)=diag(Sigma_logit_thin_12c)+0.000001  #adding a  small nugget to  prevent numerical problems


#repeat for 2013
Sigma_thin = vector("list",4)
# put in disturbance
n_surveyed = n.13
mu_d = 23/24
var_d = mu_d*(1-mu_d)/24
Mu_dp = rep(mu_d,n_surveyed)
Var_dp = diag(var_d,n_surveyed)
Mu_spotted = as.numeric(HO_out$Mu_sd_13)*Mu_dp
X = Y = matrix(0,n_surveyed,n_surveyed)
for(i in 1:n_surveyed){
  for(j in 1:n_surveyed){
    X[i,j]= HO_out$Mu_sd_13[i]*HO_out$Mu_sd_13[j]
    Y[i,j]= Mu_dp[i]*Mu_dp[j]
  }
}
Sigma_thin[[1]] = Y*HO_out$Var_sd_13 +
  X*diag(Var_dp)-
  HO_out$Var_sd_13*diag(Var_dp)

# ribbon seals - no extra uncertainty other than availabilty
Mu_ribbon = as.numeric(HO_out$Mu_rn_13)
Sigma_thin[[2]] = HO_out$Var_rn_13

# bearded seals
mu_d = 48/51
var_d = mu_d*(1-mu_d)/51
Mu_dp = rep(mu_d,n_surveyed)
Var_dp = diag(var_d,n_surveyed)
Mu_bearded = as.numeric(HO_out$Mu_bd_13)*Mu_dp
X = Y = matrix(0,n_surveyed,n_surveyed)
for(i in 1:n_surveyed){
  for(j in 1:n_surveyed){
    X[i,j]= HO_out$Mu_bd_13[i]*HO_out$Mu_bd_13[j]
    Y[i,j]= Mu_dp[i]*Mu_dp[j]
  }
}
Sigma_thin[[3]] = Y*HO_out$Var_bd_13 +
  X*diag(Var_dp)-
  HO_out$Var_bd_13*diag(Var_dp)

#ringed
mu_d = 1-57/189
var_d = mu_d*(1-mu_d)/189
Mu_d = rep(mu_d,n_surveyed)
Var_d = diag(var_d,nrow=n_surveyed)
Mu_ringed_13_bering = Mu_d * HO_out$Mu_rd_13_bering
Mu_ringed_13_chukchi = Mu_d * HO_out$Mu_rd_13_chukchi
X13b = X13c = Y = matrix(0,n_surveyed,n_surveyed)
for(i in 1:n_surveyed){
  for(j in 1:n_surveyed){
    X13b[i,j]= HO_out$Mu_rd_13_bering[i]*HO_out$Mu_rd_13_bering[j]
    X13c[i,j]= HO_out$Mu_rd_13_chukchi[i]*HO_out$Mu_rd_13_chukchi[j]
    Y[i,j]= Mu_d[i]*Mu_d[j]
  }
}
Sigma_thin[[4]] = Y*HO_out$Var_rd_13_chukchi +
  X13c*diag(Var_d)-
  HO_out$Var_rd_13_chukchi*diag(Var_d)
Sigma_thin_13c = Sigma_thin

Sigma_thin[[4]] = Y*HO_out$Var_rd_13_bering +
  X13b*diag(Var_d)-
  HO_out$Var_rd_13_bering*diag(Var_d)
Sigma_thin_13b = Sigma_thin

Thin_13c = c(Mu_spotted, Mu_ribbon, Mu_bearded,Mu_ringed_13_chukchi)
Thin_13b = c(Mu_spotted, Mu_ribbon, Mu_bearded,Mu_ringed_13_bering)

Sigma_thin_13b = as.matrix(bdiag(Sigma_thin_13b))
Sigma_thin_13b = as(Sigma_thin_13b,"dgTMatrix")
Sigma_thin_13c = as.matrix(bdiag(Sigma_thin_13c))
Sigma_thin_13c = as(Sigma_thin_13c,"dgTMatrix")

#compute logit scale distribution using delta method
diff_logit <- function(x) -1/(x*(x-1))
Thin_logit_13b = log(Thin_13b/(1-Thin_13b))
Diff_logit = diag(diff_logit(Thin_13b))
Sigma_logit_thin_13b = Diff_logit %*% Sigma_thin_13b %*% t(Diff_logit)
diag(Sigma_logit_thin_13b)=diag(Sigma_logit_thin_13b)+0.000001  #adding a  small nugget to  prevent numerical problems

diff_logit <- function(x) -1/(x*(x-1))
Thin_logit_13c = log(Thin_13c/(1-Thin_13c))
Diff_logit = diag(diff_logit(Thin_13c))
Sigma_logit_thin_13c = Diff_logit %*% Sigma_thin_13c %*% t(Diff_logit)
diag(Sigma_logit_thin_13c)=diag(Sigma_logit_thin_13c)+0.000001  #adding a  small nugget to  prevent numerical problems



n_species=4
n_transects_2012 = sum(BOSS_effort_sp$Year==2012)
n_transects_2013 = sum(BOSS_effort_sp$Year==2013) 
Mapping_2012 = Grid_sampled[BOSS_effort_sp$Year==2012]
Mapping_2013 = Grid_sampled[BOSS_effort_sp$Year==2013]

#for extra day effect on availability...
X_day = matrix(0,n_transects_2012,2)
X_day[,1] = DayHour_2012[,1]
X_day[,2] = (X_day[,1])^2 
X_day[,1]=(X_day[,1]-mean(X_day[,1]))/sqrt(var(X_day[,1]))
X_day[,2]=(X_day[,2]-mean(X_day[,2]))/sqrt(var(X_day[,2]))
X_day_2012 = kronecker(diag(n_species),X_day)

X_day = matrix(0,n_transects_2013,2)
X_day[,1] = DayHour_2013[,1]
X_day[,2] = (X_day[,1])^2 
X_day[,1]=(X_day[,1]-mean(X_day[,1]))/sqrt(var(X_day[,1]))
X_day[,2]=(X_day[,2]-mean(X_day[,2]))/sqrt(var(X_day[,2]))
X_day_2013 = kronecker(diag(n_species),X_day)

S_i_2012 = (DayHour_2012[,1]-1)*n_cells+Mapping_2012
S_i_2013 = (DayHour_2013[,1]-1)*n_cells+Mapping_2013

# compile habitat covariate data for each cell
n_hab_col=ncol(BOSS_grid@data)

Hab_cov = data.frame(matrix(0,n_cells*t_steps_2012,n_hab_col+2))
colnames(Hab_cov)=c(colnames(BOSS_grid@data),"ice2","depth2")
counter=1
for(it in 1:t_steps_2012){
  Tmp_data = BOSS_data2$y2012[[it]]@data
  Tmp_data$ice2 = Tmp_data$ice_conc^2
  Tmp_data$depth2 = Tmp_data$depth^2
  Hab_cov[counter:(counter+n_cells-1),]=Tmp_data
  counter=counter+n_cells
}
Hab_cov$Ecoregion=factor(Hab_cov$Ecoregion)
Hab_cov$depth=Hab_cov$depth/mean(Hab_cov$depth)  #note: this is the 2nd time this has been standardized (1st was wrt US + Russia grid)
Hab_cov$depth2=Hab_cov$depth^2
Hab_cov$sqrt_edge = sqrt(Hab_cov$dist_edge)
Northing = coordinates(BOSS_data2$y2012[[1]])[,2]
Northing = Northing/mean(Northing)
Easting = coordinates(BOSS_data2$y2012[[1]])[,1]
Easting = Easting/mean(Easting)
Hab_cov$easting = rep(Easting,t_steps_2012)
Hab_cov$easting2 = rep(Easting^2,t_steps_2012)
Hab_cov$northing = rep(Northing,t_steps_2012)
Hab_cov$northing2 = rep(Northing^2,t_steps_2012)

#attach fast ice
date.text.to.jday.leap <- function(date_text){  #convert april/may dates in text string to integer julian day - leap years (e.g. 2012)
  month<- as.numeric(substr(date_text,5,6))
  day<- as.numeric(substr(date_text,7,8))
  jday <- rep(0,length(date_text))
  jday[month==4]=day[month==4]+91
  jday[month==5]=day[month==5]+121
  jday
}
Hab_cov$fi_nic = 0

Jday_nic <- date.text.to.jday.leap(fast_ice_nic_sf$date_text)
#nic more difficult as we'll need to union different regions together and available dates are different for each region
Jdays_BW = sort(unique(Jday_nic[fast_ice_nic_sf$region=="Bering West"]))  #bering - west
Jdays_K = sort(unique(Jday_nic[fast_ice_nic_sf$region=="Kamchatka"]))  #kamchatka
W_remove = date.text.to.jday.leap(c("20120412","20120531","20120528"))
Jdays_BW = Jdays_BW[-which(Jdays_BW %in% W_remove)]

counter=0
cur_jday = 110  #surveys start April 20 in 2012
Grid_sf = st_as_sf(Grid_wBS)
cell_area = sf::st_area(Grid_sf[1,])

for(it in 1:t_steps_2012){
  #find fast ice date with julian day closest to day of survey
  Jday_diff_BW = abs(Jdays_BW - cur_jday)
  Jday_diff_K = abs(Jdays_K - cur_jday)
  
  which_BW = which(Jday_diff_BW==min(Jday_diff_BW))[1]
  which_K = which(Jday_diff_K==min(Jday_diff_K))[1]
  
  Rows = c(which(fast_ice_nic_sf$region=="Bering West" & Jday_nic==Jdays_BW[which_BW]),
           which(fast_ice_nic_sf$region=="Kamchatka" & Jday_nic==Jdays_K[which_K])
  )
  
  fi_union = sf::st_union(fast_ice_nic_sf[Rows,])
  
  fi_intersects = sf::st_intersects(fi_union,Grid_sf)
  n.inter = length(fi_intersects[[1]]==1)
  for(icell in 1:n.inter){
    cur_cell = fi_intersects[[1]][icell]
    inter_sf = sf::st_intersection(Grid_sf[cur_cell,],fi_union)
    Hab_cov$fi_nic[counter+fi_intersects[[1]][icell]]=sf::st_area(inter_sf)/cell_area
  }
  
  cur_jday = cur_jday+1
  counter = counter + n_cells
}


X_s_2012 = Hab_cov
#X_s_2012 = Hab_cov[,c(2,4,5,6,7,8,10,11,12,13,14)]  #choose in model script for now
X_s_2012[,"depth"]=sign(X_s_2012[,"depth"])*abs(X_s_2012[,"depth"])^(1/3)
X_s_2012[,"depth2"]=(X_s_2012[,"depth"])^2




Hab_cov = data.frame(matrix(0,n_cells*t_steps_2013,n_hab_col+2))
colnames(Hab_cov)=c(colnames(BOSS_grid@data),"ice2","depth2")
counter=1
for(it in 1:t_steps_2013){
  Tmp_data = BOSS_data2$y2013[[it]]@data
  Tmp_data$ice2 = Tmp_data$ice_conc^2
  Tmp_data$depth2 = Tmp_data$depth^2
  Hab_cov[counter:(counter+n_cells-1),]=Tmp_data
  counter=counter+n_cells
}
Hab_cov$Ecoregion=factor(Hab_cov$Ecoregion)
Hab_cov$depth=Hab_cov$depth/mean(Hab_cov$depth)  #note: this is the 2nd time this has been standardized (1st was wrt US + Russia grid)
Hab_cov$depth2=Hab_cov$depth^2
Hab_cov$sqrt_edge = sqrt(Hab_cov$dist_edge)
Northing = coordinates(BOSS_data2$y2013[[1]])[,2]
Northing = Northing/mean(Northing)
Easting = coordinates(BOSS_data2$y2013[[1]])[,1]
Easting = Easting/mean(Easting)
Hab_cov$easting = rep(Easting,t_steps_2013)
Hab_cov$easting2 = rep(Easting^2,t_steps_2013)
Hab_cov$northing = rep(Northing,t_steps_2013)
Hab_cov$northing2 = rep(Northing^2,t_steps_2013)

#attach fast ice
date.text.to.jday <- function(date_text){  #convert april/may dates in text string to integer julian day - non-leap years 
  month<- as.numeric(substr(date_text,5,6))
  day<- as.numeric(substr(date_text,7,8))
  jday <- rep(0,length(date_text))
  jday[month==4]=day[month==4]+90
  jday[month==5]=day[month==5]+120
  jday
}

Hab_cov$fi_nic = 0

Jday_nic <- date.text.to.jday(fast_ice_nic_sf$date_text)
#nic more difficult as we'll need to union different regions together and available dates are different for each region
Jdays_BW = sort(unique(Jday_nic[fast_ice_nic_sf$region=="Bering West"]))  #bering - west
Jdays_K = sort(unique(Jday_nic[fast_ice_nic_sf$region=="Kamchatka"]))  #kamchatka
nic_remove = date.text.to.jday(c("20130425"))
Jdays_BW = Jdays_BW[-which(Jdays_BW %in% nic_remove)]

counter=0
cur_jday = 100  #surveys start April 11 in 2013
Grid_sf = st_as_sf(Grid_wBS)
cell_area = sf::st_area(Grid_sf[1,])

for(it in 1:t_steps_2013){
  #find fast ice date with julian day closest to day of survey
  Jday_diff_BW = abs(Jdays_BW - cur_jday)
  Jday_diff_K = abs(Jdays_K - cur_jday)

  which_BW = which(Jday_diff_BW==min(Jday_diff_BW))[1]
  which_K = which(Jday_diff_K==min(Jday_diff_K))[1]

  Rows = c(which(fast_ice_nic_sf$region=="Bering West" & Jday_nic==Jdays_BW[which_BW]),
           which(fast_ice_nic_sf$region=="Kamchatka" & Jday_nic==Jdays_K[which_K])
  )
  
  fi_union = sf::st_union(fast_ice_nic_sf[Rows,])
  
  fi_intersects = sf::st_intersects(fi_union,Grid_sf)
  n.inter = length(fi_intersects[[1]]==1)
  for(icell in 1:n.inter){
    cur_cell = fi_intersects[[1]][icell]
    inter_sf = sf::st_intersection(Grid_sf[cur_cell,],fi_union)
    Hab_cov$fi_nic[counter+fi_intersects[[1]][icell]]=sf::st_area(inter_sf)/cell_area
  }
  
  cur_jday = cur_jday+1
  counter = counter + n_cells
}

X_s_2013=Hab_cov
#X_s_2013 = Hab_cov[,c(2,4,5,6,7,8,10,11,12,13,14)] #choose in model run R script for now
X_s_2013[,"depth"]=sign(X_s_2013[,"depth"])*abs(X_s_2013[,"depth"])^(1/3)
X_s_2013[,"depth2"]=(X_s_2013[,"depth"])^2


## format observed data
C_i_2012 = BE_2012@data[,c(10,11,9,8,12)]
C_i_2013 = BE_2013@data[,c(10,11,9,8,12)]

#for cases where seals observed in 0 sea ice, resample sea ice as uniform (0.05,0.30)
Ice_i_2012 = X_s_2012$ice_conc[S_i_2012]
Seals_i_2012 = rowSums(C_i_2012)
Which_ice0 = which(Ice_i_2012<0.01 & Seals_i_2012>0)
X_s_2012$ice_conc[S_i_2012[Which_ice0]]=runif(length(Which_ice0),0.05,0.30)  #8 instances in 2012
X_s_2012$ice2=X_s_2012$ice_conc^2

Ice_i_2013 = X_s_2013$ice_conc[S_i_2013]
Seals_i_2013 = rowSums(C_i_2013)
Which_ice0 = which(Ice_i_2013<0.01 & Seals_i_2013>0)  #no instances in 2013


#compute proportion of area surveyed using effective strip widths from Column 3, Table 1 of Chernook et al.
Width = rep(0,nrow(BE_2012))
Width[BE_2012$Date=="20/4"]=471.74
Width[BE_2012$Date=="21/4"]=453.22
Width[BE_2012$Date=="23/4"]=524.34
Width[BE_2012$Date=="27/4"]=417.60
Width[BE_2012$Date=="2/5"]=528.94
Width[BE_2012$Date=="3/5"]=331.98
Width[BE_2012$Date=="4/5"]=423.76
Width[BE_2012$Date=="5/5"]=673.84
Area_trans_2012 = Width*BE_2012@data$Length*1000/gArea(BOSS_grid[1,])


BE_2013$Width = 0
Width[BE_2013$Date=="11/4"]=447.72
Width[BE_2013$Date=="12/4"]=611.52
Width[BE_2013$Date=="19/4"]=667.48
Width[BE_2013$Date=="22/4"]=578.86
Width[BE_2013$Date=="23/4"]=630.36
Width[BE_2013$Date=="25/4"]=626.46
Width[BE_2013$Date=="30/4"]=437.16
Area_trans_2013 = Width*BE_2013@data$Length*1000/gArea(BOSS_grid[1,])

#add in pseudo_zeros
n_0 = 50

Which_ice0 = which(X_s_2012[,"ice_conc"]==0)
S_i_2012 = c(S_i_2012,sample(Which_ice0,n_0))
Which_ice0 = which(X_s_2013[,"ice_conc"]==0)
S_i_2013 = c(S_i_2013,sample(Which_ice0,n_0))

MisID_pos_rows = c(rep(1,4),rep(2,4),rep(3,4),rep(4,4))
MisID_pos_cols = c(2:5,1,3:5,1,2,4,5,1:3,5)

Pup_props = c(0.12,0.14,0.12,0.12)
Pup_prop_sigma = (diag(0.2*Pup_props))^2
Diff_logit = diag(diff_logit(Pup_props))
Sigma_pup_logit = Diff_logit %*% Pup_prop_sigma %*% t(Diff_logit)

Pups_2012 = BE_2012@data[,13] 
Pups_2012[Pups_2012=="-"]=0
Pups_2012=as.numeric(as.character(Pups_2012))
Photo_i = rowSums(BE_2012@data[,8:12])+Pups_2012
Tot_i = BE_2012[["IR"]]
Nophoto_i = Tot_i - Photo_i
Prop_photo_i = Photo_i/Tot_i
Which_undefined = which(is.na(Prop_photo_i) | Prop_photo_i==Inf)
Prop_photo_i[Which_undefined]=mean(Prop_photo_i[-Which_undefined])
Which.error = which(Nophoto_i<0)    #errors in reported survey data where # identified to species exceeds # detected by IR
Nophoto_i[Which.error]=0
Prop_photo_i[Which.error]=1.0
Nophoto_2012=Nophoto_i
Prop_photo_2012=Prop_photo_i

Pups_2013 = BE_2013@data[,13] 
Pups_2013[Pups_2013=="-"]=0
Pups_2013=as.numeric(as.character(Pups_2013))
Photo_i = rowSums(BE_2013@data[,8:12])+Pups_2013
Tot_i = BE_2013[["IR"]]
Nophoto_i = Tot_i - Photo_i
Prop_photo_i = Photo_i/Tot_i
Which_undefined = which(is.na(Prop_photo_i) | Prop_photo_i==Inf)
Prop_photo_i[Which_undefined]=mean(Prop_photo_i[-Which_undefined])
Which.error = which(Nophoto_i<0)    #errors in reported survey data where # identified to species exceeds # detected by IR
Nophoto_i[Which.error]=0
Prop_photo_i[Which.error]=1.0
Nophoto_2013=Nophoto_i
Prop_photo_2013=Prop_photo_i

Data_2012 = list("C_i"=C_i_2012, "Pups_i"=Pups_2012,"Nophoto_i" = Nophoto_2012, "Prop_photo_i"=Prop_photo_2012,"P_i"=Area_trans_2012,"A_s"=rep(1-Hab_cov[,"land_cover"],t_steps_2012),"S_i"=S_i_2012-1,"X_s"=X_s_2012,
                 "thin_mu_logit_bering"=Thin_logit_12b,"thin_mu_logit_chukchi"=Thin_logit_12c,
                 "Sigma_logit_thin_bering"=as(Sigma_logit_thin_12b,'dgTMatrix'),"Sigma_logit_thin_chukchi"=as(Sigma_logit_thin_12c,'dgTMatrix'),
                 "X_day"=X_day_2012,"MisID_mu" = misID_list$MisID_Mu, "MisID_Sigma"=misID_list$MisID_Sigma,"MisID_pos_rows"=MisID_pos_rows-1,"MisID_pos_cols"=MisID_pos_cols-1,"n_s" = n_cells, "n_sp" = n_species,"Pup_prop_mu" = qlogis(Pup_props),"Pup_prop_sigma" = diag(Sigma_pup_logit),"n_0"=n_0,"n_t"=t_steps_2012)
Data_2013 = list("C_i"=C_i_2013, "Pups_i"=Pups_2013,"Nophoto_i" = Nophoto_2013, "Prop_photo_i"=Prop_photo_2013,"P_i"=Area_trans_2013,"A_s"=rep(1-Hab_cov[,"land_cover"],t_steps_2013),"S_i"=S_i_2013-1,"X_s"=X_s_2013,
                 "thin_mu_logit_bering"=Thin_logit_13b,"thin_mu_logit_chukchi"=Thin_logit_13c,
                 "Sigma_logit_thin_bering"=as(Sigma_logit_thin_13b,'dgTMatrix'),"Sigma_logit_thin_chukchi"=as(Sigma_logit_thin_13c,'dgTMatrix'),
                 "X_day"=X_day_2013,"MisID_mu" = misID_list$MisID_Mu, "MisID_Sigma"=misID_list$MisID_Sigma,"MisID_pos_rows"=MisID_pos_rows-1,"MisID_pos_cols"=MisID_pos_cols-1,"n_s" = n_cells, "n_sp" = n_species,"Pup_prop_mu" = qlogis(Pup_props),"Pup_prop_sigma" = diag(Sigma_pup_logit),"n_0"=n_0,"n_t"=t_steps_2013)
Data_2012$h_mean = c(mean(HO_out$Mu_sd_12), mean(HO_out$Mu_rn_12), mean(HO_out$Mu_bd_12),0.65)
Data_2013$h_mean = c(mean(HO_out$Mu_sd_13), mean(HO_out$Mu_rn_13), mean(HO_out$Mu_bd_13),0.65)

save(Data_2012,file = "wBS_2012_Data_TMB.RData")
save(Data_2013,file = "wBS_2013_Data_TMB.RData")

