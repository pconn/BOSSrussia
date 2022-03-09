#produce_CIs

log_CI <- function(N,SE_N){
  CV=SE_N/N
  SE = sqrt(log(1+CV^2))
  C = exp(1.96*SE)
  CI = c(N/C,N*C)
  CI
}

#2012 wBS

load('boot_out_wBS2012.RData')
cat('bearded 2012 wBS \n')
Base$N[3]
log_CI(Base$N[3],SE.bd)

cat('ribbon 2012 wBS \n')
Base$N[2]
log_CI(Base$N[2],SE.rn)

cat('ringed 2012 wBS \n')
Base$N[4]
log_CI(Base$N[4],SE.rd)

cat('spotted 2012 wBS \n')
Base$N[1]
log_CI(Base$N[1],SE.sd)


#2013 wBS

load('boot_out_wBS2013.RData')
cat('bearded 2013 wBS \n')
Base$N[3]
log_CI(Base$N[3],SE.bd)

cat('ribbon 2013 wBS \n')
Base$N[2]
log_CI(Base$N[2],SE.rn)

cat('ringed 2013 wBS \n')
Base$N[4]
log_CI(Base$N[4],SE.rd)

cat('spotted 2013 wBS \n')
Base$N[1]
log_CI(Base$N[1],SE.sd)