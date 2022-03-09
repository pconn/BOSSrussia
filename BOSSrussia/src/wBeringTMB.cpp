// Space time 
#include <TMB.hpp>
//#include <atomic_math.hpp>  //for D_lgamma


template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}

template<class Type>
bool isNA(Type x) {
	return R_IsNA(asDouble(x));
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace Eigen;
  using namespace density;
  
  // options vec
  //DATA_FACTOR( Options_vec );
  // Slot 0: compute SE? 
  
  // Data
  DATA_MATRIX(C_i);       	// Matrix of responses (counts) of each species at each sampled location 
  DATA_VECTOR(Pups_i);      // vector holding observations of pups in each sampled cell
  DATA_VECTOR(Nophoto_i);   // vector holding observations of number of animals without a photograph
  DATA_VECTOR(Prop_photo_i);  // proportion of animals in cell i that were photographed 
  DATA_VECTOR( P_i );      // Proportion of survey unit that is surveyed for each observation i
  DATA_VECTOR( A_s );      // Relative area of each survey unit s (can set = 1.0 if all the same size)
  DATA_IVECTOR(S_i); // Site-time index for each sample
  DATA_MATRIX( X_sd );  //design matrix for fixed effects
  DATA_MATRIX( X_rn );  //design matrix for fixed effects
  DATA_MATRIX( X_bd );  //design matrix for fixed effects
  DATA_MATRIX( X_rd );  //design matrix for fixed effects
  DATA_IVECTOR(Which_zero); // set abundance to zero for cells w/ less than 1% ice
  DATA_VECTOR(Ice_s);  //proportion of grid cell containing ice

  DATA_VECTOR(thin_mu_logit);
  DATA_SPARSE_MATRIX(Sigma_logit_thin);
  DATA_MATRIX(X_day);  // design matrix for extra day and day^2 effects
  DATA_VECTOR(MisID_mu);
  DATA_MATRIX(MisID_Sigma);  //variance covariance matrix for MisID_pars on mlogit scale
  DATA_IVECTOR(MisID_pos_rows);   //row indices for positive entries of MisID matrix
  DATA_IVECTOR(MisID_pos_cols);   //ibid, columns
  DATA_INTEGER(n_s); //number of cells
  DATA_INTEGER(n_sp); //number of species
  DATA_INTEGER(n_0);  //number of pseudo-zeroes (included in S_i)
  DATA_INTEGER(n_t);  //number of time steps
  DATA_VECTOR(Pup_prop_mu);  //prior mean (logit scale) for pup proportion
  DATA_VECTOR(Pup_prop_sigma);  //prior sd for pup proportion (logit scale)
  DATA_MATRIX(Beta_VCsd); //prior covariance matrix for regression parms
  DATA_MATRIX(Beta_VCrn); //prior covariance matrix for regression parms
  DATA_MATRIX(Beta_VCbd); //prior covariance matrix for regression parms
  DATA_MATRIX(Beta_VCrd); //prior covariance matrix for regression parms
  DATA_VECTOR(Melt_i);  // includes 
  DATA_VECTOR(h_mean);  //mean of haulout distributions for each species - use in penalty
  DATA_SCALAR(phi_max); //upper bound for Tweedie phi parameter


  // Parameters 
  PARAMETER_VECTOR(log_N);
  PARAMETER_VECTOR(beta_sd);              // fixed effects on density
  PARAMETER_VECTOR(beta_rn);              // fixed effects on density
  PARAMETER_VECTOR(beta_bd);              // fixed effects on density
  PARAMETER_VECTOR(beta_rd);              // fixed effects on density
  PARAMETER_VECTOR(thin_beta_day);    // for potential time adjustment to estimated detection probs
  PARAMETER_VECTOR(phi_logit);  //tweedie phi  (impose (0,5) bounds)
  PARAMETER_VECTOR(p_logit);  // tweedie power (impose (1,2) bounds)
  PARAMETER_VECTOR(thin_logit_i);         // thinning "parameter" for each surveyed location (assumed MVN on logit scale)
  PARAMETER_VECTOR(MisID_pars);   //confusion matrix estimates (mlogit scale)
  PARAMETER_VECTOR(logit_Pup_prop);
  PARAMETER(beta0_ringed);
  PARAMETER(beta1_ringed);
  //

  // derived sizes
  //int n_b = X_s.row(0).size();
  int n_i = C_i.col(0).size();
  int n_st = X_sd.col(0).size();
  //int n_i_model = Which_counts_model.size();
  int n_obs_types = C_i.row(0).size();
  int n_misID_par = MisID_pars.size();
  int n_zero = Which_zero.size();
  Type zero = 0;
 
  // global stuff
  vector<Type> jnll_comp(5);
  jnll_comp.setZero();
  
  vector<Type> N(n_sp);
  for (int isp = 0; isp<n_sp; isp++)N(isp) = exp(log_N(isp));
  vector<Type> Pup_prop(n_sp);
  for(int isp=0;isp<n_sp;isp++){
    Pup_prop(isp)=1.0/(1.0+exp(-logit_Pup_prop(isp)));
  }
  MVNORM_t<Type>neg_log_density_misID(MisID_Sigma);
  MVNORM_t<Type>neg_log_density_thin(Sigma_logit_thin);
  MVNORM_t<Type>neg_log_density_beta_sd(Beta_VCsd);
  MVNORM_t<Type>neg_log_density_beta_rn(Beta_VCrn);
  MVNORM_t<Type>neg_log_density_beta_bd(Beta_VCbd);
  MVNORM_t<Type>neg_log_density_beta_rd(Beta_VCrd);
  matrix<Type> Thin_i(n_sp,n_i);  
  matrix<Type> Thin_trans(n_sp, n_i);
  vector<Type> Day_effect(n_i);
  vector<Type> h_mean_obs(n_sp);
  Day_effect = X_day * thin_beta_day;
  h_mean_obs = h_mean_obs.setZero();
  Thin_i.setZero();
  for (int isp = 0; isp<(n_sp-1); isp++) {
	  for (int i = 0; i<n_i; i++) {
		  Thin_trans(isp, i) = 1 / (1 + exp(-thin_logit_i(n_i*isp + i) - Day_effect(n_i*isp + i)));
	    if(Ice_s(S_i(i))>0)Thin_i(isp, i) = (P_i(i)/(Ice_s(S_i(i))*A_s(S_i(i))))*Thin_trans(isp, i);
		  h_mean_obs(isp) += Thin_trans(isp, i);
	  }
	  h_mean_obs(isp) = h_mean_obs(isp) / n_i;
  }
  int isp=3;
  for( int i=0;i<n_i;i++){  //ringed
    Thin_trans(isp,i)=1/((1+exp(-(beta0_ringed + beta1_ringed*Melt_i(i))))*(1.0+exp(-thin_logit_i(n_i*isp+i))));//detection * availability
    if(Ice_s(S_i(i))>0)Thin_i(isp,i)=(P_i(i)/(Ice_s(S_i(i))*A_s(S_i(i))))*Thin_trans(isp,i);
  }

  vector<Type> phi(n_obs_types+2);
  vector<Type> power(n_obs_types+2);
  for (int i = 0; i<(n_obs_types+2); i++) {
	  phi(i) = phi_max / (1 + exp(-phi_logit(i)));
	  power(i) = 1.0 + 1 / (1 + exp(-p_logit(i)));
  }
  
  // Transform confusion matrix
  matrix<Type> Psi(n_sp, n_obs_types);
  Psi.fill(0.0);
  for (int ipar = 0; ipar<n_misID_par; ipar++) {
	  Psi(MisID_pos_rows(ipar), MisID_pos_cols(ipar)) = MisID_pars(ipar);
  }
  Type tmp_sum;
  for (int isp = 0; isp<n_sp; isp++) {
	  tmp_sum = 0.0;
	  for (int itype = 0; itype<n_obs_types; itype++) {
		  tmp_sum = tmp_sum + exp(Psi(isp, itype));
	  }
	  for (int itype = 0; itype<n_obs_types; itype++) {
		  Psi(isp, itype) = exp(Psi(isp, itype)) / tmp_sum;
	  }
  }
  //std::cout<<"Psi "<<Psi<<'\n'<<'\n';

  
  // Predicted densities
  matrix<Type> Pi_s(n_sp, n_st);
  matrix<Type> Z_s(n_sp, n_st);
  matrix<Type> E_count_sp(n_sp, n_i);
  matrix<Type> E_count_obs(n_i, n_obs_types);
  matrix<Type> E_count_zero(n_sp, n_0);
  vector<Type> E_count_pup(n_i);
  vector<Type> E_count_nophoto(n_i);
  matrix<Type> linpredZ_s(n_sp,n_st);
  
  Type cur_sum;
	linpredZ_s.row(0) = X_sd * beta_sd;
	linpredZ_s.row(1) = X_rn * beta_rn;
	linpredZ_s.row(2) = X_bd * beta_bd;
	linpredZ_s.row(3) = X_rd * beta_rd;
	
	for(int isp=0;isp<n_sp;isp++){
	  for (int ist = 0; ist<n_st; ist++) {
		  Pi_s(isp, ist) = A_s(ist)*Ice_s(ist)*exp(linpredZ_s(isp,ist));
	  }
    for(int izero = 0; izero<n_zero;izero++){
	    Pi_s(isp,Which_zero(izero))=0.0;
	  }
	  for (int it = 0; it<n_t; it++) {
		  cur_sum = 0.0;
		  for (int is = 0; is<n_s; is++) {
			  cur_sum += Pi_s(isp, it*n_s + is);
		  }
		  for (int is = 0; is<n_s; is++) {
			  Pi_s(isp, it*n_s + is) = Pi_s(isp, it*n_s + is) / cur_sum;
			  Z_s(isp, it*n_s + is) = N(isp)*Pi_s(isp, it*n_s + is);
		  }
	  }
  }
  //std::cout<<"Beta "<<Beta<<"\n";
  //std::cout<<"Z_s "<<Z_s<<"\n";
  //std::cout<<"linpred "<<linpredZ_s<<"\n";
  //std::cout<<"Beta_tmp "<<Beta_tmp<<"\n";
  
  // Probability of counts
  E_count_obs = E_count_obs.setZero();
  E_count_pup = E_count_pup.setZero();
  E_count_zero = E_count_zero.setZero();
  E_count_nophoto = E_count_nophoto.setZero();
  for (int i = 0; i<n_i; i++) {
	  for (int isp = 0; isp<n_sp; isp++) {
		  E_count_sp(isp, i) = Z_s(isp, S_i(i))*(1.0-Pup_prop(isp))*Thin_i(isp, i)*Prop_photo_i(i);
		  for (int itype = 0; itype<n_obs_types; itype++) {
			  E_count_obs(i, itype) += E_count_sp(isp, i)*Psi(isp, itype);
		  }
		  E_count_pup(i) += Z_s(isp, S_i(i))*Pup_prop(isp)*Thin_i(isp, i)*Prop_photo_i(i);
		  E_count_nophoto(i) += Z_s(isp, S_i(i))*Thin_i(isp, i)*(1.0 - Prop_photo_i(i));
	  }
  }
  for (int i = 0; i<n_0; i++) {
    for (int isp = 0; isp<n_sp; isp++) {
      E_count_zero(isp, i) = Z_s(isp, S_i(n_i+i));
    }
  }
  //std::cout<<"E_count_sp: "<<E_count_sp<<"\n";
  
  for (int i = 0; i<n_i; i++) {
	  for (int itype = 0; itype<n_obs_types;itype++) {
		  //if (!isNA(C_i(i,itype)) && E_count_obs(i,itype)>0.0) jnll_comp(0) -= dtweedie(C_i(i, itype), E_count_obs(i, itype), phi(itype), power(itype), true);
		  jnll_comp(0) -= dtweedie(C_i(i, itype), E_count_obs(i, itype), phi(itype), power(itype), true);
	  }
	  //if (!isNA(Pups_i(i)) && E_count_pup(i)>0.0) jnll_comp(0) -= dtweedie(Pups_i(i), E_count_pup(i), phi(n_obs_types), power(n_obs_types), true);
	  jnll_comp(0) -= dtweedie(Pups_i(i), E_count_pup(i), phi(n_obs_types), power(n_obs_types), true);

	  //if (!isNA(Nophoto_i(i)) && E_count_nophoto(i)>0.0) jnll_comp(0) -= dtweedie(Nophoto_i(i), E_count_nophoto(i), phi(n_obs_types+1), power(n_obs_types+1), true);
	  jnll_comp(0) -= dtweedie(Nophoto_i(i), E_count_nophoto(i), phi(n_obs_types + 1), power(n_obs_types + 1), true);
  }
  for (int i = 0; i<n_0; i++){
    for (int isp = 0; isp<n_sp;isp++) {
      jnll_comp(0) -= dtweedie(zero, E_count_zero(isp,i), phi(isp), power(isp), true);
    }
  }
  
  //thinning prior
  //jnll_comp(1) = neg_log_density_thin(thin_logit_i - thin_mu_logit - Day_effect);
  
  jnll_comp(1) += 10000.0 * pow(1.0 - exp(beta0_ringed + beta1_ringed*0.0) / (1.0+exp(beta0_ringed + beta1_ringed*0.0)),2.0);
  
  //for (int ispeff = 0; ispeff<(n_sp * 2); ispeff++) jnll_comp(1) += pow(thin_beta_day(ispeff), 2.0);  //ridge reg
  
  //for (int isp = 0; isp < n_sp; isp++)jnll_comp(1) += 1000000000.0 * pow(h_mean_obs(isp) - h_mean(isp), 2.0);

  
  //misID prior
  //jnll_comp(2) = neg_log_density_misID(MisID_pars - MisID_mu);

  
  //pup proportion prior
  for (int isp = 0; isp < n_sp; isp++) {
	  jnll_comp(3) -= dnorm(logit_Pup_prop(isp), Pup_prop_mu(isp), Pup_prop_sigma(isp), true);
  }
  
  //beta prior - this doesn't work well with the strata covariates
  //jnll_comp(4) = neg_log_density_beta_sd(beta_sd);
  //jnll_comp(4) += neg_log_density_beta_rn(beta_rn);
  //jnll_comp(4) += neg_log_density_beta_bd(beta_bd);
  //jnll_comp(4) += neg_log_density_beta_rd(beta_rd);
  
  // Total objective
  Type jnll = jnll_comp.sum();
  //std::cout<<"jnll_comp: "<<jnll_comp<<"\n";

  // Reporting
  REPORT( Z_s );
  REPORT( N );
  REPORT( beta_sd );
  REPORT( beta_rn );
  REPORT( beta_bd );
  REPORT( beta_rd );
  REPORT(E_count_obs);
  REPORT(E_count_sp);
  REPORT(E_count_pup);
  REPORT(E_count_nophoto);
  REPORT(Pup_prop);
  REPORT(Thin_i);
  REPORT(Psi);
  REPORT(phi);
  REPORT(power);
  REPORT(beta0_ringed);
  REPORT(beta1_ringed);
  
  //REPORT( linpredZ_s );
  REPORT( jnll_comp );
  REPORT( jnll );
  
  ADREPORT( N);
  
  return jnll;
}
