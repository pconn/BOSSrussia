// Space time 
#include <TMB.hpp>
//#include <atomic_math.hpp>  //for D_lgamma


/** Precision matrix for the anisotropic case, eqn (20) in Lindgren et al. (2011) */    
namespace R_inla_generalized {
using namespace Eigen;
using namespace tmbutils;
using namespace R_inla;

template<class Type>
  SparseMatrix<Type> Q_spde_generalized(spde_t<Type> spde, Type kappa, int alpha=2){
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  	
  if( alpha==1 ) return kappa_pow2*spde.M0 + spde.M1;
  if( alpha==2 ) return kappa_pow4*spde.M0 + Type(2.0)*kappa_pow2*spde.M1 + spde.M2;
}

template<class Type>
  SparseMatrix<Type> Q_spde_generalized(spde_aniso_t<Type> spde, Type kappa, matrix<Type> H, int alpha=2){

  int i;
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  
  int n_s = spde.n_s;
  int n_tri = spde.n_tri;
  vector<Type> Tri_Area = spde.Tri_Area;
  matrix<Type> E0 = spde.E0;
  matrix<Type> E1 = spde.E1;
  matrix<Type> E2 = spde.E2;
  matrix<int> TV = spde.TV;
  SparseMatrix<Type> G0 = spde.G0;
  SparseMatrix<Type> G0_inv = spde.G0_inv;
	  	  
  //Type H_trace = H(0,0)+H(1,1);
  //Type H_det = H(0,0)*H(1,1)-H(0,1)*H(1,0);
  SparseMatrix<Type> G1_aniso(n_s,n_s); 
  SparseMatrix<Type> G2_aniso(n_s,n_s); 
  // Calculate adjugate of H
  matrix<Type> adj_H(2,2);
  adj_H(0,0) = H(1,1);
  adj_H(0,1) = -1 * H(0,1);
  adj_H(1,0) = -1 * H(1,0);
  adj_H(1,1) = H(0,0);
  // Calculate new SPDE matrices

  // Calculate G1 - pt. 1
  array<Type> Gtmp(n_tri,3,3);
  for(i=0; i<n_tri; i++){    
    // 1st line: E0(i,) %*% adjH %*% t(E0(i,)), etc.    
    Gtmp(i,0,0) = (E0(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E0(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,1) = (E1(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E1(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,2) = (E2(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E2(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,1) = (E1(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E1(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,2) = (E2(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E2(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,2,2) = (E2(i,0)*(E2(i,0)*adj_H(0,0)+E2(i,1)*adj_H(1,0)) + E2(i,1)*(E2(i,0)*adj_H(0,1)+E2(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
  }
  // Calculate G1 - pt. 2
  for(i=0; i<n_tri; i++){
    G1_aniso.coeffRef(TV(i,1),TV(i,0)) = G1_aniso.coeffRef(TV(i,1),TV(i,0)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i,0),TV(i,1)) = G1_aniso.coeffRef(TV(i,0),TV(i,1)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i,2),TV(i,1)) = G1_aniso.coeffRef(TV(i,2),TV(i,1)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i,1),TV(i,2)) = G1_aniso.coeffRef(TV(i,1),TV(i,2)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i,2),TV(i,0)) = G1_aniso.coeffRef(TV(i,2),TV(i,0)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i,0),TV(i,2)) = G1_aniso.coeffRef(TV(i,0),TV(i,2)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i,0),TV(i,0)) = G1_aniso.coeffRef(TV(i,0),TV(i,0)) + (Gtmp(i,0,0));  
    G1_aniso.coeffRef(TV(i,1),TV(i,1)) = G1_aniso.coeffRef(TV(i,1),TV(i,1)) + (Gtmp(i,1,1));  
    G1_aniso.coeffRef(TV(i,2),TV(i,2)) = G1_aniso.coeffRef(TV(i,2),TV(i,2)) + (Gtmp(i,2,2));  
  }
  G2_aniso = G1_aniso * G0_inv * G1_aniso; 

  if( alpha==1 ) return kappa_pow2*G0 + G1_aniso;
  if( alpha==2 ) return kappa_pow4*G0 + Type(2.0)*kappa_pow2*G1_aniso + G2_aniso;
}
} // end namespace R_inla

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}

template<class Type>
Type dbern(Type x, Type prob, int give_log=1){
  Type logres;
  if( x==0 ) logres = log( 1-prob );
  if( x==1 ) logres = log( prob );
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace R_inla_generalized;
  using namespace Eigen;
  using namespace density;
  
  // options vec
  //DATA_FACTOR( Options_vec );
  // Slot 0: compute SE? 
  
  // Data
  DATA_MATRIX(C_i);       	// Matrix of responses (counts) of each species at each sampled location 
  DATA_VECTOR(Pups_i);      // vector holding observations of pups in each sampled cell
  DATA_VECTOR( P_i );      // Proportion of survey unit that is surveyed for each observation i
  DATA_VECTOR( A_s );      // Relative area of each survey unit s (can set = 1.0 if all the same size)
  DATA_IVECTOR(S_i); // Site-time index for each sample
  DATA_MATRIX( X_s );  //design matrix for fixed effects

  DATA_VECTOR(thin_mu_logit);
  DATA_SPARSE_MATRIX(Sigma_logit_thin);
  //DATA_MATRIX(X_day);  // design matrix for extra day and day^2 effects
  DATA_VECTOR(MisID_mu);
  DATA_MATRIX(MisID_Sigma);  //variance covariance matrix for MisID_pars on mlogit scale
  DATA_IVECTOR(MisID_pos_rows);   //row indices for positive entries of MisID matrix
  DATA_IVECTOR(MisID_pos_cols);   //ibid, columns
  DATA_INTEGER(n_s); //number of cells
  DATA_INTEGER(n_sp); //number of species
  DATA_VECTOR(Pup_prop_mu);  //prior mean (logit scale) for pup proportion
  DATA_VECTOR(Pup_prop_sigma);  //prior sd for pup proportion (logit scale)


  // spatial covariance object
  DATA_STRUCT(spde, spde_t);

  // Parameters 
  PARAMETER_MATRIX(beta);              // fixed effects on density
  PARAMETER_VECTOR(logtau_z);      
  PARAMETER_VECTOR(logkappa_z);
  PARAMETER_VECTOR(phi_log);  //tweedie phi
  PARAMETER_VECTOR(p_logit);  // tweedie power
  PARAMETER_VECTOR(thin_logit_i);         // thinning "parameter" for each surveyed location (assumed MVN on logit scale)
  PARAMETER_VECTOR(MisID_pars);   //confusion matrix estimates (mlogit scale)
  PARAMETER_VECTOR(logit_Pup_prop);
  //
  //// Random effects
  PARAMETER_MATRIX( Etainput_s );           // Spatial variation in abundance

  // derived sizes
  int n_b = X_s.col(0).size();
  int n_i = C_i.col(0).size();
  //int n_i_model = Which_counts_model.size();
  int n_obs_types = C_i.row(0).size();
  int n_misID_par = MisID_pars.size();
  int n_re = Etainput_s.row(0).size();
 
  // global stuff
  vector<Type> jnll_comp(5);
  jnll_comp.setZero();
  vector<Type> MargSD_z(n_sp);
  vector<Type> Range_z(n_sp);
  for( int z=0; z<n_sp; z++){
    MargSD_z(z) = 1 / sqrt(4*M_PI) / exp(logtau_z(z)) / exp(logkappa_z(z));
    Range_z(z) = sqrt(8) / exp( logkappa_z(z) );
  }
  vector<Type> Pup_prop(n_sp);
  for(int isp=0;isp<n_sp;isp++){
    Pup_prop(isp)=1.0/(1.0+exp(-logit_Pup_prop(isp)));
  }
  MVNORM_t<Type>neg_log_density_misID(MisID_Sigma);
  MVNORM_t<Type>neg_log_density_thin(Sigma_logit_thin);
  matrix<Type> Thin_i(n_sp,n_i);  
  matrix<Type> Thin_trans(n_sp, n_i);
  for (int i = 0; i<n_i; i++) {
	  for (int isp = 0; isp<n_sp; isp++) {
		  //Thin_trans(isp,i)=1/(1+exp(-thin_mu_logit(n_i*isp+i)-Day_effect(n_i*isp+i)));
		  Thin_trans(isp, i) = 1 / (1 + exp(-thin_logit_i(n_i*isp + i)));
		  Thin_i(isp, i) = P_i(i)*A_s(S_i(i))*Thin_trans(isp, i);
	  }
  }
  vector<Type> phi(n_obs_types+1);
  vector<Type> power(n_obs_types+1);
  for (int i = 0; i<(n_obs_types+1); i++) {
	  phi(i) = exp(phi_log(i));
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

  // Transform random effects
  matrix<Type> Eta_s( n_sp,n_s );
  vector<Type> Tau_z(n_sp);
  for (int isp = 0; isp < n_sp; isp++) {
	  Tau_z(isp) = exp(logtau_z(isp));
	  for (int ire = 0; ire < n_s; ire++) {
		  Eta_s(isp, ire) = Etainput_s(isp, ire) / Tau_z(isp);
	  }
  }

  //random effects priors
  SparseMatrix<Type> Q;
  vector<Type> Cur_eta(n_re);
  for (int isp = 0; isp < n_sp; isp++) {
	  Q = Q_spde_generalized(spde, exp(logkappa_z(isp)), 2);  //default smoothing is alpha = 2
     for (int ieta = 0; ieta<n_re; ieta++) {
	    Cur_eta(ieta) = Etainput_s(isp, ieta);
     }
     jnll_comp(1) += GMRF(Q)(Cur_eta);
  }
  
  // Predicted densities
  matrix<Type> Z_s(n_sp, n_s);
  matrix<Type> E_count_sp(n_sp, n_i);
  matrix<Type> E_count_obs(n_i, n_obs_types);
  vector<Type> E_count_pup(n_i);
  vector<Type> linpredZ_s(n_s);
  vector<Type> Beta_tmp(n_b);
  for (int isp = 0; isp<n_sp; isp++) {
	  Beta_tmp = beta.row(isp);
	  linpredZ_s = X_s * Beta_tmp;
	  for (int is = 0; is < n_s; is++) {
		  Z_s(isp, is) = exp(linpredZ_s(is)+Eta_s(isp, is));
	  }
  }

  // Probability of counts
  E_count_obs = E_count_obs.setZero();
  E_count_pup = E_count_pup.setZero();
  for (int i = 0; i<n_i; i++) {
	  for (int isp = 0; isp<n_sp; isp++) {
		  E_count_sp(isp, i) = Z_s(isp, S_i(i))*(1.0-Pup_prop(isp))*Thin_i(isp, i);
		  for (int itype = 0; itype<n_obs_types; itype++) {
			  E_count_obs(i, itype) += E_count_sp(isp, i)*Psi(isp, itype);
		  }
		  E_count_pup(i) += Z_s(isp, S_i(i))*Pup_prop(isp)*Thin_i(isp, i);
	  }
  }
  for (int i = 0; i<n_i; i++) {
	  for (int itype = 0; itype<n_obs_types;itype++) {
		  if (!isNA(C_i(i,itype))) jnll_comp(0) -= dtweedie(C_i(i, itype), E_count_obs(i, itype), phi(itype), power(itype), true);
	  }
	  if (!isNA(Pups_i(i))) jnll_comp(0) -= dtweedie(Pups_i(i), E_count_pup(i), phi(n_obs_types), power(n_obs_types), true);
  }

  //thinning prior
  jnll_comp(2) = neg_log_density_thin(thin_logit_i - thin_mu_logit);

  //misID prior
  jnll_comp(3) = neg_log_density_misID(MisID_pars - MisID_mu);

  //pup proportion prior
  for (int isp = 0; isp < n_sp; isp++) {
	  jnll_comp(4) -= dnorm(logit_Pup_prop(isp), Pup_prop_mu(isp), Pup_prop_sigma(isp),true);
  }

  //Type total_abundance = Z_s.sum();
  vector<Type> total_abundance = Z_s.rowwise().sum();

  // Total objective
  Type jnll = jnll_comp.sum();

  // Reporting
  REPORT( Z_s );
  REPORT( total_abundance );
  REPORT( Range_z );
  REPORT( MargSD_z );
  REPORT( beta );
  REPORT( Eta_s );
  REPORT(E_count_obs);
  REPORT(Pup_prop);
  REPORT(Thin_i);
  REPORT(Z_s);
  REPORT(Psi);
  REPORT(phi);
  REPORT(power);
  //REPORT( linpredZ_s );
  REPORT( jnll_comp );
  REPORT( jnll );
  
  // Bias correction output
  ADREPORT( total_abundance);

  return jnll;
}
