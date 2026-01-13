//
// This Stan program defines a two-component Gaussian mixture model,
// where each  logit transformed WL ratio (WL_odd_ratios = WLratio/(1-WLratio))is modeled as coming 
// from either a herbaceous-derived or woody-derived charcoal distribution.
// The mixture proportions are defined by the simplex 'theta',
// the component means are 'muH' and 'muW',
// and the standard deviations are in 'sigma'.
// Log-likelihood is computed using log_mix.
//
data {
  int<lower=1> N;                 // Number of observations
  vector<lower = 0>[N] WL_odd_ratios;  // Observed data (logit transformed WLratio)
}
parameters {
  real lambda; // parameter for BoxCox
  simplex[2] theta;   // Mixing proportions for the two components (herbaceous vs woody)
  
  // Mean of the herbaceous charcoal component, constrained to be ≤ -0.4
  real<upper=-0.5> muH;
  
  // Mean of the woody charcoal component, constrained to be ≥ 0.2
  real<lower=0.4> muW;

  // Standard deviations for both components, constrained between 0.4 and 1.1
  vector<lower=0.4, upper=1.1>[2] sigma;
  
}
transformed parameters {
  vector[2] mu;             // Vector of component means
  mu[1] = muH;          // Index 1: grass-derived mean
  mu[2] = muW;          // Index 2: tree-derived mean
  vector[N] BoxCoxlogitWL;
  if (abs(lambda) < 1e-6)
    BoxCoxlogitWL = log(WL_odd_ratios);
  else
    BoxCoxlogitWL = (pow(WL_odd_ratios, lambda) - 1) / lambda;

}
model {
// Prior: weakly informative Dirichlet for mixing proportions
  theta ~ dirichlet(rep_vector(1.0, 2));
  
// Priors on component means 
  muH ~ normal(-1.7, 0.3); 
  muW ~ normal(0.7, 0.3); 
  
// Prior on standard deviations
  sigma ~ normal(0.7, 0.2);
// Prior on lambda (Box-Cox logit parameter)
  lambda ~ normal(-0.07, 0.3);
// Mixture likelihood using log_mix for numerical stability
  for (n in 1:N) {
    target += log_mix(theta[1],
                      normal_lpdf(BoxCoxlogitWL[n] | mu[1], sigma[1]),
                      normal_lpdf(BoxCoxlogitWL[n] | mu[2], sigma[2]));
  }
// Jacobian term linked to BoxCox variable transform
  target += (lambda - 1) * sum(log(WL_odd_ratios));
}
