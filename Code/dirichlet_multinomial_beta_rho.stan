data {
  int<lower=1> K;
  int<lower=1> C;
  array[K, C] int<lower=0> X;
  int<lower=1> m;
}
parameters {
  simplex[C] pi;
  // Model rho (ICC) instead of concentration directly
  real<lower=0, upper=1> rho; 
}
transformed parameters {
  // Transform rho back to concentration for the DM function
  // Avoid division by zero: if rho is 0, concentration is infinite (handled by large number)
  real concentration = (rho < 1e-5) ? 1e5 : (1 - rho) / rho;
}
model {
  pi ~ dirichlet(rep_vector(1.0, C));
  
  // Shrinkage Prior: Puts mass near 0 (Multinomial), penalizes complexity (Overdispersion)
  // Beta(1, 10) is conservative (favors Multinomial). 
  // Beta(2, 10) peaks slightly away from 0 but still low.
  rho ~ beta(1, 10); 

  for (k in 1:K) {
    X[k] ~ dirichlet_multinomial(pi * concentration);
  }
}
generated quantities {
  array[C] int y_pred;
  y_pred = dirichlet_multinomial_rng(pi * concentration, m);
}
