// Stan model using a Gamma prior for concentration
data {
  int<lower=1> K;
  int<lower=1> C;
  array[K, C] int<lower=0> X;
  int<lower=1> m;
}
parameters {
  simplex[C] pi;
  real<lower=0> concentration;
}
model {
  // Priors
  pi ~ dirichlet(rep_vector(1.0, C));
  concentration ~ gamma(2, 0.1); // Gamma Prior

  // Likelihood
  for (k in 1:K) {
    X[k] ~ dirichlet_multinomial(pi * concentration);
  }
}
generated quantities {
  array[C] int y_pred;
  y_pred = dirichlet_multinomial_rng(pi * concentration, m);
}
