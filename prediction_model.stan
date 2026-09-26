data{
  int<lower=1> T;    // istanti temporali (SOLO pre)
  int<lower=1> M;    // numero stazioni
  int<lower=1> K;    // numero covariate

  array[T] vector[M] y;            // y[t] = vector[M]
  array[T] matrix[M, K] X;         // X[t] = matrix[M,K]

  matrix[M, M] dist_mat;           // distanze tra stazioni (km)
  real<lower=0> l;                 // range (es. max_dist/3)

  cov_matrix[M] Sigma0;            // cov empirica del pre (sulla y che passi a Stan)
}

transformed data{
  vector[M] s0;
  for (m in 1:M)
    s0[m] = sqrt(Sigma0[m, m] + 1e-12);  // sd marginali empiriche (con eps)
}

parameters{
  real<lower=-1, upper=1> alpha;
  vector[K] beta;

  real<lower=0> sigma_w;
  vector[M] z_w;                   // non-centered per w

  vector<lower=0>[M] tau;          // sd marginali di Sigma
  cholesky_factor_corr[M] L_Omega; // Cholesky correlazioni
}

transformed parameters{
  matrix[M, M] Sw;
  matrix[M, M] Lw;
  vector[M] w;

  for (i in 1:M) {
    for (j in 1:M) {
      Sw[i, j] = square(sigma_w) * exp(-dist_mat[i, j] / l);
    }
    Sw[i, i] += 1e-6; // jitter numerico
  }

  Lw = cholesky_decompose(Sw);
  w  = Lw * z_w;
}

model{
  // ---- Priors ----
  alpha   ~ uniform(-1, 1);
  beta    ~ normal(0, 1);

  sigma_w ~ normal(0, 1);  // se y è standardizzata, ok
  z_w     ~ normal(0, 1);

  // Prior informata da Sigma0 sulle sd marginali (tau)
  // 0.35 = moderatamente informativa; aumenta a 0.5-0.7 se vuoi più weak.
  tau ~ lognormal(log(s0), 0.35);

  L_Omega ~ lkj_corr_cholesky(2);

  // ---- Likelihood ----
  {
    matrix[M, M] L_Sigma;
    L_Sigma = diag_pre_multiply(tau, L_Omega);

    y[1] ~ multi_normal_cholesky(X[1] * beta + w, L_Sigma);

    for (t in 2:T) {
      vector[M] mu;
      mu = alpha * y[t-1] + X[t] * beta + w;
      y[t] ~ multi_normal_cholesky(mu, L_Sigma);
    }
  }
}


