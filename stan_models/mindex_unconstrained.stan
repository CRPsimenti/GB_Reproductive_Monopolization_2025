// adapted from https://github.com/ctross/SkewCalc/R/StanModel.R
// by Christof Neumann

functions{
  real Mraw_index(vector r, vector t){
    int N = rows(r);
    
    real R; 
    real T; 
    vector[N] s; 
    real Mraw;                       
    
    R = sum(r);
    T = sum(t);  
    
    for(i in 1:N)
    s[i] = ((r[i]/R)-(t[i]/T))^2;  
    
    Mraw = N * sum(s); 
    return Mraw;            
  }
  
  vector pow2(vector x, real y){
    vector[rows(x)] z;
    for (i in 1:rows(x)) 
    z[i] = x[i]^y;
    return(z);
  }
}

data{
  int N;
  array[N] int r;
  vector[N] t; 
  vector[N] t0;  
}

transformed data{
  int R;
  R = sum(r);
}

parameters{
  vector[N] alpha;
  real<lower=0> xsd;
  real<lower=0> gamma;
}

model{ 
  real T;
  real T_star;
  
  vector[N] t_hat;
  vector[N] t_hat_star;
  
  gamma ~ normal(1, 0.1);

  T = sum(t-t0);
  t_hat = (t-t0)/T;
  
  T_star = sum(pow2(t,gamma) - pow2(t0,gamma));
  t_hat_star = (pow2(t,gamma) - pow2(t0,gamma))/T_star;
  
  alpha ~ normal(t_hat_star, xsd);
  
  xsd ~ exponential(1);
  r ~ multinomial_logit(alpha); 
}

generated quantities{
  array[N] int r_rep;
  real M_raw;
  real M;

  r_rep = multinomial_logit_rng(alpha, R);
  
  { 
    real Bias;
    real T;
    real T_star;
    
    vector[N] t_hat;
    vector[N] t_hat_star;
    
    T = sum(t-t0);
    t_hat = (t-t0)/T;
    
    T_star = sum(pow2(t,gamma) - pow2(t0,gamma));
    t_hat_star = (pow2(t,gamma) - pow2(t0,gamma))/T_star;
    
    M_raw =  Mraw_index(to_vector(multinomial_logit_rng(alpha,R)),t_hat);

    Bias = Mraw_index(to_vector(multinomial_logit_rng(t_hat,R)),t_hat);
    M = M_raw - Bias;
  }
} 
