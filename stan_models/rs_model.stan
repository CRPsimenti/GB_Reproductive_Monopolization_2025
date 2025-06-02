functions{
  real rating_priors(
                   int mig_mode, // 1=inherit, 2=average, 3=bottom, 4=top
                   array[] int startindex,
                   vector start_ratings,
                   vector k,
                   int N, // n interaction
                   int n_ind, // n individuals
                   array[] int winner, // winner index
                   array[] int loser, // loser index
                   array[] int group, // winner index
                   matrix presence,
                   array[] int group_size
  ) {
    vector[n_ind] out = rep_vector(0.0, n_ind);
    if (mig_mode == 1) {
      for (i in 1:n_ind) {
        if (startindex[i] == 1) {
          out[i] = 0.0;
        }
        if (startindex[i] > 1) {
          vector[n_ind] current_ratings = start_ratings;
          vector[startindex[i]] result;
          // real result;
          for (j in 1:startindex[i]) {
            array[group_size[j]] int sel = group_index_foo(presence[j, ], group_size[j], n_ind);

            current_ratings[sel] = current_ratings[sel] - dot_product(presence[j, sel], current_ratings[sel]) / sum(presence[j, sel]);
            // likelihood:
            result[j] = 1/(1 + exp((current_ratings[loser[j]] - current_ratings[winner[j]])));
            // update:
            real toAdd = (1 - result[j]) * k[group[j]];
            current_ratings[winner[j]] = current_ratings[winner[j]] + toAdd;
            current_ratings[loser[j]] = current_ratings[loser[j]] - toAdd;
          }
          out[i] = current_ratings[i - 1];
        }
      }
    }
    if (mig_mode == 2) {
      for (i in 1:n_ind) {
        out[i] = 0.0;
      }
    }
    // bottom
    if (mig_mode == 3) {
      for (i in 1:n_ind) {
        if (startindex[i] == 1) {
          out[i] = 0.0;
        }
        if (startindex[i] > 1) {
          out[i] = -1.5;
        }
      }
    }
    // top
    if (mig_mode == 4) {
      for (i in 1:n_ind) {
        if (startindex[i] == 1) {
          out[i] = 0.0;
        }
        if (startindex[i] > 1) {
          out[i] = 1.5;
        }
      }
    }
    return(normal_lpdf(start_ratings | out, 1.0));
  }

  // function returns an index vector of individuals that were present for a given interaction
  // i.e. the index of group members
  // taken from a row of the presence matrix
  // in R this could be avoided via a rugged array (a list)
  array[] int group_index_foo(row_vector prow, int gs, int n_ind) {
    array[gs] int g_index = rep_array(0, gs);
    int cnt = 0;
    for (i in 1:n_ind) {
      if (prow[i] == 1) {
        cnt = cnt + 1;
        g_index[cnt] = i;
      }
    }
    return(g_index);
  }  
  
  // obtain winning probabilities for a given k and set of start ratings
  // code is largely based on Goffe et al 2018
  // result is a vector with winning probabilities (of the actual winner) for each interaction
  vector win_probs_from_seq(vector start_ratings,
                            vector k,
                            int N, // n interaction
                            int n_ind, // n individuals
                            // int n_grp, // n individuals
                            array[] int winner, // winner index
                            array[] int loser, // loser index
                            array[] int group, // index for group
                            matrix presence, // presence matrices for groups by interaction
                            array[] int group_size
                            ) {
    vector[N] result;
    real toAdd;
    vector[n_ind] current_ratings = start_ratings;
    for (i in 1:N) {
      array[group_size[i]] int sel = group_index_foo(presence[i, ], group_size[i], n_ind);
      current_ratings[sel] = current_ratings[sel] - dot_product(presence[i, sel], current_ratings[sel]) / sum(presence[i, sel]);
      result[i] = 1/(1 + exp((current_ratings[loser[i]] - current_ratings[winner[i]])));
      toAdd = (1 - result[i]) * k[group[i]];
      current_ratings[winner[i]] = current_ratings[winner[i]] + toAdd;
      current_ratings[loser[i]] = current_ratings[loser[i]] - toAdd;
    }
    return result;
  }
  
  // specific *rating* of an individual at a given time point
  // time point is the index position of an interaction in the sequence
  // i.e. the rating for a given observation of the RS data
  vector up_to(vector start_ratings,
               vector k,
               int N, // n interaction
               int n_ind_elo, // n individuals
               array[] int winner, // winner index
               array[] int loser, // loser index
               array[] int group, // index for group
               matrix presence, // presence matrices for groups by interaction
               array[] int group_size,
               array[] int stop_here, // the index position of the interaction sequence that corresponds to a data point in the RS table
               int n_model_obs,
               array[] int model_male_index_elo
                            ) {
                              
    matrix[N, n_ind_elo] ratmat; // interaction by individual matrix (rating matrix)
    vector[n_model_obs] out;
    vector[N] result;
    vector[n_ind_elo] current_ratings = start_ratings;
    for (i in 1:N) { // iterate through interactions
      real toAdd;
      array[group_size[i]] int sel = group_index_foo(presence[i, ], group_size[i], n_ind_elo);
      current_ratings[sel] = current_ratings[sel] - dot_product(presence[i, sel], current_ratings[sel]) / sum(presence[i, sel]);
      result[i] = 1/(1 + exp((current_ratings[loser[i]] - current_ratings[winner[i]])));
      toAdd = (1 - result[i]) * k[group[i]];
      current_ratings[winner[i]] = current_ratings[winner[i]] + toAdd;
      current_ratings[loser[i]] = current_ratings[loser[i]] - toAdd;
      ratmat[i, ] = to_row_vector(current_ratings);
    }
    // from the ratmat: extract the relevant ratings into the obs-level ratings vector
    for (j in 1:n_model_obs) {
      out[j] = ratmat[stop_here[j], model_male_index_elo[j]];
    }
    
    return out;
  } 
  
  matrix up_to_ratmat(vector start_ratings,
               vector k,
               int N, // n interaction
               int n_ind_elo, // n individuals
               array[] int winner, // winner index
               array[] int loser, // loser index
               array[] int group, // index for group
               matrix presence, // presence matrices for groups by interaction
               array[] int group_size,
               array[] int stop_here, // the index position of the interaction sequence that corresponds to a data point in the RS table
               int n_model_obs,
               array[] int model_male_index_elo
                            ) {
                              
    matrix[N, n_ind_elo] ratmat; // interaction by individual matrix (rating matrix)
    vector[n_model_obs] out;
    vector[N] result;
    vector[n_ind_elo] current_ratings = start_ratings;
    for (i in 1:N) { // iterate through interactions
      real toAdd;
      array[group_size[i]] int sel = group_index_foo(presence[i, ], group_size[i], n_ind_elo);
      current_ratings[sel] = current_ratings[sel] - dot_product(presence[i, sel], current_ratings[sel]) / sum(presence[i, sel]);
      result[i] = 1/(1 + exp((current_ratings[loser[i]] - current_ratings[winner[i]])));
      toAdd = (1 - result[i]) * k[group[i]];
      current_ratings[winner[i]] = current_ratings[winner[i]] + toAdd;
      current_ratings[loser[i]] = current_ratings[loser[i]] - toAdd;
      ratmat[i, ] = to_row_vector(current_ratings);
    }
    // from the ratmat: extract the relevant ratings into the obs-level ratings vector
    for (j in 1:n_model_obs) {
      out[j] = ratmat[stop_here[j], model_male_index_elo[j]];
    }
    return ratmat;
  }
}

data {
  int do_model; // allow prior-only

  int<lower=1> n_ind_elo; // number of individuals for Elo, wide format (= migrants occur multiple times)
  int<lower=1> n_int; // number of interactions
  int<lower=1> n_party; // number of parties
  array[n_int] int<lower=1,upper=n_ind_elo> winner_index; // winner's index
  array[n_int] int<lower=1,upper=n_ind_elo> loser_index; // losers's index
  array[n_int] int<lower=1,upper=n_party> party_index; // party's index
  
  array[n_ind_elo] int<lower=1,upper=n_int> start_index; // index positions (wrt interactions) at which individuals switch groups

  int<lower=0,upper=4> mig_mode; //0='classic', 1 = inherit; 2 = average; 3 = bottom; 4 = top

  matrix[n_int, n_ind_elo] presence; // presence matrix (again, migrants occur per group/party)
  
  array[n_int] int<lower=1> party_size; // group/party size for each interaction
  
  // model part (RS data)
  int<lower=1> n_model_obs; // number of observations for actual model (RS data)
  array[n_model_obs] int<lower=0> model_fems; // number of females (the response)
  array[n_model_obs] int<lower=0> rating_index; // the index to match rating sequence location to model data
  array[n_model_obs] int<lower=0> model_male_index_elo;
  int<lower=1> n_model_males; // number of males for actual model (RS data) (might be smaller than males for Elo)
  array[n_model_obs] int<lower=0> model_male_index_per_obs;
  int<lower=1> n_model_years; // number of years for actual model (RS data)
  array[n_model_obs] int<lower=0> model_year_index_per_obs;
  int<lower=1> n_model_parties; // number of parties for actual model (RS data)
  array[n_model_obs] int<lower=0> model_party_index_per_obs;
  
  vector[n_model_obs] nonprime; // indicator for nonprime 
}

transformed data {
  // create the response variable for Elo model, which is simply a 'vector' of 1's
  array[n_int] int<lower=0> y = rep_array(1, n_int);
}

parameters {
  vector[n_ind_elo] EloStart;
  vector<lower=0.0>[n_party] k;
  
  // for actual model part
  real intercept;
  real elo_slope;
  real ageclass_slope;
  real ageclass_elo_interaction;
  
  vector[n_model_obs] olre_blups_z;
  real<lower=0.0> olre_sd;

  // correlated random slopes
  // males: elo and age (but not the interaction!)
  matrix[3, n_model_males] blups_mat_males_z; // is transposed!
  cholesky_factor_corr[3] chol_males; // cholesky factor of correlation matrix for males
  vector<lower=0>[3] sds_males;
  // years: elo, age, and interaction
  matrix[4, n_model_years] blups_mat_year_z; // is transposed!
  cholesky_factor_corr[4] chol_year; // cholesky factor of correlation matrix for year
  vector<lower=0>[4] sds_year;
  // party: elo, age, and interaction
  matrix[4, n_model_parties] blups_mat_party_z; // is transposed!
  cholesky_factor_corr[4] chol_party; // cholesky factor of correlation matrix for party
  vector<lower=0>[4] sds_party;

}
transformed parameters {
  vector[n_model_obs] elo_predictor;
  matrix[n_model_males, 3] blups_mat_males_actual;
  matrix[n_model_years, 4] blups_mat_year_actual;
  matrix[n_model_parties, 4] blups_mat_party_actual;
  vector[n_model_obs] olre_blups_actual;

  blups_mat_males_actual = transpose(diag_pre_multiply(sds_males, chol_males) * blups_mat_males_z);
  blups_mat_year_actual = transpose(diag_pre_multiply(sds_year, chol_year) * blups_mat_year_z);
  blups_mat_party_actual = transpose(diag_pre_multiply(sds_party, chol_party) * blups_mat_party_z);
  olre_blups_actual = olre_blups_z * olre_sd;

  elo_predictor = up_to(EloStart, k, n_int, n_ind_elo, winner_index, loser_index, party_index, presence, party_size, rating_index, n_model_obs, model_male_index_elo);
}

model {
  vector[n_model_obs] mu_count = rep_vector(0.0, n_model_obs);
  vector[n_model_obs] mu_intercept;
  vector[n_model_obs] mu_elo_slope;
  vector[n_model_obs] mu_ageclass_slope;
  vector[n_model_obs] mu_interaction_slope;
  
  mu_intercept = intercept + 
                 olre_blups_actual +
                 blups_mat_males_actual[model_male_index_per_obs, 1] + 
                 blups_mat_year_actual[model_year_index_per_obs, 1] + 
                 blups_mat_party_actual[model_party_index_per_obs, 1];
  mu_elo_slope = elo_slope + 
                 blups_mat_males_actual[model_male_index_per_obs, 2] + 
                 blups_mat_year_actual[model_year_index_per_obs, 2] + 
                 blups_mat_party_actual[model_party_index_per_obs, 2];
  mu_ageclass_slope = ageclass_slope + 
                      blups_mat_males_actual[model_male_index_per_obs, 3] + 
                      blups_mat_year_actual[model_year_index_per_obs, 3] + 
                      blups_mat_party_actual[model_party_index_per_obs, 3];
  mu_interaction_slope = ageclass_elo_interaction + 
                         blups_mat_year_actual[model_year_index_per_obs, 4] + 
                         blups_mat_party_actual[model_party_index_per_obs, 4];

  mu_count = mu_count + mu_intercept +
             mu_elo_slope .* elo_predictor +
             mu_ageclass_slope .* nonprime +
             mu_interaction_slope .* (elo_predictor .* nonprime);
  if (do_model) {
    target += poisson_lpmf(model_fems | exp(mu_count));
  }
  target += bernoulli_lpmf(y | win_probs_from_seq(EloStart, k, n_int, n_ind_elo, winner_index, loser_index, party_index, presence, party_size));
  
  target += exponential_lpdf(k | 2);
  if (mig_mode == 0) {
    target += normal_lpdf(EloStart | 0, 1);
  }
  if (mig_mode > 0) {
    target += rating_priors(mig_mode, start_index, EloStart, k, n_int, n_ind_elo, winner_index, loser_index, party_index, presence, party_size);
  }
  
  target += normal_lpdf(intercept | 0, 1);
  target += normal_lpdf(elo_slope | 0, 1);
  target += normal_lpdf(ageclass_slope | 0, 1);
  target += normal_lpdf(ageclass_elo_interaction | 0, 1);
  target += normal_lpdf(olre_blups_z | 0, 1);
  target += exponential_lpdf(olre_sd | 5);

  // blups for random intercepts and elo slope for model males
  target += exponential_lpdf(sds_males | 5);
  for (i in 1:3) {
      target += normal_lpdf(blups_mat_males_z[i, ] | 0, 1);
  }
  target += lkj_corr_cholesky_lpdf(chol_males | 2);

  target += exponential_lpdf(sds_year | 5);
  for (i in 1:4) {
    target += normal_lpdf(blups_mat_year_z[i, ] | 0, 1);
  }
  target += lkj_corr_cholesky_lpdf(chol_year | 2);
  
  target += exponential_lpdf(sds_party | 5);
  for (i in 1:4) {
    target += normal_lpdf(blups_mat_party_z[i, ] | 0, 1);
  }
  target += lkj_corr_cholesky_lpdf(chol_party | 2);
}

generated quantities {
  array[n_model_obs] int<lower=0> femcount_rep = rep_array(0, n_model_obs);
  // predicted elos for model observations
  vector[n_model_obs] elo_pred_rep;
  
  // correlations between varying effects
  corr_matrix[3] cor_mat_males = multiply_lower_tri_self_transpose(chol_males);
  vector<lower=-1, upper=1>[3] cor_males;
  corr_matrix[4] cor_mat_year = multiply_lower_tri_self_transpose(chol_year);
  vector<lower=-1, upper=1>[6] cor_year;
  corr_matrix[4] cor_mat_party = multiply_lower_tri_self_transpose(chol_party);
  vector<lower=-1, upper=1>[6] cor_party;
  
  elo_pred_rep = up_to(EloStart, k, n_int, n_ind_elo, winner_index, loser_index, party_index, presence, party_size, rating_index, n_model_obs, model_male_index_elo);
  
  for (i in 1 : 3) {
    for (j in 1 : (i - 1)) {
      cor_males[choose(i - 1, 2) + j] = cor_mat_males[j, i];
    }
  }
  for (i in 1 : 4) {
    for (j in 1 : (i - 1)) {
      cor_year[choose(i - 1, 2) + j] = cor_mat_year[j, i];
    }
  }
  for (i in 1 : 4) {
    for (j in 1 : (i - 1)) {
      cor_party[choose(i - 1, 2) + j] = cor_mat_party[j, i];
    }
  }

  if (do_model) {
  {
    vector[n_model_obs] mu_count_sim = rep_vector(0.0, n_model_obs);
    vector[n_model_obs] mu_intercept_sim;
    vector[n_model_obs] mu_elo_slope_sim;
    vector[n_model_obs] mu_ageclass_slope_sim;
    vector[n_model_obs] mu_interaction_slope_sim;
    
    mu_intercept_sim = intercept + 
                   olre_blups_actual +
                   blups_mat_males_actual[model_male_index_per_obs, 1] + 
                   blups_mat_year_actual[model_year_index_per_obs, 1] + 
                   blups_mat_party_actual[model_party_index_per_obs, 1];
    mu_elo_slope_sim = elo_slope + 
                   blups_mat_males_actual[model_male_index_per_obs, 2] + 
                   blups_mat_year_actual[model_year_index_per_obs, 2] + 
                   blups_mat_party_actual[model_party_index_per_obs, 2];
    mu_ageclass_slope_sim = ageclass_slope + 
                        blups_mat_males_actual[model_male_index_per_obs, 3] + 
                        blups_mat_year_actual[model_year_index_per_obs, 3] + 
                        blups_mat_party_actual[model_party_index_per_obs, 3];
    mu_interaction_slope_sim = ageclass_elo_interaction + 
                           blups_mat_year_actual[model_year_index_per_obs, 4] + 
                           blups_mat_party_actual[model_party_index_per_obs, 4];
  
    
    mu_count_sim = mu_count_sim + mu_intercept_sim +
               mu_elo_slope_sim .* elo_predictor +
               mu_ageclass_slope_sim .* nonprime +
               mu_interaction_slope_sim .* (elo_predictor .* nonprime);
    femcount_rep = poisson_rng(exp(mu_count_sim));
  }
  } 
  if (do_model == 0) {
    vector[n_model_obs] mu_count_sim = rep_vector(-50.0, n_model_obs);
    femcount_rep = poisson_rng(exp(mu_count_sim));
  }
}
