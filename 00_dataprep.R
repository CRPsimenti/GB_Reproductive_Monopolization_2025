# read raw data and prepare data sets
# this script doesn't require any additional packages installed

# read data ----
# dominance interactions
intdata <- read.csv("data/domdata.csv") 
# reproductive success data 
demdata <- read.csv("data/rs_data.csv") 
# presence matrix (binary present in party on day of interaction per day)
presence <- as.matrix(read.csv("data/presence_matrix.csv", check.names = FALSE))


# create start index for each wide id (when in the interaction sequence could they occur the first time?) ----
# 'wide id': ID code with party id attached (because party association changed for some males)
start_index <- integer(length = ncol(presence)) - 99
names(start_index) <- colnames(presence)

yy <- unique(unlist(lapply(strsplit(colnames(presence), "@"), function(x)x[1]))) # males in pmat

for (i in yy) {
  aux <- presence[, grep(i, colnames(presence)), drop = FALSE]
  if (ncol(aux) == 1) {
    start_index[colnames(aux)] <- 1
  } else {
    if (ncol(aux) == 2) {
      start_index[colnames(aux)[1]] <- 1
      start_index[colnames(aux)[2]] <- min(which(aux[, 2] == 1))
    } else {
      if (ncol(aux) == 3) {
        start_index[colnames(aux)[1]] <- 1
        start_index[colnames(aux)[2]] <- min(which(aux[, 2] == 1))
        start_index[colnames(aux)[3]] <- min(which(aux[, 3] == 1))
      } else {
        stop()
      }
    }
  }
}


# combine to standat ------------
# ids matching wide presence data (= ids as character vector)
winner <- sapply(seq_len(nrow(intdata)), function(x){
  aux <- presence[x, ]
  aux <- aux[aux == 1]
  names(aux)[grepl(paste0(intdata$winner[x], "@"), names(aux))]
})
loser <- sapply(seq_len(nrow(intdata)), function(x){
  aux <- presence[x, ]
  aux <- aux[aux == 1]
  names(aux)[grepl(paste0(intdata$loser[x], "@"), names(aux))]
})

group <- intdata$party


# create index variables for ids and parties

all_ids <- colnames(presence)
winner_index <- as.numeric(sapply(winner, function(x)which(all_ids == x)))
loser_index <- as.numeric(sapply(loser, function(x)which(all_ids == x)))
names(winner_index) <- winner
names(loser_index) <- loser
n_ind <- length(all_ids)
n_int <- length(winner_index)

grp_labels <- unique(group)
n_party <- length(unique(group))
party_index <- as.numeric(sapply(group, function(x)which(c("five", "six", "sixi", "sixw", "nine") == x)))
party_size <- rowSums(presence)

# map male ids from elo data to male ids from reproductive success data
model_male_index_elo <- as.integer(sapply(demdata$widename, function(x)which(all_ids == x)))
names(model_male_index_elo) <- demdata$widename

n_model_males <- length(unique(demdata$id))
model_male_index_per_obs <- as.numeric(factor(demdata$id, levels = unique(demdata$id)))
names(model_male_index_per_obs) <- unique(demdata$id)[model_male_index_per_obs]

nonprime <- as.numeric(demdata$age_cat == "nonprime")

model_year_index_per_obs <- as.numeric(sapply(demdata$year, function(x)which(sort(unique(demdata$year)) == x)))
n_model_years <- length(unique(model_year_index_per_obs))
names(model_year_index_per_obs) <- as.character(demdata$year)

model_party_index_per_obs <- as.numeric(sapply(demdata$party, function(x)which(sort(unique(demdata$party)) == x)))
n_model_parties <- length(unique(model_party_index_per_obs))
names(model_party_index_per_obs) <- demdata$party

standat <- list(# data for Elo part of the model
                n_ind_elo = ncol(presence), 
                n_int = nrow(presence),
                n_party = 5,
                winner_index = winner_index,
                loser_index = loser_index,
                party_index = party_index,
                party_size = party_size,
                presence = presence,
                start_index = start_index, # index for migrants/group splits
                # actual model data
                n_model_obs = nrow(demdata), # n of observations in model
                model_fems = demdata$nfem, # response 
                rating_index = demdata$ratindex, 
                model_male_index_elo = model_male_index_elo,
                nonprime = nonprime,
                n_model_males = n_model_males,
                model_male_index_per_obs = model_male_index_per_obs,
                n_model_years = n_model_years,
                model_year_index_per_obs = model_year_index_per_obs,
                n_model_parties = n_model_parties,
                model_party_index_per_obs = model_party_index_per_obs
)

# add cumulative interaction numbers ----
if ("date" %in% colnames(intdata)) intdata$year <- as.numeric(substr(intdata$date, 1, 4))

irdata <- data.frame(id = names(standat$model_male_index_per_obs))
irdata$year <- as.numeric(names(standat$model_year_index_per_obs))
irdata$party <- names(standat$model_party_index_per_obs)
irdata$cumint <- 0
for (i in seq_len(nrow(irdata))) {
  m <- irdata$id[i]
  y <- irdata$year[i]
  p <- irdata$party[i]
  irdata$cumint[i] <- sum((intdata$winner == m | intdata$loser == m) & intdata$year <= y & intdata$party == p)
}
standat$cumint <- irdata$cumint

# indicator for prior-predictive checks
standat$do_model <- 1
standat$mig_mode <- 1

# descriptives for RS data ----
rs_descr <- list(all = list(n = nrow(demdata),
                            mean = sprintf(sprintf("%.2f", mean(demdata$nfem))),
                            median = median(demdata$nfem),
                            range = range(demdata$nfem)),
                 prime = list(n = sum(demdata$age_cat == "prime"),
                              mean = sprintf(sprintf("%.2f", mean(demdata$nfem[demdata$age_cat == "prime"]))),
                              median = median(demdata$nfem[demdata$age_cat == "prime"]),
                              range = range(demdata$nfem[demdata$age_cat == "prime"])),
                 nonprime = list(n = sum(demdata$age_cat == "nonprime"),
                                 mean = sprintf(sprintf("%.2f", mean(demdata$nfem[demdata$age_cat == "nonprime"]))),
                                 median = median(demdata$nfem[demdata$age_cat == "nonprime"]),
                                 range = range(demdata$nfem[demdata$age_cat == "nonprime"]))
)



# clean up workspace...
rm(model_party_index_per_obs, model_year_index_per_obs)
rm(all_ids, n_ind, n_int, n_model_males, n_party, nonprime, party_index, party_size)
rm(loser_index, winner_index, group, aux, grp_labels)
rm(model_male_index_elo, model_male_index_per_obs, n_model_parties, n_model_years)
rm(start_index, loser, winner, yy)
rm(i, m, p, y)
