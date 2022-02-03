args <- commandArgs(TRUE)

name <- args[1]
set <- args[2]
la <- args[3]
mu <- as.numeric(args[4])
beta_n <- as.numeric(args[5])
beta_phi <- as.numeric(args[6])
gamma_n <- as.numeric(args[7])
gamma_phi <- as.numeric(args[8])
age <- as.numeric(args[9])
metric <- as.numeric(args[10])
model <- as.character(args[11])
offset <- as.character(args[12])

pars <- c(la, mu, beta_n, beta_phi, gamma_n, gamma_phi)

combo <- c(pars, age, metric, model, offset)

out <- eve:::edd_sim_rep(combo = combo, nrep = 3)

eve:::save_result(out, name, set)
