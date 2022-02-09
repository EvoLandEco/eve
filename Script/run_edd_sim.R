args <- commandArgs(TRUE)

name <- args[1]
set <- as.numeric(args[2])
nrep <- as.numeric(args[3])

combo <- eve::edd_combo_maker(
  la = c(0.5, 0.8),
  mu = c(0, 0.2),
  beta_n = c(-0.01, 0),
  beta_phi = c(-0.01, 0),
  gamma_n = c(0.01, 0),
  gamma_phi = c(0.01, 0),
  age = c(20),
  model = "dsde2",
  metric = c("ed"),
  offset = "none"
)

out <- eve:::edd_sim_rep(combo = combo[[set]], nrep = nrep)

eve:::check_folder(name)
eve:::save_result(out, name, set)

