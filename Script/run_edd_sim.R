args <- commandArgs(TRUE)

name <- args[1]
set <- as.numeric(args[2])
nrep <- as.numeric(args[3])

combo <- eve::edd_combo_maker(
  la = c(0.5),
  mu = c(0, 0.1),
  beta_n = c(-0.001, 0),
  beta_phi = c(-0.001, 0.001),
  gamma_n = c(-0.001, 0.001),
  gamma_phi = c(-0.001, 0.001),
  age = c(5),
  model = "dsde2",
  metric = c("ed"),
  offset = "none"
)

out <- eve:::edd_sim_rep(combo = combo[[set]], nrep = nrep)

eve:::save_result(out, name, set)

