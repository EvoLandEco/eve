args <- commandArgs(TRUE)

name <- args[1]
set <- as.numeric(args[2])
nrep <- as.numeric(args[3])

combo1 <- eve::edd_combo_maker(
  la = c(0.5),
  mu = c(0.1),
  beta_n = c(-0.02, 0),
  beta_phi = c(-0.02),
  gamma_n = c(0.02, 0),
  gamma_phi = c(0.02),
  age = c(20),
  model = "dsde2",
  metric = c("ed"),
  offset = c("none")
)

combo2 <- eve::edd_combo_maker(
  la = c(0.5),
  mu = c(0.1),
  beta_n = c(-0.02, 0),
  beta_phi = c(-0.02),
  gamma_n = c(0.02, 0),
  gamma_phi = c(0.02),
  age = c(20),
  model = "dsde2",
  metric = c("pd"),
  offset = c("none", "simtime")
)

combo3 <- eve::edd_combo_maker(
  la = c(0.5),
  mu = c(0.1),
  beta_n = c(0),
  beta_phi = c(-0.02),
  age = c(20),
  model = "dsce2",
  metric = c("ed"),
  offset = c("none")
)

combo4 <- eve::edd_combo_maker(
  la = c(0.5),
  mu = c(0.1),
  beta_n = c(0),
  beta_phi = c(-0.02),
  age = c(20),
  model = "dsce2",
  metric = c("pd"),
  offset = c("none", "simtime")
)

combo <- c(combo1, combo2, combo3, combo4)
names(combo) <- as.character(seq(1, length(combo)))

out <- eve:::edd_sim_rep(combo = combo[[set]], nrep = nrep)

eve:::check_folder(name)
eve:::save_result(out, name, set)
