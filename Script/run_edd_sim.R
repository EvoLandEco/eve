args <- commandArgs(TRUE)

name <- args[1]
set <- as.numeric(args[2])
nrep <- as.numeric(args[3])

combo1 <- eve::edd_combo_maker(
  la = c(0.5),
  mu = c(0.1),
  beta_n = c(-0.01),
  beta_phi = c(0.02),
  gamma_n = c(0.01),
  gamma_phi = c(-0.05),
  age = c(15),
  model = "dsde2",
  metric = c("ed"),
  offset = c("none")
)

combo2 <- eve::edd_combo_maker(
  la = c(0.5),
  mu = c(0.1),
  beta_n = c(0, -0.01),
  beta_phi = c(0.01),
  gamma_n = c(0, 0.01),
  gamma_phi = c(-0.005),
  age = c(15),
  model = "dsde2",
  metric = c("pd"),
  offset = c("none")
)

combo <- c(combo1, combo2)
names(combo) <- as.character(seq(1, length(combo)))

out <- eve:::edd_sim_rep(combo = combo[[set]], nrep = nrep)

eve:::check_folder(name)
eve:::save_result(out, name, set)
