args <- commandArgs(TRUE)

name <- args[1]
set <- as.numeric(args[2])
nrep <- as.numeric(args[3])

combo1 <- eve::edd_combo_maker(
  la = c(0.5, 0.8),
  mu = c(0, 0.1, 0.2),
  beta_n = c(-0.02, -0.01, 0),
  beta_phi = c(-0.05, -0.02, 0, 0.0005, 0.001),
  age = c(6),
  model = "dsce2",
  metric = c("pd"),
  offset = c("none", "simtime", "spcount")
)

combo2 <- eve::edd_combo_maker(
  la = c(0.5, 0.8),
  mu = c(0, 0.1, 0.2),
  beta_n = c(-0.02, -0.01, 0),
  beta_phi = c(-0.05, -0.02, 0, 0.0005, 0.001),
  age = c(6),
  model = "dsce2",
  metric = c("ed"),
  offset = c("none")
)

combo <- c(combo1, combo2)

names(combo) <- as.character(seq(1, length(combo)))

out <-
  eve:::edd_sim_rep(
    combo = combo[[set]],
    history = TRUE,
    verbose = FALSE,
    nrep = nrep
  )

print(paste0("Starting parameter set ", set))

eve:::check_folder(name)
eve:::save_result(out, name, set)
