devtools::install_github("rsetienne/DDD@tianjian_Rampal")
devtools::install_github("EvoLandEco/eve")

combo <- eve::edd_combo_maker(
  la = c(0.5, 0.8),
  mu = c(0.1, 0.2),
  beta_n = c(-0.001),
  beta_phi = c(-0.001, 0.001),
  gamma_n = c(-0.001, 0.001),
  gamma_phi = c(-0.001, 0.001),
  age = c(5),
  model = "dsde2",
  metric = c("ed"),
  offset = "none"
)


bm <- microbenchmark::microbenchmark(
  eve::edd_go(combo = combo, nrep = 3, name = "test"),
  eve::edd_go(
    combo = combo,
    nrep = 3,
    name = "test",
    strategy = future::multisession,
    workers = 16
  ),
  eve::edd_go(
    combo = combo,
    nrep = 3,
    name = "test",
    strategy = future::multisession,
    workers = 32
  ),
  times = 5
)

write.csv2(summary(bm), "benchmark/bm.csv")
