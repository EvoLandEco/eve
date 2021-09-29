devtools::install_github("rsetienne/DDD@tianjian_Rampal")
devtools::install_github("EvoLandEco/eve")
devtools::install_github("DavisVaughan/furrr")
install.packages("microbenchmark")

combo <- eve::edd_combo_maker(
  la = c(0.5, 0.8),
  mu = c(0, 0.1, 0.2),
  beta_n = c(-0.001, 0),
  beta_phi = c(-0.001, 0, 0.001),
  gamma_n = c(-0.001, 0.001),
  gamma_phi = c(-0.001, 0, 0.001),
  age = c(5, 10, 15),
  model = "dsde2",
  metric = c("ed"),
  offset = "none"
)


bm <- microbenchmark::microbenchmark(
  edd_go(combo = combo, nrep = 3),
  edd_go(
    combo = combo2,
    nrep = 3,
    strategy = future::multisession,
    workers = 32
  ),
  times = 5
)

write.csv2(summary(bm), "benchmark/bm.csv")
