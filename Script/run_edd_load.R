args <- commandArgs(TRUE)

name <- args[1]

sim_data <- eve::edd_load(path = name, strategy = "multicore", workers = 8, verbose = TRUE)

print("Saving loaded data")
save(sim_data, file = file.path(name, "sim_data.RData"))

