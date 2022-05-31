args <- commandArgs(TRUE)

name <- args[1]

edd_data <- eve::edd_load(name = name, strategy = "multicore", workers = 8)

eve::test_plot_las(edd_data)
eve::test_plot_eds(edd_data)
eve::test_plot_ltt(edd_data)
