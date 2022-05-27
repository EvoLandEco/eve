#' @name calculate_CI
#' @title Calculating confidence interval for LTT plot (reversed time scale)
#' @author Thijs Janzen; Tianjian Qin
calculate_CI <- Rcpp::cppFunction('
NumericMatrix get_stats_cpp(const List& brts_list,
                            const NumericVector& tt,
                            const float alpha) {
  NumericMatrix output(tt.size(), 5);

  size_t min_index = brts_list.size() * alpha/2;
  size_t max_index = brts_list.size() * (1 - alpha/2);

  for (int i = (tt.size() - 1); i >= 0; --i) {
    float focal_t = tt[i];
    std::vector<size_t> found_lin(brts_list.size());
    float mean = 0.f;
    NumericVector focal_brts;
    for (int j = 0; j < brts_list.size(); ++j) {
      focal_brts = brts_list[j]; // temp
      size_t cnt = focal_brts.size() - 1;
      while(focal_brts[cnt] > focal_t) {
        cnt--;
      }
      found_lin[j] = cnt; // two at first timepoint, crown age.
      mean += found_lin[j];
    }
    mean *= 1.0 / brts_list.size();

    std::sort(found_lin.begin(), found_lin.end());
    float median = found_lin[ found_lin.size() / 2 + 1 ];

    if (found_lin.size() % 2 == 0) {
      median = 0.5 * (found_lin[ found_lin.size() / 2] +
                      found_lin[ found_lin.size() / 2 + 1]);
    }

    float minalpha = found_lin.at( min_index );
    float maxalpha = found_lin.at( max_index );

    NumericVector add(5);
    add(0) = focal_t;
    add(1) = median;
    add(2) = minalpha;
    add(3) = maxalpha;
    add(4) = mean;

    output(i, _) = add;
  }
  return output;
}')



calculate_aldous_beta <- function(phy = NULL, method = "treestats") {
  if (!ape::is.ultrametric.phylo(phy)) {
    stop("Ultrametric tree required, considering using tree with only extant taxa")
  }

  if (method == "treestats") {
    beta <- treestats::beta_statistic(phy)
  } else {
    stop("No such method")
  }
}
