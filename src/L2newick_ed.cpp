#include <array>
#include <Rcpp.h>
#include "L2newick_ed.h"
#include "util.h"

// [[Rcpp::export]]
std::string l_to_newick_ed_cpp(const Rcpp::NumericMatrix& ltable_R,
                           const double t,
                           bool drop_extinct) {
  auto ltable_cpp = convert_to_ltable(ltable_R);
  auto newick_string = ltable_to_newick_ed(ltable_cpp, t, drop_extinct);
  return newick_string;
}