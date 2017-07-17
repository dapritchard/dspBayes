#include "Rcpp.h"
#include "ProposalFcns.h"


double ProposalFcns::abs_unif(double val, double delta) {

    // sample from a `unif(val - delta, val + delta)` distribution.  The
    // expression below is due to the fact that the formula to sample a `unif(a,
    // b)` distribution from `u ~ unif(0, 1)` is `a + (b - a) * u`.
    double x = val - delta  + 2 * delta * R::unif_rand();

    // return the absolute value of `x`
    return (x < 0) ? -x : x;
}
