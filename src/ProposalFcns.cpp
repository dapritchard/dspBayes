#include "Rcpp.h"
#include "ProposalFcns.h"

#define LOG_ONE_HALF -0.6931471805599453




// sample from a `unif(cond - delta, cond + delta)` distribution.
//
// The expression below is due to the fact that the formula to sample a `unif(a,
// b)` distribution from `u ~ unif(0, 1)` is `a + (b - a) * u` so for `a = cond
// - delta` and `b = cond + delta` we obtain:
//
//     (cond - delta) + [(cond + delta) - (cond - delta)] * u
//
//         = (cond - delta) + (2 * delta * u)

double ProposalFcns::unif(double cond, double delta) {
    return (cond - delta) + (2.0 * delta * R::unif_rand());
}




// sample from an `abs( unif(cond - delta, cond + delta) )` distribution

double ProposalFcns::abs_unif(double cond, double delta) {

    double x = unif(cond, delta);

    // return the absolute value of `x`
    return (x < 0.0) ? -x : x;
}




// density function for a `unif(cond - delta, cond + delta)` distribution

double ProposalFcns::log_den_unif(double val, double cond, double delta) {

    // value given by `log(1 / (2 * delta))`
    return LOG_ONE_HALF - log(delta);
}
