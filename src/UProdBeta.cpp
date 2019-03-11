#include <cmath>
#include "UProdBeta.h"


UProdBeta::UProdBeta(int n_days) :
    // initialization list
    m_vals(new double[n_days]),
    m_exp_vals(new double[n_days]),
    m_n_days(n_days)
{
    // initialize `U * beta` values to 0 (i.e. `beta` values are all 0,
    // corresponding to no effect in the model)
    for (int i = 0; i < m_n_days; ++i) {
        m_vals[i] = 0.0;
        m_exp_vals[i] = 1.0;
    }
}




UProdBeta::~UProdBeta() {
    delete[] m_vals;
    delete[] m_exp_vals;
}




// change the values of the data pointed to by `ubeta_no_h` so that each element
// has the value of the corresponding element of `U_h * beta_h* added to it

void UProdBeta::add_uh_prod_beta_h(const double* U_h, double beta_h) {

    // each iteration updates the ijk-th element of `U * beta - U_h * beta_h` to
    // instead have the values given by `ubeta`.
    for (int r = 0; r < m_n_days; ++r) {

        // case: U_{ijkh} has a value of 1, so update the ijk-th element of
        // `U*beta - U_h*beta_h` to have the value of the ijk-th element of
        // `ubeta`.  If U_{ijkh} has a value of 0 then no update is needed.
        if (U_h[r]) {
            m_vals[r] += beta_h;
        }
    }
}




// update `U * beta` and `exp(U * beta)` based upon an updated value of
// `beta_h`
void UProdBeta::update(const double* U_h, double beta_h_new, double beta_h_curr) {

    const double beta_h_diff = beta_h_new - beta_h_curr;

    for (int i = 0; i < m_n_days; i++) {
        m_vals[i] += U_h[i] * beta_h_diff;
        m_exp_vals[i] = exp(m_vals[i]);
    }
}




// update `exp(U * beta)` based upon updated `U * beta`
void UProdBeta::update_exp() {
    for (int i = 0; i < m_n_days; i++) {
        m_exp_vals[i] = std::exp(m_vals[i]);
    }
}
