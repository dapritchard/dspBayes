#include <cmath>
#include "UProdBeta.h"


UProdBeta::UProdBeta(int n) :
    // initialization list
    m_vals(new double[n]),
    m_exp_vals(new double[n]),
    m_n_days(n)
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



void UProdBeta::update_exp(int* X) {
    for (int i = 0; i < m_n_days; i++) {
	m_exp_vals[i] = X[i] ?
	    std::exp(m_vals[i]) :
	    0;
    }
}
