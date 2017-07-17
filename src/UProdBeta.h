#ifndef DSP_BAYES_SRC_U_PROD_BETA_H
#define DSP_BAYES_SRC_U_PROD_BETA_H


class UProdBeta {

public:

    double* m_vals;
    double* m_exp_vals;
    const int m_n_days;

    UProdBeta(int n) : m_n_days(n) {
	m_vals = new double[m_n_days];
	m_exp_vals = new double[m_n_days];
    }
    ~UProdBeta() {
	delete[] m_vals;
	delete[] m_exp_vals;
    }
    double* vals() { return m_vals; }
    double* exp_vals() { return m_exp_vals; }
    int n_days() { return m_n_days; }
    void update_exp(double* X);
};


void UProdBeta::update_exp(double* X) {
    for (int i = 0; i < m_n_days; i++) {
	m_exp_vals[i] = X[i] ?
	    std::exp(m_vals[i]) :
	    0;
    }
}


#endif
