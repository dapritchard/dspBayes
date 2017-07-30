#ifndef DSP_BAYES_SRC_U_PROD_BETA_H
#define DSP_BAYES_SRC_U_PROD_BETA_H


class UProdBeta {

public:

    double* m_vals;
    double* m_exp_vals;
    const int m_n_days;

    UProdBeta(int n_days);
    ~UProdBeta();

    void add_uh_prod_beta_h(const double* U_h, double beta_h);

    double* vals() { return m_vals; }
    const double* vals() const { return m_vals; }
    double* exp_vals() { return m_exp_vals; }
    const double* exp_vals() const {return m_exp_vals; }
    int n_days() { return m_n_days; }

    void update_exp(int* X);
};


#endif
