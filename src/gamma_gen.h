



class GammaGen {

 public:

    // current value of gamma_h
    double val;

    // gamma_h hyperparameters
    double hyp_a, hyp_b, hyp_p, bnd_l, bnd_u;

    // whether h-th variable categorical or continuous
    bool is_categ;

    virtual double samp_gam(const vector<double>& U_prod_beta,
			    const vector<double>& U_h,
			    const vector<double>& W);
};




class GammaGenCateg : public GammaGen {

public:

    // whether h-th variable is truncated or not
    bool is_trunc;

    // C(a_h, b_h)
    double c_of_ab;


};




class GammaGenContAux : public GammaGen {

public:

    // // whether M == 1 (i.e no effect in the model)
    // bool M_is_one;


};




class GammaGenContMetr : public GammaGen {

public:

    // Metropolis-Hastings tuning parameter.  Specifies the probability of
    // selecting 1 as the proposal value.
    double mh_pi;

};
