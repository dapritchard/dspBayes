#ifndef DSP_BAYES_FW_PRIORS_H_
#define DSP_BAYES_FW_PRIORS_H_


class MDay {
public:

    double m_peak_idx;

    MDay() : m_peak_idx {3} {}
    double val() { return m_peak_idx; }

    // void sample();
};




class Mu {
public:

    double m_mu;

    Mu() : m_mu {0.44} {}
    double val() { return m_mu; }
};




class Nu {
public:

    double m_nu;

    Nu() : m_nu {1.0} {}
    double val() { return m_nu; }
};




class Delta {
public:

    double m_delta;

    Delta() : m_delta {0.5} {}
    double val() { return m_delta; }
};




// FIXME: need to call destructors for these classes?
class FWPriors {
public:

    MDay mday;
    Mu mu;
    Nu nu;
    Delta delta;

    FWPriors(Rcpp::List fw_prior_specs) :
        mday  {MDay()},
        mu    {Mu()},
        nu    {Nu()},
        delta {Delta()}
    {}
};


#endif
