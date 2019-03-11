#ifndef DSP_BAYES_FW_PRIORS_H_
#define DSP_BAYES_FW_PRIORS_H_


class MDay {
public:

    double m_peak_idx;

    MDay() {}
    double val() { return m_peak_idx; }

    // void sample();
};




class Mu {
public:

    double m_mu;

    double val() { return m_mu; }
};



class Nu {
public:

    double m_nu;

    double val() { return m_nu; }
};




class Delta {
public:

    double m_delta;

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
