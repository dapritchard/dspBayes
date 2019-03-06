#ifndef DSP_BAYES_FW_PRIORS_H_
#define DSP_BAYES_FW_PRIORS_H_


class MDay {
public:

    double m_peak_idx;

    MDay() {}
    double val() { return m_peak_idx; }
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

    MDay& mday;
    Mu& mu;
    Nu& nu;
    Delta& delta;

    FWPriors(MDay mday, Mu mu, Nu nu, Delta delta) :
	mday{mday}, mu{mu}, nu{nu}, delta{delta}
    {}
};


#endif
