#ifndef DSP_BAYES_SRC_PROPOSAL_FCNS_H
#define DSP_BAYES_SRC_PROPOSAL_FCNS_H


class ProposalFcns {

public:

    // sampling
    static double unif(double cond, double delta);
    static double abs_unif(double val, double delta);

    // density functions
    static double log_den_unif(double val, double cond, double delta);

};


#endif
