#include "CoefGen.h"
#include "UGen.h"
#include "UGenVar.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "WGen.h"
#include "XGen.h"
#include "XiGen.h"

#define U_MISS_CONTIN  0
#define U_MISS_CATEG   1

using Rcpp::as;




UGen::UGen(Rcpp::NumericMatrix& u_rcpp,
	   Rcpp::List& miss_info,
	   Rcpp::IntegerVector& miss_type,
	   Rcpp::IntegerVector& preg_map,
	   Rcpp::IntegerVector& sex_map) :
    m_vars(new UGenVar*[miss_info.size()]),
    m_n_vars(miss_info.size()) {

    // each iteration initialize one of the `UGenVar` subclasses and stores a
    // pointer to the data in `m_vars[i]`
    for (int i = 0; i < m_n_vars; ++i) {

	if (miss_type[i] == U_MISS_CATEG) {

	    Rcpp::List curr_var(miss_info[i]);

	    Rcpp::IntegerVector var_info     ( as<Rcpp::IntegerVector>(curr_var["var_info"])      );
	    Rcpp::NumericVector u_prior_probs( as<Rcpp::NumericVector>(curr_var["u_prior_probs"]) );
	    Rcpp::List var_block_list        ( as<Rcpp::List>(curr_var["var_block_list"])         );

	    m_vars[i] = new UGenVarCateg(u_rcpp,
					 var_info,
					 u_prior_probs,
					 var_block_list,
					 preg_map,
					 sex_map);
	}
	else {
	    // TODO: have to construct continuous version still
	}
    }
}




UGen::~UGen() {
    delete[] m_vars;
}




void UGen::sample(const WGen& W,
		  const XiGen& xi,
		  const CoefGen& coefs,
		  const XGen& X,
		  UProdBeta& ubeta,
		  UProdTau& utau) {

    for (int i = 0; i < m_n_vars; ++i) {
	m_vars[i]->sample(W, xi, coefs, X, ubeta, utau);
    }
}
