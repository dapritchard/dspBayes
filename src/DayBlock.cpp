#include "Rcpp.h"
#include "DayBlock.h"

using Rcpp::IntegerVector;
using Rcpp::as;




// construct an array of `DayBlock`s based upon an `Rcpp::List` that provides
// the specifications for each block.  In more detail, an array of `DayBlock`s
// is created with length given by the length of `block_list`.  Furthermore, it
// is assumed that each element in `block_list` stores integer values that can
// be accessed using names `beg_idx` and `n_days` and which are used to
// initialize the struct's member values of the same name.  The return value is
// a pointer to the beginning of the array.

DayBlock* DayBlock::list_to_arr(Rcpp::List& block_list) {

    DayBlock* block_arr = new DayBlock[block_list.size()];

    // each iteration constructs a new struct based upon the information
    // provided by the t-th element of `block_list`
    for (int t = 0; t < block_list.size(); ++t) {

	block_arr[t] = DayBlock(as<IntegerVector>(block_list[t])["beg_idx"],
				as<IntegerVector>(block_list[t])["n_days"]);
    }

    return block_arr;
}




// similar to `DayBlock::list_to_arr`, but now each element in `block_list`
// stores integer values that can be accessed using names `beg_idx`, `n_days`,
// and `subj_idx`, and which are used to initialize the struct's member values
// of the same name.  The return value is a pointer to the beginning of the
// array.

PregCyc* PregCyc::list_to_arr(Rcpp::List& block_list) {

    PregCyc* block_arr = new PregCyc[block_list.size()];

    // each iteration constructs a new struct based upon the information
    // provided by the t-th element of `block_list`
    for (int t = 0; t < block_list.size(); ++t) {

	block_arr[t] = PregCyc(as<IntegerVector>(block_list[t])["beg_idx"],
			       as<IntegerVector>(block_list[t])["n_days"],
			       as<IntegerVector>(block_list[t])["subj_idx"]);
    }

    return block_arr;
}




// similar to `DayBlock::list_to_arr` and `PregCyc::list_to_arr`, but now each
// element in `block_list` stores integer values that can be accessed using
// names `beg_idx`, `n_days`, `subj_idx`, and `n_miss`, and which are used to
// initialize the struct's member values of the same name.  The return value is
// a pointer to the beginning of the array.

MissCyc* MissCyc::list_to_arr(Rcpp::List& block_list) {

    MissCyc* block_arr = new MissCyc[block_list.size()];

    // each iteration constructs a new struct based upon the information
    // provided by the t-th element of `block_list`
    for (int t = 0; t < block_list.size(); ++t) {

	Rcpp::IntegerVector block_list_t = as<Rcpp::IntegerVector>(block_list[t]);
	block_arr[t] = MissCyc(block_list_t["beg_idx"],
			       block_list_t["n_days"],
			       block_list_t["subj_idx"],
			       block_list_t["n_miss"],
			       block_list_t["preg"]);
    }

    return block_arr;
}
