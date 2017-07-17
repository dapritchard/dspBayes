#include "Rcpp.h"
#include "DayBlock.h"


// construct an array of `DayBlock`s based upon an `Rcpp::List` that provides
// the specifications for each block.  In more detail, an array of `DayBlock`s
// is created with length given by the length of `block_list`.  Furthermore, it
// is assumed that each element in `block_list` stores integer values that can
// be accessed using names `beg_idx` and `n_days` and which are used to
// initialize the struct's member values of the same name.  The return value is
// pointer to the beginning of the array.

DayBlock* DayBlock::list_to_arr(Rcpp::List& block_list) {

    DayBlock* block_arr = new DayBlock[block_list.size()];

    // each iteration constructs a new struct based upon the information
    // provided by the t-th element of `block_list`
    for (int t = 0; t < block_list.size(); ++t) {
	block_arr[t] = DayBlock(block_list[t]["beg_idx"],
				block_list[t]["n_days"]);
    }

    return block_arr;
}




// construct an array of `PregCyc`s based upon an `Rcpp::List` that provides the
// specifications for each block.  In more detail, an array of `PregCyc`s is
// created with length given by the length of `block_list`.  Furthermore, it is
// assumed that each element in `block_list` stores integer values that can be
// accessed using names `beg_idx`, `n_days`, and `subj_idx`, and which are used
// to initialize the struct's member values of the same name.  The return value
// is pointer to the beginning of the array.

PregCyc* PregCyc::list_to_arr(Rcpp::List& block_list) {

    PregCyc* block_arr = new PregCyc[block_list.size()];

    // each iteration constructs a new struct based upon the information
    // provided by the t-th element of `block_list`
    for (int t = 0; t < block_list.size(); ++t) {
	block_arr[t] = PregCyc(block_list[t]["beg_idx"],
			       block_list[t]["n_days"],
			       block_list[t]["subj_idx"]);
    }

    return block_arr;
}
