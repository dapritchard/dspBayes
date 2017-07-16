

DayBlock* DayBlock::list_to_arr(Rcpp::List block_list) {

    DayBlock* block_arr = new DayBlock[block_list.size()];

    // create `DayBlock` structs that provide the indices in the day-specific
    // data that correspond to the t-th individual
    for (int t = 0; t < m_n_subj; ++t) {
	block_arr[t] = DayBlock(block_list[t]["beg_idx"], block_list[t]["n_days"]);
    }

    return block_arr;
}


PregCyc* PregCyc::list_to_arr(Rcpp::List block_list) {

    PregCyc* block_arr = new PregCyc[block_list.size()];

    // create `PregCyc` structs that provide the indices in the day-specific
    // data that correspond to the t-th individual
    for (int t = 0; t < m_n_subj; ++t) {
	block_arr[t] = PregCyc(block_list[t]["beg_idx"],
			       block_list[t]["n_days"],
			       block_list[t]["subj_idx"]);
    }

    return block_arr;
}
