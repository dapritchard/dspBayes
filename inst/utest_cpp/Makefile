# specify C++ flags.  The implicit rule used to link together object files
# utilizes $(CC) which by default has value cc, so this needs to be changed to
# invoke the C++ toolchain.
rcpp_loc := /home/dpritch/R/x86_64-pc-linux-gnu-library/3.4/Rcpp/include
r_loc := /usr/share/R/include/
CXXFLAGS := -Wall -g3 $(CXXFLAGS) -I$(rcpp_loc) -I$(r_loc)
CC := $(CXX)


# primary makefile targets -----------------------------------------------------

.PHONY : all

# default rule
all :




InputData.o : InputData.h read_data.h
