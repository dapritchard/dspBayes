# set variables ----------------------------------------------------------------

# Rcpp and R header files, and R math library location
r_incl_loc := /usr/share/R/include/
r_lib_loc := /usr/lib/R/lib
rcpp_incl_loc := /home/dpritch/R/x86_64-pc-linux-gnu-library/3.4/Rcpp/include
rcpp_lib_loc := /home/dpritch/R/x86_64-pc-linux-gnu-library/3.4/Rcpp/libs

# compiler directive
CC := $(CXX)

# preprocessor settings
CPPFLAGS := -D_FORTIFY_SOURCE=2 -I$(r_incl_loc) -I$(rcpp_incl_loc)

# compiler settings
CXXFLAGS := -Wall -g3 -fpic -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time

# linker settings
LDFLAGS :=




# primary makefile targets -----------------------------------------------------

.PHONY : all

# default rule
all :

clean :
	rm -f *.o




# specific build rules ---------------------------------------------------------


InputData.o : InputData.h read_data.h

test : InputData.o read_data.o

print_data_read.o : InputData.h

# print_data_read : LDFLAGS := -L$(r_lib_loc) -L$(rcpp_lib_loc) -Wl,-rpath=$(r_lib_loc) \
# 	-Wl,-rpath=$(rcpp_lib_loc)
print_data_read : LDFLAGS := -L$(rcpp_lib_loc)
print_data_read : LDLIBS := -lR -lRmath -l:Rcpp.so
print_data_read : InputData.o
