# ````````````````````````````````````````````````````` #
#  example R package build commands with Rcpp Makevars  #
# ..................................................... #

# g++ -I/usr/share/R/include -DNDEBUG `/usr/lib/R/bin/Rscript -e "Rcpp:::CxxFlags()"`            \
#     -I"/home/dpritch/R/x86_64-pc-linux-gnu-library/3.4/Rcpp/include"    -fpic  -g -O2          \
#     -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2  \
#     -g  -c dsp_sampler.cpp -o dsp_sampler.o

# g++ -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o dspBayes.so CoefGen.o      \
#     DayBlock.o GammaCateg.o GammaGen.o PhiGen.o ProposalFcns.o RcppExports.o UProdBeta.o WGen.o  \
#     XiGen.o collect_data.o dsp_sampler.o -llapack -lblas -lgfortran -lm -lquadmath               \
#     -L/usr/lib/R/lib -lR




# `````````````````````````````````` #
#  object file compilation settings  #
# .................................. #

# -fpic
#     Generate position-independent code (PIC) suitable for use in a shared library, if supported
#     for the target machine.  Such code accesses all constant addresses through a global offset
#     table (GOT).  The dynamic loader resolves the GOT entries when the program starts (the
#     dynamic loader is not part of GCC; it is part of the operating system).  If the GOT size for
#     the linked executable exceeds a machine-specific maximum size, you get an error message from
#     the linker indicating that -fpic does not work; in that case, recompile with -fPIC instead.
#     (These maximums are 8k on the SPARC and 32k on the m68k and RS/6000.  The x86 has no such
#     limit.)
#
#     Position-independent code requires special support, and therefore works only on certain
#     machines.  For the x86, GCC supports PIC for System V but not for the Sun 386i.  Code
#     generated for the IBM RS/6000 is always position-independent.
#
#     When this flag is set, the macros "__pic__" and "__PIC__" are defined to 1.


# -fstack-protector
#     Emit extra code to check for buffer overflows, such as stack smashing attacks.  This is done
#     by adding a guard variable to functions with vulnerable objects.  This includes functions
#     that call "alloca", and functions with buffers larger than 8 bytes.  The guards are
#     initialized when a function is entered and then checked when the function exits.  If a guard
#     check fails, an error message is printed and the program exits.
#
# -fstack-protector-strong
#      Like -fstack-protector but includes additional functions to be protected --- those that have
#      local array definitions, or have references to local frame addresses.


# -Wdate-time
#     Warn when macros "__TIME__", "__DATE__" or "__TIMESTAMP__" are encountered as they might
#     prevent bit-wise-identical reproducible compilations.


# -Wformat
# -Wformat=n
#     Check calls to "printf" and "scanf", etc., to make sure that the arguments supplied have
#     types appropriate to the format string specified, and that the conversions specified in the
#     format string make sense.  This includes standard functions, and others specified by format
#     attributes, in the "printf", "scanf", "strftime" and "strfmon" (an X/Open extension, not in
#     the C standard) families (or other target-specific families).  Which functions are checked
#     without format attributes having been specified depends on the standard version selected,
#     and such checks of functions without the attribute specified are disabled by -ffreestanding
#     or -fno-builtin.
#
#     The formats are checked against the format features supported by GNU libc version 2.2.
#     These include all ISO C90 and C99 features, as well as features from the Single Unix
#     Specification and some BSD and GNU extensions.  Other library implementations may not
#     support all these features; GCC does not support warning about features that go beyond a
#     particular library's limitations.  However, if -Wpedantic is used with -Wformat, warnings
#     are given about format features not in the selected standard version (but not for "strfmon"
#     formats, since those are not in any version of the C standard).
#
#     -Wformat=1
#     -Wformat
#         Option -Wformat is equivalent to -Wformat=1, and -Wno-format is equivalent to
#         -Wformat=0.  Since -Wformat also checks for null format arguments for several functions,
#         -Wformat also implies -Wnonnull.  Some aspects of this level of format checking can be
#         disabled by the options: -Wno-format-contains-nul, -Wno-format-extra-args, and
#         -Wno-format-zero-length.  -Wformat is enabled by -Wall.


# -Werror=
#     Make the specified warning into an error.  The specifier for a warning is appended; for
#     example -Werror=switch turns the warnings controlled by -Wswitch into errors.  This switch
#     takes a negative form, to be used to negate -Werror for specific warnings; for example
#     -Wno-error=switch makes -Wswitch warnings not be errors, even when -Werror is in effect.
#
#     The warning message for each controllable warning includes the option that controls the
#     warning.  That option can then be used with -Werror= and -Wno-error= as described above.
#     (Printing of the option in the warning message can be disabled using the
#     -fno-diagnostics-show-option flag.)
#
#     Note that specifying -Werror=foo automatically implies -Wfoo.  However, -Wno-error=foo does
#     not imply anything.


# -Wformat-security
#     If -Wformat is specified, also warn about uses of format functions that represent
#     possible security problems.  At present, this warns about calls to "printf" and "scanf"
#     functions where the format string is not a string literal and there are no format
#     arguments, as in "printf (foo);".  This may be a security hole if the format string came
#     from untrusted input and contains %n.  (This is currently a subset of what
#     -Wformat-nonliteral warns about, but in future warnings may be added to
#     -Wformat-security that are not included in -Wformat-nonliteral.)


# Note: with FORTIFY_SOURCE whenever possible GCC tries to uses buffer-length aware replacements
# for functions like strcpy, memcpy, memset, etc.




# ````````````````````````````` #
#  shared object link settings  #
# ............................. #

# -Wl,option
#     Pass option as an option to the linker.  If option contains commas, it is split into
#     multiple options at the commas.  You can use this syntax to pass an argument to the option.
#     For example, -Wl,-Map,output.map passes -Map output.map to the linker.  When using the GNU
#     linker, you can also get the same effect with -Wl,-Map=output.map.

# -Bsymbolic-functions
#     When creating a shared library, bind references to global function symbols to the definition
#     within the shared library, if any.  This option can also be used with the --export-dynamic
#     option, when creating a position independent executable, to bind references to global
#     function symbols to the definition within the executable.  This option is only meaningful on
#     ELF platforms which support shared libraries and position independent executables.

# -z keyword
#     The recognized keywords are:
#
#     ...
#
#     relro
#         Create an ELF "PT_GNU_RELRO" segment header in the object.




# ```````````````````````` #
#  begin Makefile section  #
# ........................ #

# set variables ----------------------------------------------------------------

# cppUnit, R, and Rcpp header file locations
cppunit_incl_loc := /usr/local/include
r_incl_loc := /usr/share/R/include/
rcpp_incl_loc := /home/dpritch/R/x86_64-pc-linux-gnu-library/3.4/Rcpp/include

# cppUnit and R library locations
cppunit_lib_loc := /usr/local/lib
r_lib_loc := /usr/lib/R/lib

# compiler directive
CC := $(CXX)

# preprocessor settings.  Note: -DNDEBUG should be included for the final
# build.
CPPFLAGS := -D_FORTIFY_SOURCE=2 -I$(r_incl_loc) -I$(rcpp_incl_loc)

# compiler settings.  Note: -O2 should be included for the final build
CXXFLAGS := -Wall -g3 -std=c++11 -fpic -fstack-protector-strong -Wformat  \
            -Werror=format-security -Wdate-time

# linker settings
LDFLAGS := -shared -Wl,-Bsymbolic-functions -Wl,-z,relro -L$(r_lib_loc)

# libraries to link against.  Note: the default R build links against liblapack,
# libblas, libgfortran, libm, libquadmath, and libR.  The `addprefix` directive
# is used in case we want to add more libraries later.
LDLIBS := $(addprefix -l, R)




# primary makefile targets -----------------------------------------------------

utests := $(addsuffix .o, $(basename $(wildcard UTest*.cpp)))
targets := $(filter-out $(utests), $(addsuffix .o, $(basename $(wildcard *.cpp))))

.PHONY : all clean clobber print-%

# default rule
all : $(targets) $(utests) dspBayes.so

clean :
	rm -f *.o

clobber :
	rm -f *.o dspBayes.so

# for debugging the makefile.  Print out a variable name `varname' by running
# `make print-varname'
print-%:
	@echo '$*=$($*)'




# package build rules ----------------------------------------------------------

# -l library
#     Search the library named library when linking.  (The second alternative with the library as
#     a separate argument is only for POSIX compliance and is not recommended.)

#     It makes a difference where in the command you write this option; the linker searches and
#     processes libraries and object files in the order they are specified.  Thus, foo.o -lz bar.o
#     searches library z after file foo.o but before bar.o.  If bar.o refers to functions in z,
#     those functions may not be loaded.

dspBayes.so : override LDFLAGS += -L$(cppunit_lib_loc) -Wl,-rpath=$(cppunit_lib_loc)
dspBayes.so : override LDLIBS += -lcppunit
dspBayes.so : $(targets) $(utests)
	$(CC) $(targets) $(utests) $(LDFLAGS) $(LDLIBS) -o dspBayes.so

CoefGen.o : CoefGen.h

DayBlock.o : DayBlock.h

GammaCateg.o : GammaGen.h global_vars.h

GammaGen.o : GammaGen.h

PhiGen.o : PhiGen.h ProposalFcns.h XiGen.h

ProposalFcns.o : ProposalFcns.h

RcppExports.cpp : dsp_sampler.cpp UTestDriver.cpp
	Rscript -e 'Rcpp::compileAttributes("..")'
RcppExports.o : RcppExports.cpp

XiGen.o : XiGen.h PhiGen.h DayBlock.h UProdBeta.h

WGen.o : WGen.h XiGen.h DayBlock.h UProdBeta.h





# testing build rules ----------------------------------------------------------

utests : override CPPFLAGS += $(cpp_incl_loc)

UTestDriver.o : UTestPhiGen.h

UTestFactory.o : UTestFactory.h XiGen.h WGen.h PhiGen.h UProdBeta.h

UTestPhiGen.o : PhiGen.h UTestPhiGen.h XiGen.h

UTestXiGen.o : UTestXiGen.h XiGen.h WGen.h PhiGen.h UProdBeta.h

UTestWGen.o : UTestWGen.h WGen.h XiGen.h UProdBeta.h UTestFactory.h
