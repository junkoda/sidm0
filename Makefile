#
# Zero dimensional sidm
#

EXEC      = sidm0

all: $(EXEC)

CXXFLAGS  = 


OBJS1     = sidm0.o
LIBS1     = -lboost_program_options -lgsl -lgslcblas
DIRS      = $(GSL_DIR) $(BOOST_DIR)

INCLDIRS  = $(foreach dir, $(DIRS), -I$(dir)/include)
LIBDIRS   = $(foreach dir, $(DIRS), -L$(dir)/lib)
CXXFLAGS := $(INCLDIRS)

sidm0: $(OBJS1)
	$(CXX) $(OBJS1) $(LIBS1) -o $@

.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS1)
