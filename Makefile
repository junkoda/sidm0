#
# Zero dimensional sidm
#

EXEC      = sidm0

all: $(EXEC)

CXX       := c++
CXXFLAGS  := -O2


OBJS1     = sidm0.o scatter.o
LIBS1     = -lboost_program_options -lgsl -lgslcblas
DIRS      = $(GSL_DIR) $(BOOST_DIR)

INCLDIRS  = $(foreach dir, $(DIRS), -I$(dir)/include)
LIBDIRS   = $(foreach dir, $(DIRS), -L$(dir)/lib)
CXXFLAGS += $(INCLDIRS)

sidm0: $(OBJS1)
	$(CXX) $(OBJS1) $(LIBS1) -o $@

scatter.o: scatter.cpp scatter.h particle.h
sidm0.o: sidm0.cpp particle.h scatter.h

.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS1)
