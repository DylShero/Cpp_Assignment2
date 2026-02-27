CXX 		:= g++
CXXFLAGS 	:= --std=c++23 -Wextra -Wall  
LDFLAGS 	:= 


SRCS 		:= assignment2.cc

EXES 		:= ${SRCS:.cc=}

all: ${EXES} 

${EXES}: %: %.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)


.PHONY: clean all

clean:
	$(RM) $(EXES) 