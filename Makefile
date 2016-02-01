SRCS=main.cc dmaps_util_fns.cc
OBJECTS=$(SRCS:.cc=.o)

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x -leigensolvers -lutil_fns -L/home/alexander/local/lib -I/home/alexander/workspace/newton_gmres -I/home/alexander/workspace/util_fns -I/home/alexander/local/eigen -O3

all: dmaps

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

dmaps: $(OBJECTS)
	$(CXX) -o $@ $^ $(CXXFLAGS)

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^ > ./.depend

clean:
	$(RM) *.o 

include .depend
