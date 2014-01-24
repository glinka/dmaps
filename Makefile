SRCS=main.cc util_fns.cc #my_allocator.cc
OBJECTS=$(SRCS:.cc=.o)

CXX = g++

CXXFLAGS = -g -Wall -std=c++0x #-O3

all: dmaps

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

dmaps: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^ > ./.depend

clean:
	$(RM) *.o 

include .depend
