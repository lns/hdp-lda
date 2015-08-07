VERSION = 1

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O2 

all: exp

main: main.cpp *.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

exp: main
	./main -data ap/ap.dat -ndoc 2246 -nword 10473
	R -q -f print_topic.R

clean:
	$(RM) main

