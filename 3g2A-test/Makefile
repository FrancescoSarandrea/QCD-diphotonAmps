LOCAL?=$(HOME)/local
GIT?=$(HOME)/git
NJET_LOCAL?=$(LOCAL)/njet

CXX=LD_LIBRARY_PATH=$(NJET_LOCAL)/lib g++
CXXFLAGS=-O0 -std=c++17 -Wall -Wextra
CPPFLAGS=-I$(NJET_LOCAL)/include -I$(GIT)/njet-develop
LDFLAGS=-L$(NJET_LOCAL)/lib
DEFS= -DUSE_DD -DUSE_QD
DEFS+= -DDEBUG
LIBS= -lnjet2 -lnjet2tools -lqd

SRCS = $(wildcard *.cpp)
PRGS = $(basename $(SRCS))

all: $(PRGS)

$(PRGS): %: %.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(CPPFLAGS) $(DEFS)

clean:
	rm -f *.o $(PRGS)
