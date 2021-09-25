ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
CXXFLAGS  += $(ROOTCFLAGS) -I$(P56xxLIBS) -Wall -O3
LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS) -Wl,-rpath,$(P56xxLIBS) -L$(P56xxLIBS) -lP56xx
GXX	   = g++ $(CXXFLAGS)

SRCS = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)
EXES = $(SRCS:%.cpp=%)
dep = $(OBJ:.o=.d)  # one dependency file for each source


all: $(EXES)

$(EXES):  $(SRCS)
	$(GXX) $(CXXFLAGS) -MMD -c $@.cpp
	$(GXX) $@.cpp -o$@ $(LDFLAGS)

clean:
	rm -f *.o *.so *.d *~ $(EXES) $(dep)


