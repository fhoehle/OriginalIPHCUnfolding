### makfile for unfolding macros
SOURCES=classy_unfold1d.cpp classy_unfold.cpp datamig.cpp findbinning.cpp lincheck.cpp migrationmatrix.cpp plotsforpseudo.cpp pseudoexp.cpp reweighted_pseudoexp.cpp seleffmatrix.cpp
#OBJECTS=$(SOURCES:.cpp=.o)
MYTARGETS=$(SOURCES:.cpp=)
.SUFFIXES: .exec .cpp

EXECUTABLES=$(SOURCES:.cpp=.exec)

CFLAGS = $(shell root-config --cflags)
LIBS   = $(shell root-config --libs)

#all: $(EXECUTABLES)

all: $(MYTARGETS)

#%.exec : %.cpp
#	g++ $< -g -O2 -fPIC -o $@ $(CFLAGS) $(LIBS)

#.cpp.exec:
#       g++ $< -g -O2 -fPIC -o $@ $(CFLAGS) $(LIBS)
$(MYTARGETS): 
	g++ $@.cpp -g -O2 -fPIC -o $@.exec $(CFLAGS) $(LIBS)

clean :
	-rm $(EXECUTABLES)
