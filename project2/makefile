CXX = c++
# Flags for optimized execs
CXXFLAGS = 
# CXXFLAGS = 
# Flags for debugging
#CXXFLAGS = -ftemplate-depth-30 -g -DBZ_DEBUG

LIBS = -lm -larmadillo

EXEC = md
SRC = md.cpp Atom.cpp System.cpp
OBJ = md.o Atom.o System.o

$(SCR): 
	$(CXX) $(CXXFLAGS) -c $(SRC)  

$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJ) $(LIBS)

clean:
	rm *.o
