CXX = c++
# Flags for optimized execs
# CXXFLAGS = O3
CXXFLAGS = 
# Flags for debugging
#CXXFLAGS = -ftemplate-depth-30 -g -DBZ_DEBUG

LIBS = -lm -larmadillo

EXEC = saw
SRC = saw.cpp
OBJ = saw.o lib.o

$(SCR): 
	$(CXX) $(CXXFLAGS) -c $(SRC)  

$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJ) $(LIBS)

clean:
	rm *.o
