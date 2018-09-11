DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)rhic/src
DIR_H          = $(DIR_MAIN)rhic/include
DIR_BUILD      = $(DIR_MAIN)build/
DIR_OBJ        = $(DIR_BUILD)rhic

DEBUG =
OPTIMIZATION = -O3
FLOWTRACE =

#Different options for different compilers and optimizations 
#OPTIONS = -fopenmp -march=native -fopt-info-vec #-funroll-loops #-static-libstdc++ # for g++ with vector info 
#OPTIONS = -fopenmp -march=native -std=c++11 # for g++
OPTIONS = -qopenmp -std=c++11 -lhdf5 -lhdf5_cpp # for icpc

#Need to link against libconfig, googletest should be unecessary now? 
LINK_OPTIONS = -L/home/everett.165/libconfig-1.5/lib/.libs -lconfig -L/home/everett.165/googletest-master/googletest/mybuild/ -lgtest

CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS)

#choose a compiler
#COMPILER = g++
COMPILER = icpc

#choose libs for different compilers 
#LIBS = -lm -lgsl -lgslcblas -lconfig -lgtest -lgomp # for g++ 
LIBS = -lgsl -lgslcblas -lconfig # for icpc  

INCLUDES = -I rhic/include -I /home/everett.165/libconfig-1.5/lib/ -I /home/everett.165/googletest-master/googletest/include/ -I freezeout

CPP := $(shell find $(DIR_SRC) -name '*.cpp')
CPP_OBJ  = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)
OBJ = $(CPP_OBJ)

EXE =\
cpu-vh

$(EXE): $(OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) $(OPTIONS) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES)

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	@[ -d $(DIR_OBJ) ] || find rhic -type d -exec mkdir -p ./build/{} \;
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -rf $(EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); rmdir $(DIR_BUILD); fi

.SILENT:
