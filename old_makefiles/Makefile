DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)rhic/src
DIR_H          = $(DIR_MAIN)rhic/include
DIR_BUILD      = $(DIR_MAIN)build/
DIR_OBJ        = $(DIR_BUILD)rhic

DEBUG =
OPTIMIZATION = -O3
FLOWTRACE =
OPTIONS = -fopenmp -march=native -fopt-info-vec #-funroll-loops #-static-libstdc++
#OPTIONS = -qopenmp -xavx #-static-libstdc++
LINK_OPTIONS = -L/home/everett.165/libconfig-1.5/lib/.libs -lconfig -L/home/everett.165/googletest-master/googletest/mybuild/ -lgtest
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS)
COMPILER = g++
#COMPILER = icpc
LIBS = -lm -lgsl -lgslcblas -lconfig -lgtest -lgomp
#LIBS = -lconfig -lgtest
INCLUDES = -I rhic/include -I /home/everett.165/libconfig-1.5/lib/ -I /home/everett.165/googletest-master/googletest/include/ -I freezeout

CPP := $(shell find $(DIR_SRC) -name '*.cpp')
CPP_OBJ  = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)
OBJ = $(CPP_OBJ)

EXE =\
cpu-vh

$(EXE): $(OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES)

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	@[ -d $(DIR_OBJ) ] || find rhic -type d -exec mkdir -p ./build/{} \;
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -rf $(EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); rmdir $(DIR_BUILD); fi

.SILENT:
