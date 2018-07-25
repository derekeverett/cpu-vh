DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)rhic/src
DIR_H          = $(DIR_MAIN)rhic/include
DIR_BUILD      = $(DIR_MAIN)build/
DIR_OBJ        = $(DIR_BUILD)rhic

DEBUG =
OPTIMIZATION = -O3
FLOWTRACE =
OPTIONS = -fopenmp -march=native -std=c++11 #-funroll-loops #-static-libstdc++
#OPTIONS = -qopenmp -xavx #-static-libstdc++
LINK_OPTIONS = -L/home/du.458/libconfig/lib/.libs -lconfig -L/home/du.458/googletest/googletest/mybuild/ -lgtest
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS)
COMPILER = g++
#COMPILER = icpc
LIBS = -lm -lgsl -lgslcblas -lconfig -lgtest -lgomp
#LIBS = -lconfig -lgtest
INCLUDES = -I rhic/include -I /home/du.458/libconfig/lib/ -I /home/du.458/googletest/googletest/include/ -I freezeout

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
