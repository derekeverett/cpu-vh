DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)rhic
DIR_H          = $(DIR_MAIN)include/
DIR_BUILD      = $(DIR_MAIN)build/
DIR_OBJ        = $(DIR_BUILD)rhic

DEBUG =
OPTIMIZATION = -O3 
FLOWTRACE =
OPTIONS = -std=c++11
LINK_OPTIONS = -pthread
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS)
COMPILER = c++
LIBS = -lm -lgsl -lgslcblas -lconfig -lgtest
INCLUDES = -I rhic/rhic-core/src/include -I rhic/rhic-harness/src/main/include \
					 -I rhic/rhic-trunk/src/include -I rhic/rhic-harness/src/include \
					 -I rhic/freezeout -I rhic/freezeout/cornelius-c++-1.3

CPP := $(shell find $(DIR_SRC) -name '*.cpp')
CPP_OBJ  = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)
OBJ = $(CPP_OBJ) 

EXE =\
cpu-vh

$(EXE): $(OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES)

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	@[ -d $(DIR_OBJ) ] || find rhic/freezeout rhic/rhic-core rhic/rhic-harness rhic/rhic-trunk -type d -exec mkdir -p ./build/{} \;
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -rf $(EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); rmdir $(DIR_BUILD); fi

.SILENT:
