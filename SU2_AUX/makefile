#
# Overall makefile for SU2_AUX
# Author: F. Palacios
# Version: 1.0
# 

MAIN=			SU2_AUX

include $(SU2_HOME)/$(MAIN)/makefile.in

INC_COMMON=		$(SU2_HOME)/Common/include
SRC_COMMON=		$(SU2_HOME)/Common/src
INCLUDES_COMMON=	$(INC_COMMON)/config_structure.hpp $(INC_COMMON)/dual_grid_structure.hpp \
			$(INC_COMMON)/geometry_structure.hpp  $(INC_COMMON)/primal_grid_structure.hpp  \
			$(INC_COMMON)/grid_movement_structure.hpp $(INC_COMMON)/matrix_structure.hpp \
			$(INC_COMMON)/vector_structure.hpp $(INC_COMMON)/linear_solvers_structure.hpp \
                       $(INC_COMMON)/option_structure.hpp
SOURCES_COMMON=		$(SRC_COMMON)/config_structure.cpp $(SRC_COMMON)/dual_grid_structure.cpp \
			$(SRC_COMMON)/geometry_structure.cpp $(SRC_COMMON)/primal_grid_structure.cpp \
			$(SRC_COMMON)/grid_movement_structure.cpp $(SRC_COMMON)/matrix_structure.cpp \
			$(SRC_COMMON)/vector_structure.cpp $(SRC_COMMON)/linear_solvers_structure.cpp

SRC_AUX=		$(SU2_HOME)/$(MAIN)/src
INC_AUX=		$(SU2_HOME)/$(MAIN)/include
INCLUDES_AUX=		$(wildcard $(INC_AUX)/*.hpp)
SOURCES_AUX=		$(wildcard $(SRC_AUX)/*.cpp)

OBJ=			$(SU2_HOME)/$(MAIN)/obj
OBJECTS_COMMON=		$(patsubst %.cpp, $(OBJ)/%.o, $(notdir $(SOURCES_COMMON)))
OBJECTS_AUX=		$(patsubst %.cpp, $(OBJ)/%.o, $(notdir $(SOURCES_AUX)))

BIN=			$(SU2_HOME)/$(MAIN)/bin
BINARIES=		$(BIN)/$(MAIN)

.SECONDARY:
.SUFFIXES:
.PHONY: default help clean all

default: help

help:
	@echo 'You must specify a target between: make all, make clean '

all: $(BINARIES)

$(OBJ)/%.o: $(SRC_COMMON)/%.cpp $(INCLUDES_COMMON) makefile
	@echo "Compiling \""$@"\" from \""$<"\""
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_AUX) $(INC_CGNS) -o $@ -c $<

$(OBJ)/%.o: $(SRC_AUX)/%.cpp $(INCLUDES_AUX) makefile
	@echo "Compiling \""$@"\" from \""$<"\""
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_AUX) $(INC_CGNS) -o $@ -c $<

$(BIN)/%: $(OBJECTS_COMMON) $(OBJECTS_AUX)
	@echo "Creating "$@
	@$(CXX) $(CXXFLAGS) -I./$(INC_COMMON) -I./$(INC_AUX) $(INC_CGNS) $(LIB_CGNS) -o $@ $(OBJECTS_COMMON) $(OBJECTS_AUX) $(LDFLAGS)

clean:
	@echo "Deleting temporary files"
	@-rm -f core 
	@-rm -f `(find . | grep "~")`
	@echo "Deleting objects"
	@-rm -f *.o $(OBJECTS_AUX)
	@-rm -f *.o $(OBJECTS_COMMON)
	@echo "Deleting binaries"
	@-rm -f $(BINARIES)
