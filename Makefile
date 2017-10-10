#Makefile for hsol 
#C++ compiler 
#CCP      = g++ -D JAJH_dev -D JAJH_rp
CCP      = g++ -D SCIP_dev 

# For up-to-date gnu compiler (mingw-32 on Fetteresk)
CCSTDFLAGS= -std=c++11 -fPIC

# Default C++ macros
# For debugging, need -g, but other flags like -Wall -Werror can be valuable
CCFLAGS = $(CCSTDFLAGS) -ggdb3 
LDFLAGS = 
LINK    = $(CCP)

# Define the name for the output with alternative extensions for special
# run modes
# Default mode is debug (made by make all)
EXECUTABLE       = hsol
DEBUG            = $(EXECUTABLE:%=./%_deb.exe)
# Default executable (for test) 
EXE=$(DEBUG) 
# Default executable (for make all) (exe/deb)
all: deb test

SOURCE_CPP_H=HDual.cpp HDualRHS.cpp HDualRow.cpp HFactor.cpp HMatrix.cpp HModel.cpp HPrimal.cpp HTester.cpp HVector.cpp \
	HCrash.cpp \
	HPresolve.cpp HPreData.cpp \
	KktCheck.cpp KktChStep.cpp \
	HinOut.cpp

SOURCE_CPP=HApp.cpp HDualMulti.cpp #$(SOURCE_CPP_H:%=HApp.cpp HDualMulti.cpp)

OBJECTS_H= $(SOURCE_CPP_H:%.cpp=%.o) 
OBJECTS_NO_HEADER=  $(SOURCE_CPP:%.cpp=%.o) 

deb : $(OBJECTS_H) $(OBJECTS_NO_HEADER) 
	$(LINK) $(CCFLAGS) $(OBJECTS_H) $(OBJECTS_NO_HEADER) -o $(DEBUG)
	EXE=$(DEBUG)

exe: $(OBJECTS_H) $(OBJECTS_NO_HEADER) 
	CCFLAGS=$(CCSTDFLAGS)
	EXE=$(EXECUTABLE)
	$(LINK) $(CCFLAGS) $(OBJECTS_H) $(OBJECTS_NO_HEADER) -o $(EXECUTABLE:%=%.exe)

lib: $(OBJECTS_H) $(OBJECTS_NO_HEADER)
	-rm -f $(LIBFILE)
	ar cr -o lib$(EXECUTABLE).a $(OBJECTS_H) $(OBJECTS_NO_HEADER)
	ranlib lib$(EXECUTABLE).a

# Report some environment variables
$(info $$DEBUG      is [${DEBUG}])
$(info $$LINK       is [${LINK}])
$(info $$CCFLAGS    is [${CCFLAGS}])
#$(info $$OBJECTS    is [${OBJECTS}])
#$(info $$SOURCE_CPP is [${SOURCE_CPP}])

# State how to compile the objects
#$(OBJECTS):  $(OBJECTS_H) # $(OBJECTS_NO_HEADER)

$(OBJECTS_H): %.o : %.cpp %.h
	$(CCP) $(CCFLAGS) -c $(DIR)$(@F:.o=.cpp)

$(OBJECTS_NO_HEADER): %.o : %.cpp
	$(CCP) $(CCFLAGS) -c $(DIR)$(@F:.o=.cpp)

test:
	$(EXE) -f 25fv47.mps

# "Cleaning the Directory": for .PHONY logic see GNU make documentation
# Make sure that rm is known: not in mingw-32!
.PHONY: clean cleanall
clean:
	-rm -f $(OBJECTS_H) $(OBJECTS_NO_HEADER) 

cleanall: 
	rm -f $(EXE) $(DEBUG) $(OBJECTS_H) $(OBJECTS_NO_HEADER)  $(EXECUTABLE:%=%.exe) *.o 

.DONE:          # performed last, unless errors occur so .FAILED defined
.FAILED:        # performed instead of .DONE if error occurs
		-@echo 
