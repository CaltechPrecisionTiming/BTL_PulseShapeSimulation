CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)
INC = $(shell pwd)

CPPFLAGS := $(shell root-config --cflags) -I$(INC)/include
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR) -lRooFit -lRooFitCore
CPPFLAGS += -g -std=c++14 -I$(INC)/include

TARGET = SimulatePulseShape

SRC = app/SimulatePulseShape.cc src/PulseShape.cc
OBJ = $(SRC:.cc=.o)


TARGETS = $(TARGET) $(TARGET1) $(TARGET2)

all : $(TARGETS)

$(TARGET) : $(OBJ)
	@echo $@
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)


%.o : %.cc
	@echo $@
	$(CXX) $(CPPFLAGS) -o $@ -c $<
clean :
	rm -rf *.o app/*.o src/*.o $(TARGET) *~ *.dSYM
