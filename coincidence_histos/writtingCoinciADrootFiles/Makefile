USER_SRCS = $(wildcard *.cc)
EXE = $(patsubst %.cc,%, $(wildcard *.cc))

.PHONY: all depend clean

ROOT_CXXFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOT_LDFLAGS     = $(shell $(ROOTSYS)/bin/root-config --ldflags)
ROOT_LDFLAGS    += $(shell $(ROOTSYS)/bin/root-config --libs)

DEBUG         =  -O
CXXFLAGS      = -ggdb3 -Wall -fPIC -D_REENTRANT
LD            = g++
LDFLAGS       = -ggdb3
SOFLAGS       = -shared
EXPLLINKLIBS  = $(ROOTGLIBS) -lMinuit
# C++ flags
ifndef CXX
  CXX           = g++
endif

# Rec stuff
#L_RECINC= -L$(ADSTROOT)/lib -lRecEventKG
L_RECINC= -L$(CDASHOME)/lib -lEc -lIoAuger -lIoSd
#I_RECINC= -I$(ADSTROOT)/include/adst
I_RECINC+= -I$(CDASHOME)/include/ 
I_RECINC+= -I$(FD_ROOT)/include/


## Local libraries
HEADER = include/*.h
CFLAGS = -Iinclude

all: $(EXE)

%: %.cc src/*.cc $(HEADER)
	$(CXX) $^ $(ROOT_CXXFLAGS) $(I_RECINC) $(L_RECINC) $(ROOT_LDFLAGS) $(CFLAGS) -lMinuit -o $@


depend: Make-depend

Make-depend: $(USER_SRCS)
	$(CPP) $(ROOT_CXXFLAGS) $(I_RECINC) $(CFLAGS) -MM $^ > $@

clean:
	- rm -f *.o  *.so *.ps core Make-depend $(EXE)

-include Make-depend
