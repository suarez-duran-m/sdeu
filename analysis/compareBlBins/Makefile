SRC_DIR=$(PWD)/
OBJ_DIR=$(PWD)/
LIB_DIR=$(PWD)/
BIN_DIR=$(PWD)/

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --nonew --glibs)

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
I_RECINC= -I$(CDASHOME)/include/ 
I_RECINC+= -I$(FD_ROOT)/include/


all:   $(BIN_DIR)/analyseSDEUbins

OBJS  = $(patsubst $(SRC_DIR)/%LinkDef.h, $(OBJ_DIR)/%Dict.o, $(wildcard $(SRC_DIR)/*LinkDef.h))
OBJS += $(patsubst $(SRC_DIR)/%.cc, $(OBJ_DIR)/%.o, $(wildcard $(SRC_DIR)/*.cc))



$(BIN_DIR)/analyseSDEUbins: $(OBJS)
	@mkdir -p $(BIN_DIR)
	@echo -n "making exe..." 
	$(LD) $(LDFLAGS) $(OBJS) -o $@ $(EXPLLINKLIBS) -I$(ROOTSYS)/include/  $(L_RECINC) -L$(ROOTSYS)/lib -lMinuit


$(SRC_DIR)/%Dict.cc: $(SRC_DIR)/%.h  $(SRC_DIR)/%LinkDef.h
	@(echo generating $@ dictionary)
	rootcint -f $@ -c $(ROOTCFLAGS) $(I_RECINC) $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc 
	@mkdir -p $(OBJ_DIR)
	$(CXX) -c $(CXXFLAGS) $(ROOTCFLAGS)  $(I_RECINC) \
		-c $^ -o $@

#dependencies
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(SRC_DIR)/%.h

clean:
	@rm -f *~ $(OBJS) $(BIN_DIR)/analyseSDEUbins


