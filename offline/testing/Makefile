# General Makefile for the OfflineUser package
#
# If the environment variable AUGEROFFLINEROOT is not set
# AND the executable 'auger-offline-config' is not in your PATH
# the definition of the following variable is required
#
# AUGEROFFLINEROOT := @prefix@
#
# Replace the wildcard expression with .cc file list if you do
# not want to compile all .cc files in this directory
#
USER_SRCS := $(wildcard *.cc)
#
# All .xml.in files will be transformed into .xml with correct
# config and schema locations
#
USER_XMLS := $(patsubst %.xml.in,%.xml,$(wildcard *.xml.in))
#
# Give your executable a name
#
EXE := userAugerOffline
#
#############################################################

## You should not need to change anything below this line ###

# Authors: T. Paul, S. Argiro, L. Nellen, D. Veberic
# $Id: Makefile.in 27673 2015-05-16 11:42:53Z darko $
# Send bug reports to http://www.auger.unam.mx/bugzilla/

.PHONY: all depend clean

ifdef AUGEROFFLINEROOT
  AUGEROFFLINECONFIG := $(AUGEROFFLINEROOT)/bin/auger-offline-config
else
  AUGEROFFLINECONFIG := auger-offline-config
  AUGEROFFLINEROOT := $(shell $(AUGEROFFLINECONFIG) --install)
endif

OBJS := $(USER_SRCS:.cc=.o)

CPPFLAGS    := $(shell $(AUGEROFFLINECONFIG) --cppflags)
CXXFLAGS    := $(shell $(AUGEROFFLINECONFIG) --cxxflags)
LDFLAGS     := $(shell $(AUGEROFFLINECONFIG) --ldflags) -lMinuit
MAIN        := $(shell $(AUGEROFFLINECONFIG) --main)
CONFIGFILES := $(shell $(AUGEROFFLINECONFIG) --config)
XMLSCHEMALOCATION := $(shell $(AUGEROFFLINECONFIG) --schema-location)

ifndef INPUTDATA
  INPUTTYPE := IoAuger 
  INPUTDATA := $(shell find $(HOME)/Documents/augerData/2020/12 -iname "ad_*.root")
endif

DATABASESERVERS := $(if $(wildcard databaseServers.*),<configLink id="databaseServers" type="XML" xlink:href="databaseServers.xml"/>,)

export INPUTTYPE
export INPUTDATA
export DATABASESERVERS
export CONFIGDIR
export CONFIGFILES
export XMLSCHEMALOCATION

all: $(EXE) $(USER_XMLS)

config: $(USER_XMLS)

$(EXE): $(OBJS)
	$(CXX) -o $@ $^ $(MAIN) $(CXXFLAGS) $(LDFLAGS)

%: %.in
	@echo "Generating $@ file..."
	@perl -pe 's/\@(INPUTDATA|INPUTTYPE|INPUTUPTIME)\@/$$ENV{$$1}/g;s/@(CONFIGDIR)@/$$ENV{CONFIGFILES}/g;s/@(SCHEMALOCATION)@/$$ENV{XMLSCHEMALOCATION}/g;s/<!--\s*PLACEHOLDER\s+databaseServers\s*-->/$$ENV{DATABASESERVERS}/g;' $< >$@

clean:
	- rm -f *.o  *.so *.ps core $(USER_XMLS) Make-depend

#############################################################
# 'make run' will run the thing

run: $(EXE) $(USER_XMLS)
	./$(EXE) -b bootstrap.xml

#############################################################
# gcc can generate the dependency list

depend: Make-depend

Make-depend: $(USER_SRCS)
	$(CPP) $(CPPFLAGS) -MM $^ >$@

-include Make-depend
