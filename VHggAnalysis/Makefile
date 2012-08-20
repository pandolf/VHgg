# $Id: Makefile,v 1.3 2012/08/20 12:41:45 pandolf Exp $
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

BINFILES = $(wildcard *.cpp)

PROGRAMS = $(patsubst %.cpp,%,$(BINFILES))


# --- External configuration ---------------------------------
CC         = g++ -Wall
CCFLAGS    =  -g 
MFLAGS     = -MM
INCLUDES   =
WORKDIR    = $(PWD)
LIBDIR     = $(WORKDIR)
OBJDIR=$(WORKDIR)/objects/
# -------------------------------------------------------------

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')
ROOFIT_LIBS := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIB= | sed 's|LIB=||')
ROOFIT_LIBS += $(shell cd $(CMSSW_BASE); scram tool info roofit | grep LIB= | sed 's|LIB=||') 


INCLUDES += -I.  -I.. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/
ROOTSYS  ?= ERROR_RootSysIsNotDefined

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit -lRIO -lTMVA

# CC files excluding the binaries
#CCFILES=$(filter-out $(BINFILES),$(wildcard $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/*.cc))
CCFILES=$(filter-out $(BINFILES),$(wildcard *.cc))

# List of all object files to build
OLIST = $(OBJDIR)/RedNtpFinalizer.o $(OBJDIR)/AnalysisJet.o $(OBJDIR)/BTagSFUtil.o $(OBJDIR)/SFlightFuncs.o $(OBJDIR)/MistagFuncs.o
OLIST +=$(patsubst %.cc,$(OBJDIR)/%.o,$(notdir $(CCFILES)))

$(OBJDIR)/RedNtpFinalizer.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/RedNtpFinalizer.cc
	@echo "Compiling $<"
	mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) -c $< -o $@ $(INCLUDES)

$(OBJDIR)/AnalysisJet.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisJet.cc
	@echo "Compiling $<"
	mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) -c $< -o $@ $(INCLUDES)

$(OBJDIR)/BTagSFUtil.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/BTagSFUtil.cc
	@echo "Compiling $<"
	mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) -c $< -o $@ $(INCLUDES)

$(OBJDIR)/SFlightFuncs.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/SFlightFuncs.cc
	@echo "Compiling $<"
	mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) -c $< -o $@ $(INCLUDES)

$(OBJDIR)/MistagFuncs.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/MistagFuncs.cc
	@echo "Compiling $<"
	mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) -c $< -o $@ $(INCLUDES)

# Implicit rule to compile all classes
$(OBJDIR)/%.o : %.cc
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) -c $< -o $@ $(INCLUDES)


$(PROGRAMS) : $(OLIST)
	@echo "Linking $@"
	@$(CC) $(CCFLAGS)  $(INCLUDES) $(OLIST) \
	$(ROOTLIBS) $(EXTRALIBS) -o $(WORKDIR)/$@  $(patsubst %,%.cpp,$@)

default : analysis 

all : ${PROGRAMS}

clean:
	rm -Rf $(OBJDIR)/*.o


#.PHONY : shlib lib default clean
