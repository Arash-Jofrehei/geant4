CC = g++
CFLAGS = -Wall -c -g 

INCDIR   = interface
SRCDIR   = src
OBJDIR   = lib
BINDIR   = bin
TESTDIR  = codes/makeBin


SOURCES  := $(wildcard $(SRCDIR)/*.cc)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

INCLUDES = -I. -I$(ROOTSYS)/include -I../ -I./interface# -I/shome/micheli/H4_2015/CommonTools/interface 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib -Wl,--no-as-needed $(ROOTGLIBS) -lHtml -lMathCore -lGenVector -lMinuit -lTMVA -lRooFit `xml2-config --cflags` `xml2-config --libs` -std=c++11


CPP =  $(wildcard codes/makeBin/*.cpp)
CPPOBJ = $(patsubst codes/makeBin/%,%,$(CPP))
BINOBJ	=$(patsubst %.cpp,bin/%,$(CPPOBJ))
COMMONTOOLSDIR=/shome/micheli/H4_2015/CommonTools/src/
COMMONTOOLSDIRINC=/shome/micheli/H4_2015/CommonTools/interface
COMMONTOOLS :=$(wildcard /shome/micheli/H4_2015/CommonTools/src/*.cc)
COMMONTOOLSOBJ =$(patsubst /shome/micheli/H4_2015/CommonTools/src/%.cc,lib/%.o,$(COMMONTOOLS))


all: $(BINOBJ) 

$(BINDIR)/% : codes/makeBin/%.cpp $(OBJECTS)  $(COMMONTOOLSOBJ)  
	@echo "Compiling $@"
	@echo $(CC) $(INCLUDES) -o $@ $<  $(ROOTFLAG) $(EXTRALIBS) $(OBJECTS) $(COMMONTOOLSOBJ)	
	@$(CC) $(INCLUDES) -o $@ $<  $(ROOTFLAG) $(EXTRALIBS) $(OBJECTS) $(COMMONTOOLSOBJ)	

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cc  $(INCDIR)/%.h  
	@echo "Compiling "$<" "
	@echo $(CC) $(CFLAGS) $(INCLUDES) -o $@  -c $< $(ROOTFLAG) $(EXTRALIBS) 
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@  -c $< $(ROOTFLAG) $(EXTRALIBS)	


$(COMMONTOOLSOBJ): $(COMMONTOOLSDIR)/*.cc 
	@echo "Compiling "$@" "
	@echo $(CC) $(CFLAGS)  -o $@  -c $(patsubst lib/%.o,$(COMMONTOOLSDIR)/%.cc,$@) $(ROOTFLAG) $(EXTRALIBS)
	@$(CC) $(CFLAGS)  -o $@  -c $(patsubst lib/%.o,$(COMMONTOOLSDIR)/%.cc,$@) $(ROOTFLAG) $(EXTRALIBS)	





clean:
	rm bin/*.o;rm lib/*.o
