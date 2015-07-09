ROOT_LIBS=`root-config --libs --glibs` -lRooFit -lRooFitCore
ROOT_INCS=`root-config --cflags`

INCDIR=./include/
SRCDIR=./src/
BINDIR=./bin/
VPATH=$(SRCDIR):$(INCDIR)
CC = g++ $(CFLAGS) $(ROOT_INCS) -I$(INCDIR)

all: bin/ATrack.exe
all: bin/combine_AmupiTrack_1D.exe
#all: bin/combine_AmupiTrack_2D.exe

clean:
	rm bin/*.exe
purge:
	rm -rf bin	
##Rule to build %.exe                      
bin/%.exe: src/%.cpp
	$(CC) -o $@ $(ROOT_LIBS) $(ROOT_INCS) $^

#### GENERAL RULE TO BUILD LIBRARIES ######
libs/%.o: %.cpp %.h
	$(CC) -c $< -o $@
###########################################
