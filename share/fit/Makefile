all : fitscan

LIBS = `root-config --libs` -L$(HOME)/local/lib
CXXFLAGS = `root-config --cflags` -I$(HOME)/local/include

BINDIR=../../bin
fitscan : FitOniumRCompactLib.o fitscan.o interference.o
		g++ -o $(BINDIR)/fitscan FitOniumRCompactLib.o interference.o fitscan.o   $(LIBS) -lMinuit

FitOniumRCompactLib.o : 
		g++ -o FitOniumRCompactLib.o $(CXXFLAGS) -c FitOniumRCompactLib.cc

fitscan.o : FitScan.cc
		g++ -o fitscan.o  $(CXXFLAGS) -c FitScan.cc 

interference.o :
		g++ -o interference.o  $(CXXFLAGS) -c ../interference.cpp
	
	
clean :
	rm $(BINDIR)/fitscan *.o
