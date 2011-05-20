all : fitscan

LIBS = `root-config --libs` -L$(HOME)/local/lib
CXXFLAGS = `root-config --cflags` -I$(HOME)/local/include

BINDIR=../../bin
fitscan : FitOniumRCompactLib.o fitscan.o
		g++ -o $(BINDIR)/fitscan FitOniumRCompactLib.o fitscan.o $(LIBS) -lMinuit

FitOniumRCompactLib.o : 
		g++ -o FitOniumRCompactLib.o $(CXXFLAGS) -c FitOniumRCompactLib.cc

fitscan.o : FitScan.cc
		g++ -o fitscan.o  $(CXXFLAGS) -c FitScan.cc 
	
	
clean :
	rm $(BINDIR)/fitscan *.o
