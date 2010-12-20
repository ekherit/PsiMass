all : fitscan

LIBS = `root-config --libs`
CXXFLAGS = `root-config --cflags`

BINDIR=../../bin
fitscan : FitOniumRCompactLib.o fitscan.o
		g++ -o $(BINDIR)/fitscan FitOniumRCompactLib.o fitscan.o $(LIBS) -lMinuit

FitOniumRCompactLib.o : 
		g++ -o FitOniumRCompactLib.o $(CXXFLAGS) -c FitOniumRCompactLib.cc

fitscan.o : FitScan.cc
		g++ -o fitscan.o  $(CXXFLAGS) -c FitScan.cc 
	
	
clean :
	rm $(BINDIR)/fitscan *.o
