all : fitscan

LIBS = `root-config --libs`
CXXFLAGS = `root-config --cflags`

fitscan : FitOniumRCompactLib.o fitscan.o
		g++ -o fitscan FitOniumRCompactLib.o fitscan.o $(LIBS) -lMinuit

FitOniumRCompactLib.o : 
		g++ -o FitOniumRCompactLib.o $(CXXFLAGS) -c FitOniumRCompactLib.cc

fitscan.o :
		g++ -o fitscan.o  $(CXXFLAGS) -c FitScan.cc 
	
	
clean :
	rm fitscan *.o
