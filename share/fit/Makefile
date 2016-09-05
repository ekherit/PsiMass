all : fitscan draw

LIBS = `root-config --libs` -lMinuit -L$(HOME)/local/lib  -lboost_program_options
CXXFLAGS = `root-config --cflags` -I$(HOME)/local/include  -std=c++0x 

BINDIR=../../bin
fitscan : FitOniumRCompactLib.o fitscan.o interference.o
		g++ -o $(BINDIR)/fitscan  FitOniumRCompactLib.o interference.o fitscan.o  $(LIBS)

draw : FitOniumRCompactLib.o draw.o 
		g++ -o $(BINDIR)/draw FitOniumRCompactLib.o interference.o draw.o   $(LIBS)

FitOniumRCompactLib.o : 
		g++ -o FitOniumRCompactLib.o $(CXXFLAGS) -c FitOniumRCompactLib.cc

fitscan.o : FitScan.cc
		g++ -o fitscan.o  $(CXXFLAGS) -c FitScan.cc 

interference.o :
		g++ -o interference.o  $(CXXFLAGS) -c ../interference.cpp
	
	
clean :
	rm $(BINDIR)/fitscan *.o $(BINDIR)/draw
