SYSTEM     = x86-64_osx
LIBFORMAT  = static_pic
FILENAME = ccdual-decomp2

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /Users/administrator/Applications/IBM/ILOG/CPLEX_Studio1251/cplex
CONCERTDIR    = /Users/administrator/Applications/IBM/ILOG/CPLEX_Studio1251/concert
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = clang++ -O0
#CC  = gcc -O0
CC = clang -O0
JAVAC = javac 

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -stdlib=libstdc++
COPT  = -m64 -fPIC 
JOPT  = -classpath $(CPLEXDIR)/lib/cplex.jar -O

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CLUSTERDIR = /Users/administrator/Documents/ysong-luedtke/ccpartition/cluster-1.52/src

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CLNDIRS   = -L$(CPLEXLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -m64 -lm -lpthread -framework CoreFoundation -framework IOKit -stdlib=libstdc++
CLNFLAGS  = -lcplex -m64 -lm -lpthread -framework CoreFoundation -framework IOKit
JAVA      = java -d64 -Djava.library.path=$(CPLEXDIR)/bin/x86-64_darwin -classpath $(CPLEXJARDIR):

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 
JCFLAGS = $(JOPT)


#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------

all_cpp: $(FILENAME)

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf $(C_EX) $(CX_EX) $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp
EIGENDIR = /Users/administrator/Documents/ysong-luedtke/ccpartition/Eigen
CLUSTERDIR = /Users/administrator/Documents/ysong-luedtke/ccpartition/cluster-1.52/src
$(FILENAME): cluster.o $(FILENAME).o 
	$(CCC) $(CCFLAGS) $(CCLNDIRS) cluster.o $(FILENAME).o  -o $(FILENAME) $(CCLNFLAGS)
$(FILENAME).o: $(FILENAME).cpp $(CLUSTERDIR)/cluster.h
	$(CCC) -c $(CCFLAGS) -I$(EIGENDIR) -I$(CLUSTERDIR) $(FILENAME).cpp -o $(FILENAME).o
cluster.o: $(CLUSTERDIR)/cluster.c $(CLUSTERDIR)/cluster.h
	$(CC) -c -Wall $(CLUSTERDIR)/cluster.c -o cluster.o
# Local Variables:
# mode: makefile
# End:

#$(FILENAME): $(FILENAME).o cluster.o
#	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $(FILENAME) $(FILENAME).o $(CCLNFLAGS)
#$(FILENAME).o: $(FILENAME).cpp cluster.h
#	$(CCC) -c $(CCFLAGS) -I$(EIGENDIR) $(FILENAME).cpp -o $(FILENAME).o

