OPT=/opt/local/
INCDIR=-I$(PHD)/dev/include -I$(HOME)/local/include -I$(OPT)/include/ -I$(HOME)/local/arpack++/include
LIBDIR=-L$(HOME)/local/lib -L$(OPT)/lib/
LIBSSTD=-lgsl -lgslcblas -lm -larpack -llapack -lsuperlu_5.1 -lcblas -lgfortran -lconfig++
LIBSOMP=$(LIBSSTD) -lgomp

CC=g++
WARN=-Wall
CFLAGSSTD=$(WARN) -ansi -pedantic-errors -O3
CFLAGSOMP=$(CFLAGSSTD) -fopenmp
LDFLAGSSTD=-undefined dynamic_lookup
LDFLAGSGDB=$(LDFLAGSSTD) -gdb

LDFLAGS=$(LDFLAGSSTD)
#LDFLAGS=$(LDFLAGSGDB)

#EXE=transferCZ
EXE=spectrumCZ
CFLAGS=$(CFLAGSOMP)
LIBS=$(LIBSOMP)

all:$(EXE)

$(EXE).o:$(EXE).cpp
	$(CC) -c $(CFLAGS) $(INCDIR) $(EXE).cpp

$(EXE).out:$(EXE).o
	$(CC) $(LDFLAGS) -o $(EXE).out $(LIBDIR) $(EXE).o $(LIBS)

clean:
	rm -f $(EXE).o $(EXE).out
