INCDIR=-I$(HOME)/.local/include/
LIBDIR=-L$(HOME)/.local/lib/
LIBS=-lergopack -larpack -lgsl -lgslcblas -lconfig++ -lgomp -lgfortran -lopenblas

CC=g++
WARN=-Wall -Wformat=0
CFLAGS=$(WARN) -O3
LDFLAGSSTD=

LDFLAGS=$(LDFLAGSSTD)

EXE=spectrumCZ

all:$(EXE).out
	mkdir -p ../results/spectrum/eigval
	mkdir -p ../results/spectrum/eigvec
	mkdir -p ../results/spectrum/conditionNumbers


$(EXE).o:$(EXE).cpp 
	$(CC) -c $(CFLAGS) $(INCDIR) -o $@ $<

$(EXE).out:$(EXE).o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBDIR) $(LIBS)

clean:
	rm -f $(EXE).o $(EXE).out

cleanall:
	rm -f $(EXE).o $(EXE).out
