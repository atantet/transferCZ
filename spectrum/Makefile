INCDIR=-I$(HOME)/.local/include/ -I/opt/local/include/
LIBDIR=-L$(HOME)/.local/lib/ -L/opt/local/lib/
LIBS=-lgsl -lconfig++ -lopenblas -larpack -lergopack

CC=g++
WARN=-Wall -Wformat=0
CFLAGS=$(WARN) -pedantic-errors -O3
LDFLAGSSTD=-undefined dynamic_lookup
LDFLAGSGDB=$(LDFLAGSSTD) -gdb

LDFLAGS=$(LDFLAGSSTD)
#LDFLAGS=$(LDFLAGSGDB)

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
