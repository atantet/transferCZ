INCDIR=-I$(HOME)/.local/include/
LIBDIR=-L$(HOME)/.local/lib/
LIBS=-lgsl -lconfig++ -lergopack -lgslcblas -larpack

CC=g++
WARN=-Wall -Wformat=0
CFLAGS=$(WARN) -pedantic-errors -O3
LDFLAGSSTD=
LDFLAGSGDB=$(LDFLAGSSTD) -gdb
WITH_OMP=1

# Use OpenMP?
ifeq ($(WITH_OMP),1)
  LIBS +=-L/opt/local/lib/libgcc/ -lgomp
  CFLAGS +=-fopenmp -DWITH_OMP=$(WITH_OMP)
endif

LDFLAGS=$(LDFLAGSSTD)
#LDFLAGS=$(LDFLAGSGDB)

EXE=transferCZ

all:$(EXE).out
	mkdir -p ../results/transfer/transition
	mkdir -p ../results/transfer/initDist
	mkdir -p ../results/transfer/gridMem
	mkdir -p ../results/transfer/mask
	mkdir -p ../results/grid/

$(EXE).o:$(EXE).cpp 
	$(CC) -c $(CFLAGS) $(INCDIR) -o $@ $<

$(EXE).out:$(EXE).o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBDIR) $(LIBS)

clean:
	rm -f $(EXE).o $(EXE).out

cleanall:
	rm -f $(EXE).o $(EXE).out
