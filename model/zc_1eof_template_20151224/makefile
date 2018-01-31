# 
# make file for u54
CC=ifort
CFLAGS=-ipo -O3 -no-prec-div -r8 # -fast without -static and double precision 
#
#source files
SRC1=parameters.f
SRC2=continue.f
SRC3=operators.f
SRC4=output.f
SRC5=solve.f
SRC6=basis.f
COM=const.com coeff.com mat.com tau.com outfld.com
# object files
OBJ1=parameters.o
OBJ2=continue.o
OBJ3=operators.o
OBJ4=output.o
OBJ5=solve.o
OBJ6=basis.o
# NAG library location
NAG_KUSARI_FILE=/usr/local/NAG/naglm.ctl
#LIBNAG=/usr/local/NAG/fll6i22dcl/fll6i22dcl/lib
LIBNAG=/usr/local/NAG/fll6i24dcl/fll6i24dcl/lib/libnag_nag.a
#LIBS=-L$(LIBNAG) -lnag_nag
LIBS=$(LIBNAG)
EXE=run
#
# link (also link NAG libraries)
$(EXE):           $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6)
	$(CC) $(CFLAGS) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(LIBS) -o $(EXE)

$(OBJ1):        $(SRC1) $(COM)
	$(CC) $(CFLAGS) -c  $(SRC1)
$(OBJ2):        $(SRC2) $(COM) 
	$(CC) $(CFLAGS) -c $(SRC2) 
$(OBJ3):        $(SRC3) $(COM)
	$(CC) $(CFLAGS) -c  $(SRC3)
$(OBJ4):        $(SRC4) $(COM)
	$(CC) $(CFLAGS) -c  $(SRC4)
$(OBJ5):        $(SRC5) $(COM)
	$(CC) $(CFLAGS) -c  $(SRC5)
$(OBJ6):        $(SRC6) $(COM)
	$(CC) $(CFLAGS) -c  $(SRC6)

clean:
	rm *.o $(EXE)
