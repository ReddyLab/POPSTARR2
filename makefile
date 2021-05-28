CC		= g++
DEBUG		= -g
OPTIMIZE	= -O
CFLAGS		= $(OPTIMIZE) -fpermissive -w
LDFLAGS		= $(OPTIMIZE)
BOOM		= BOOM
OBJ		= obj
LIBS		= -LBOOM -lBOOM -lgsl -lm -lgslcblas
#---------------------------------------------------------
$(OBJ)/make-haplotypes.o:\
		make-haplotypes.C
	$(CC) $(CFLAGS) -o $(OBJ)/make-haplotypes.o -c \
		make-haplotypes.C
#---------------------------------------------------------
make-haplotypes: \
		$(OBJ)/make-haplotypes.o
	$(CC) $(LDFLAGS) -o make-haplotypes \
		$(OBJ)/make-haplotypes.o \
		$(LIBS)
#---------------------------------------------
