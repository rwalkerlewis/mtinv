CC = clang
CFLAGS = -Bstatic -I../../include
INSTALL = ../../bin

LIBS = -lget -lm -lglib
LIBDIR = -L../../lib

all : mtcomp mtdecomp mtdecomp2

SUBS = ../mtinv_subs.o ../tred2.o ../tqli.o ../trans1.o ../tpdss.o \
       ../mtinv_gmtsubs.o ../indexx.o ../eig2lune.o ../ascii_mech.o 

mtcomp : mtcomp.o $(SUBS)
	$(CC) $(CFLAGS) -o $@ mtcomp.o $(SUBS) $(LIBDIR) $(LIBS)
	mv $@ $(INSTALL)

mtdecomp : mtdecomp.o $(SUBS)
	$(CC) $(CFLAGS) -o $@ mtdecomp.o $(SUBS) $(LIBDIR) $(LIBS)
	mv $@ $(INSTALL)

mtdecomp2 : mtdecomp2.o $(SUBS)
	$(CC) $(CFLAGS) -o $@ mtdecomp2.o $(SUBS) $(LIBDIR) $(LIBS)
	mv $@ $(INSTALL)

clean :
	/bin/rm -f *.o mtdecomp
