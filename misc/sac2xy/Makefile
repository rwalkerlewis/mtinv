include ../../Makefile.inc

#CC = clang
#CC = gcc

INCLUDE = ../../include
CFLAGS = -I$(INCLUDE)
LIBS = -L../../lib -lget -lm
INSTALL = ../../bin

# all : sac2xy yoffset sacsource
all : sac2xy

.c.o :
	$(CC) $(CFLAGS) -c $< -o $@

sac2xy : sac2xy.o ../../src/timesubs.o version.h
	$(CC) $(CFLAGS) -o sac2xy sac2xy.o ../../src/timesubs.o $(LIBS)
	mv $@ $(INSTALL)

yoffset : yoffset.o
	$(CC) $(CFLAGS) -o yoffset yoffset.o $(LIBS)
	mv $@ $(INSTALL)

sacsource : sacsource.o
	$(CC) $(CFLAGS) -o $@ sacsource.o $(LIBS)
	mv $@ $(INSTALL)

clean :
	rm -f *.o sac2xy yoffset sacsource
