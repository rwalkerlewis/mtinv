CC	= gcc
FC	= gfortran

INCLUDE = -I../../include
CFLAGS  = $(INCLUDE)
INSTALL = ../../bin

.c.o :
	$(CC) $(CFLAGS) -c $< -o $@

all : wiggins

wiggins :  wiggins.o wiggins_sub.o ../sacio/sacio.o ../timesubs.o
	$(CC) $(CFLAGS) -o $@ wiggins.o wiggins_sub.o ../sacio/sacio.o ../timesubs.o -lm -L../../lib -lget

clean :
	rm -rf *.o wiggins
