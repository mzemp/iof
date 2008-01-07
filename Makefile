# Makefile for IOfunctions

BASE	= iof
VERSION = 1.0

CC	= gcc
CFLAGS	= -O3 -Wall
LIBS	=

# Rules

all:	$(BASE).h Makefile
	$(CC) $(CFLAGS) -c -o $(BASE).o $(BASE).c
	ar rcs lib$(BASE).a $(BASE).o

clean:
	-rm -f *~ *.o *.a

tar:
	cd ..; tar cvf - $(BASE)/*.c $(BASE)/*.h $(BASE)/Makefile > $(BASE)-$(VERSION).tar
