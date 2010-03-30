# Makefile for iof library

NAME	= iof
SOURCES = iof_auxiliary.c iof_tipsy.c iof_gadget.c iof_art.c iof_array.c	
VERSION = 1.9

CC	= gcc
CFLAGS	= -O3 -mcmodel=medium -Wall -pedantic
LIBS	=

# Rules

iof:	$(SOURCES:.c=.o) $(SOURCES:.c=.h) Makefile
	ar rcs lib$(NAME).a $(SOURCES:.c=.o) $(LIBS)

clean:
	-rm -f *~ *.o *.a

install:
	cd ../include; ln -sf ../$(NAME)/$(NAME).h .
	cd ../lib; ln -sf ../$(NAME)/lib$(NAME).a .

tar:
	cd ..; tar cvf - $(NAME)/*.c $(NAME)/*.h $(NAME)/Makefile > $(NAME)-$(VERSION).tar
