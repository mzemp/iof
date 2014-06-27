# Makefile for iof library

NAME	= iof
VERSION	= $(shell git describe --tags --long)

CC		= gcc
CFLAGS	= -O3 -mcmodel=medium -Wall -pedantic
LIBS	=

SRCS	= $(wildcard *.c)
INCS	= $(wildcard *.h)
OBJS	= $(SRCS:.c=.o)

# Rules

$(NAME): $(OBJS) $(INCS) Makefile
	ar rcs lib$(NAME).a $(OBJS) $(LIBS)

clean:
	rm -f *~ *.o *.a

install:
	cd ../include; ln -sf ../$(NAME)/$(NAME).h .
	cd ../lib; ln -sf ../$(NAME)/lib$(NAME).a .
