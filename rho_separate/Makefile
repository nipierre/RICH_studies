CC      = g++
CFLAGS  = -Wall -Wextra -O5 -g
LDFLAGS = `root-config --glibs --cflags` -lRooFit -lRooFitCore

all: fit_table

fit_table: fit_table.cc fit_table.h
		$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)
