CC      = g++
CFLAGS  = -Wall -Wextra -O5 -g
LDFLAGS = `root-config --glibs --cflags`

all: draw_rich_fig

draw_rich_fig: draw_rich_fig.cc
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)  -lHistPainter

create_table: create_table.cc
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) 
