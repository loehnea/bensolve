CC = gcc
CFLAGS = -std=c99 -O3
LDLIBS = -lm -lglpk

SRC_C = bslv_main.c bslv_lists.c bslv_vlp.c bslv_lp.c bslv_algs.c bslv_poly.c
OBJS = $(SRC_C:.c=.o)

.PHONY: all
all:	bensolve

bensolve: $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDLIBS)

bslv_main.o: bslv_main.h bslv_vlp.h bslv_lp.h bslv_algs.h
bslv_lists.o: bslv_main.h bslv_lists.h
bslv_vlp.o: bslv_main.h bslv_vlp.h bslv_algs.h bslv_lists.h
bslv_lp.o: bslv_main.h bslv_lists.h bslv_lp.h
bslv_algs.o: bslv_main.h bslv_lists.h bslv_poly.h bslv_vlp.h bslv_lp.h bslv_algs.h
bslv_poly.o: bslv_poly.h

TAGS: $(SRC_C)
	etags $(SRC_C)

.PHONY : clean
clean :
	rm -f bensolve $(OBJS)
