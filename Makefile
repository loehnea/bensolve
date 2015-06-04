OBJS = bslv_algs.o bslv_lists.o bslv_lp.o bslv_main.o bslv_poly.o bslv_vlp.o

all:	bensolve

bensolve:	$(OBJS)
	gcc -o $@ $^ -lm -lglpk

%.o:	%.c %.h
	gcc -std=c99 -c $<
