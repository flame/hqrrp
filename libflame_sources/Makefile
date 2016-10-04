
# Defintions of variables.

LIBFLAME_PATH = /usr/local/libflame-current-gc

CC      = gcc
CCFLAGS = -O -fopenmp -I$(LIBFLAME_PATH)/include-x86_64-r -msse3
LD      = gcc
LDFLAGS = -O -fopenmp -I$(LIBFLAME_PATH)/include-x86_64-r -msse3

# Defintions of rules.

simple_test.x : simple_test.o FLA_HQRRP_UT_blk_var2.o
	$(LD) $(LDFLAGS) \
            -o simple_test.x \
            simple_test.o \
            FLA_HQRRP_UT_blk_var2.o \
            -L$(LIBFLAME_PATH)/lib -lflame-x86_64-r  \
            /usr/local/lapack/mt_openblas/lib/libopenblas_haswellp-r0.2.14.a \
            -lm

simple_test.o : simple_test.c
	$(CC) $(CCFLAGS) -c simple_test.c

%.o : %.c
	$(CC) $(CCFLAGS) -c $< -o $@

clean: 
	rm -f a.out *.x *.o *~ core

