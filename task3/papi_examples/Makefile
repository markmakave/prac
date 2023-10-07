CC = gcc
LIBS = -lpapi
CFLAGS = -Wall -Werror

PROGS = example1-array-vs-list example2-array-sum-flops matrix-mul

all: $(PROGS)

% : %.c
	$(CC) $(CFLAGS) $< $(LIBS) -o $@

clean:
	rm -f $(PROGS)
