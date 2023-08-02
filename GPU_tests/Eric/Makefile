CC = gcc
OMPFLAGS = -foffload=nvptx-none -foffload=-misa=sm_75 -foffload=-lm -fcf-protection=none -fno-stack-protector
eric-gpu-test: eric-gpu-test.c
	$(CC) $(CFLAGS) -o eric-gpu-test -fopenmp $(OMPFLAGS) eric-gpu-test.c -lm

.PHONY: clean

clean:
	-rm parallel_region_test parallel_performance.txt
	-rm *.o
