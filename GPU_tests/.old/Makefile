CC = gcc
CFLAGS =
OMPFLAGS = 

ifeq ($(OMP_OPT),cpu)
	OMPFLAGS += -fopenmp
	CFLAGS += -D OMP_CPU
else ifeq ($(OMP_OPT),gpu-nvidia)
	OMPFLAGS += -fopenmp -foffload=nvptx-none \
	                     -foffload=-misa=sm_75 \
			     -foffload=-lm \
			     -fcf-protection=none -fno-stack-protector \
			     -no-pie 
	CFLAGS += -D OMP_GPU
else ifeq ($(OMP_OPT),gpu-amd)
	OMPFLAGS += -fopenmp
	CFLAGS += -D OMP_GPU
endif

parallel_region: parallel_region_test.c
	$(CC) $(CFLAGS) -o parallel_region_test $(OMPFLAGS) parallel_region_test.c

eric-gpu-test: eric-gpu-test.c
	$(CC) $(CFLAGS) -o eric-gpu-test $(OMPFLAGS) eric-gpu-test.c -lm

.PHONY: clean

clean:
	-rm parallel_region_test parallel_performance.txt test-gpu-openmp
	-rm *.o
