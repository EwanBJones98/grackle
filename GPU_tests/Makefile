CC = gcc

GPU_OMPFLAGS = -foffload=nvptx-none -foffload=-misa=sm_75 -foffload=-lm -fcf-protection=none -fno-stack-protector

parallel_region_cpu: parallel_region_test.c
	$(CC) -D OMP_CPU -o parallel_region_test_cpu -fopenmp parallel_region_test.c

parallel_region_gpu: parallel_region_test.c
	$(CC) -D OMP_GPU -o parallel_region_test_gpu -fopenmp $(GPU_OMPFLAGS) -no-pie parallel_region_test.c -lm

parallel_region_serial: parallel_region_test.c
	$(CC) -o parallel_region_test_serial parallel_region_test.c

.PHONY: clean

clean:
	-rm parallel_region_test parallel_performance.txt
	-rm *.o
