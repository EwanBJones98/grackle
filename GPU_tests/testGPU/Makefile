CC = gcc

#GPU_OMPFLAGS = -foffload=nvptx-none -foffload=-misa=sm_30 -foffload=-lm -fcf-protection=none -fno-stack-protector
GPU_OMPFLAGS = -foffload=nvptx-none -fcf-protection=none -no-pie

testGPU: testGPU.c
	$(CC) -D GPU -o testGPU -fopenmp $(GPU_OMPFLAGS) testGPU.c -lm

testCPU: testGPU.c
	$(CC) -D CPU -o testGPU -fopenmp testGPU.c

testSerial: testGPU.c
	$(CC) -D SERIAL -o testGPU testGPU.c

.PHONY: clean

clean:
	-rm testGPU
	-rm *.o
	-rm output.log
	-rm slurm-*