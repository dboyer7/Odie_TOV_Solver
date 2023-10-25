CC = gcc
CFLAGS = -O2 -march=native -fopenmp -Wall -Wno-unused-variable -std=gnu99

INCLUDEDIRS=
all: PolyTOV TOV_Driver_main.o
TOV_Driver_main.o: TOV_Driver_main.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c TOV_Driver_main.c -o TOV_Driver_main.o

PolyTOV:	TOV_Driver_main.o
	$(CC)	TOV_Driver_main.o -o PolyTOV -lm -lgomp

clean:
	rm -f *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png PolyTOV
