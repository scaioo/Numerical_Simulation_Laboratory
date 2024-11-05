CC = g++
CFLAGS = -Wall -O3 --std=c++17
AFLAGS = -larmadillo

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS) $(AFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS) $(AFLAGS)
clean :
	rm *.o main.exe seed.out
