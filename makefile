CC = g++ -O3 -Wall -o
INC = c-3po.cpp tr.cpp

default: main
all: main test
main:
	$(CC) c3po main.cpp $(INC)
test:
	$(CC) test test.cpp $(INC)
clean:
	rm test c3po

.PHONY: main test
