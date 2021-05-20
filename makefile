CC = g++ -O3 -Wall -o
INC = src/c-3po.cpp src/tr.cpp

default: main
all: main test
main:
	$(CC) c3po src/main.cpp $(INC)
test:
	$(CC) test src/test.cpp $(INC)
clean:
	rm test c3po

.PHONY: main test
