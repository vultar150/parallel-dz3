all: compile

compile: main.cpp
	mpicxx -O3  -std=c++11  main.cpp

clean:
	rm -rf a.out
