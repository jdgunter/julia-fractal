all:    compile
	$(MAKE) link
seq:	compile
	$(MAKE) seq_link
compile:
	g++ -g -std=c++11 -c julia.cpp `libpng-config --cflags`
link:
	g++ -o julia julia.o -fopenmp `libpng-config --ldflags`
seq_link:
	g++ -o julia julia.o `libpng-config --ldflags`
