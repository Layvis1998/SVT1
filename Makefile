prog: main.cpp
	g++  main.cpp lib/libblas.so lib/libumfpack.a lib/libamd.a -O0 -g -o prog
clean:
	rm prog
