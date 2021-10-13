prog: main.cpp
	g++  main.cpp lib/libblas.so lib/libumfpack.a lib/libamd.a -O3 -o prog
clean:
	rm prog
