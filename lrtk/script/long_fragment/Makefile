default: construct_fragment
construct_fragment : fragment.o
	g++ -O3 -o construct_fragment fragment.o -pthread -lhts
fragment.o : fragment.cpp cmdline.h ThreadPool/ThreadPool.h
	g++ -O3 -c fragment.cpp
clean :
	rm -f fragment.o construct_fragment
