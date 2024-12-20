all : dstgs

.PHONY : all
CXX = g++
CFLAGS = -lpthread -static-libstdc++ -std=c++11


dstgs : dstgs.o
	$(CXX) -o dstgs dstgs.o $(CFLAGS)

dstgs.o : main.cpp QueryFunction.h DSTGS.h HashFunction.h params.h
	$(CXX) -o dstgs.o -c main.cpp -D HINT
	

.PHONY:clean
clean:
	-$(RM) dstgs
	-$(RM) dstgs.o