all: clean test_matrix
	./test_matrix

clean:
	rm -f ./test_matrix

test_matrix: test_matrix.cpp matrix.h
	g++ --std=c++11 test_matrix.cpp -Wall -o test_matrix
