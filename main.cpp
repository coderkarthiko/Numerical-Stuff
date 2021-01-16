#include "numcpp.h"

// test program
int main() {
	mat M = get_mat();
	std::cout << std::fixed << "det - " << det(M) << std::endl;
	print_mat(mul(M, inv(M)));
}
