#include "numcpp.h"

int main() {
	mat M = get_mat();
	std::cout << std::fixed << "det - " << det(M) << std::endl;
	print_mat(inv(M));
}