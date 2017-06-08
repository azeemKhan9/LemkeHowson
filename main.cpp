#include "stdafx.h"
#include <iostream>
#include <Eigen\Dense>

using namespace Eigen;

typedef Eigen::Matrix<float, 3, 4> Matrix34f;

int main()
{
	Matrix34f A, B;
	A << 4, 12, 8, 6,
		16, 8, 12, 8,
		10, 8, 10, 9;
	B << 25, 5, 5, 8,
		1, 15, 8, 4,
		17, 10, 13, 9;
	std::cout << A << std::endl;
	std::cout << B << std::endl;

	return 0;
}