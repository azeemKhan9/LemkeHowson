#include "stdafx.h"
#include <iostream>
#include <Eigen\Dense>
#include <vector>
#include <time.h>

using namespace Eigen;

MatrixXf pivot(MatrixXf M, int row, int column);

void isDegenerate(VectorXf x, VectorXf y);

typedef Eigen::Matrix<float, 4, 4> MatrixXXf; // Will have to be adjusted depending on input matrix sizes

int main()
{
	int k0 = 1; //Initial pivot label
	int maxPivots = 10000;
	int player;
	int numPivots = 0;
	MatrixXf A, B;
	MatrixXf LP;
	//MatrixXXf A, B;
	
	// Matrix initialisation - adjust depending on matrix size
	std::srand((unsigned int)time(0));
	A = 10 * MatrixXf::Random(4, 4);
	B = 10 * MatrixXf::Random(4, 4);
	for (int i = 0; i < A.rows(); i++) {
		for (int j = 0; j < A.cols(); j++) {
			A(i, j) = round(abs(A(i, j)));
			//B(i, j) = round(abs(B(i, j)));
			if (i == j) {
				B(i, j) = 1.0;
			}
			else {
				B(i, j) = 0.0;
			}
		}
	}

	/*A << 4, 4, 0, 2,
		8, 9, 3, 1,
		5, 6, 2, 5,
		8, 1, 7, 10;
	B << 9, 9, 10, 8,
		8, 3, 2, 5,
		1, 10, 7, 2,
		1, 1, 1, 9;*/
	std::cout << A << std::endl;
	std::cout << B << std::endl;
	int m = A.rows();
	int n = A.cols();

	int k = k0;
	if (k0 <= m) {
		player = 0;
	}
	else if (k0 < 1 || k0 > m + n) {
		std::cerr << "ERROR : Initial pivot is outside range." << std::endl;
		return 1;
	}
	else {
		player = 1;
	}

	//Construct tableaus
	MatrixXf tab1(n, m + n + 1), tab2(m, m + n + 1);
	tab1.block(0, 0, n, m) = B.transpose();
	tab1.block(0, m, n, n) = MatrixXf::Identity(n, n);
	for (int i = 0; i < n; ++i) {
		tab1(i, m + n) = 1;
	}
	tab2.block(0, 0, m, m) = MatrixXf::Identity(m, m);
	tab2.block(0, m, m, n) = A;
	for (int i = 0; i < m; ++i) {
		tab2(i, m + n) = 1;
	}
	//std::cout << tab1 << std::endl;
	//std::cout << tab2 << std::endl;

	std::vector<int> tab1Labels(m), tab2Labels(n);
	for (int i = 0; i < m; i++) {
		tab1Labels[i] = i + 1;
		//std::cout << tab1Labels[i];
	}
	for (int i = 0; i < n; i++) {
		tab2Labels[i] = m + i + 1;
		//std::cout << tab2Labels[i];
	}
	MatrixXf tables[2] = { tab1, tab2 };
	std::vector<int> labels[2] = { tab2Labels, tab1Labels };

	while (numPivots < maxPivots) { //Pivoting loop
		LP = tables[player];
		int r = LP.rows();
		int c = LP.cols();
		int max = 0;
		int ind = -1;
		float t;
		for (int i = 0; i < r; i++) {
			t = LP(i, k - 1) / LP(i, m + n);
			if (t >= max) {
				ind = i;
				max = t;
			}
		}
		tables[player] = pivot(LP, ind, k - 1);
		//std::cout << tables[player] << std::endl;
		numPivots += 1;
		int temp = labels[player][ind];
		labels[player][ind] = k;
		k = temp;
		if (k == k0) {
			std::cout << "A Nash equilibrium has been found in " << numPivots << " pivots." << std::endl;
			break;
		}
		if (player == 0) {
			player = 1;
		}
		else {
			player = 0;
		}
	}

	if (numPivots >= maxPivots) {
		std::cerr << "ERROR: Maximum number of allowed pivots has been reached!" << std::endl;
		return 1;
	}

	VectorXf nashEqbm[2];
	VectorXf x(3);
	for (int p = 0; p < 2; p++) { //p refers to players
		if (p == 0) {
			x.resize(m);
		}
		else if (p == 1) {
			x.resize(n);
		}
		std::vector<int> rows = labels[p];
		float sum = 0.0;

		/*for (int i = 0; i < rows.size(); i++) {
			std::cout << rows[i] << std::endl;
		}*/

		LP = tables[p];
		//std::cout << tables[p] << std::endl;
		for (int i = 0; i < rows.size(); i++) {
			if (p == 0 && rows[i] <= m) {
				x(rows[i] - 1) = LP(i, m + n) / LP(i, rows[i] - 1);
			}
			else if (p == 1 && rows[i] > m) {
				x(rows[i] - 1 - m) = LP(i, m + n) / LP(i, rows[i] - 1);
			}
		}
		for (int i = 0; i < x.size(); i++) {
			if (x(i) <= 1.0e-6) {
				x(i) = 0;
			}
			else {
				sum = sum + x(i);
			}
		}
		x = (1.0 / sum) * x; //Normalise vector to make it stochastic.
		nashEqbm[p] = x;
		//std::cout << nashEqbm[p] << std::endl;
	}
	isDegenerate(nashEqbm[0], nashEqbm[1]);

	return 0;
}


MatrixXf pivot(MatrixXf M, int row, int column) { //Function that pivots tableau on given element (row, column)
	MatrixXf X = M;
	int m = M.rows();
	int n = M.cols();
	for (int i = 0; i < M.rows(); i++) {
		if (i == row) {
			continue;
		}
		else {
			X.row(i) = X.row(i) - (M(i, column) / M(row, column) * M.row(row));
		}
	}
	X.row(row) = X.row(row) / X(row, column);
	return X;
}

void isDegenerate(VectorXf x, VectorXf y) {
	int elem = 0;
	int elem2 = 0;
	for (int i = 0; i < x.size(); i++) {
		if (x(i) == 0) {
			continue;
		}
		else {
			elem += 1;
		}
	}
	for (int i = 0; i < y.size(); i++) {
		if (y(i) == 0) {
			continue;
		}
		else {
			elem2 += 1;
		}
	}
	if (elem != elem2) {
		std::cout << "WARNING: This game is degenerate" << std::endl;
	}
}