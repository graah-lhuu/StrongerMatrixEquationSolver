#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

const double EPS = 1e-10;

bool isZero(double val) {
	return fabs(val) < EPS;
}

void printMatrix(const vector<vector<double>>& matrix, const string& name = "") {
	if (!name.empty()) {
		cout << name << ":" << endl;
	}
	for (const auto& row : matrix) {
		for (double val : row) {
			if (isZero(val)) {
				cout << setw(12) << "0";
			} else {
				cout << setw(12) << fixed << setprecision(6) << val;
			}
		}
		cout << endl;
	}
	cout << endl;
}

void printVector(const vector<double>& vec, const string& name = "") {
	if (!name.empty()) {
		cout << name << ": ";
	}
	cout << "[ ";
	for (double val : vec) {
		if (isZero(val)) {
			cout << "0 ";
		} else {
			cout << fixed << setprecision(6) << val << " ";
		}
	}
	cout << "]" << endl;
}

int matrixRank(const vector<vector<double>>& matrix) {
	int rows = matrix.size();
	if (rows == 0) return 0;
	int cols = matrix[0].size();
	int rank = 0;

	vector<vector<double>> temp = matrix;

	for (int col = 0, row = 0; col < cols && row < rows; col++) {
		int pivotRow = row;
		for (int i = row + 1; i < rows; i++) {
			if (fabs(temp[i][col]) > fabs(temp[pivotRow][col])) {
				pivotRow = i;
			}
		}

		if (isZero(temp[pivotRow][col])) {
			continue;
		}

		if (pivotRow != row) {
			swap(temp[row], temp[pivotRow]);
		}

		for (int i = row + 1; i < rows; i++) {
			double factor = temp[i][col] / temp[row][col];
			for (int j = col; j < cols; j++) {
				temp[i][j] -= factor * temp[row][j];
			}
		}

		row++;
		rank++;
	}

	return rank;
}

vector<vector<double>> matrixMultiply(const vector<vector<double>>& A, const vector<vector<double>>& B) {
	int m = A.size();
	int n = A[0].size();
	int p = B[0].size();

	vector<vector<double>> result(m, vector<double>(p, 0.0));

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) {
			for (int k = 0; k < n; k++) {
				result[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return result;
}

vector<vector<double>> matrixTranspose(const vector<vector<double>>& matrix) {
	int rows = matrix.size();
	int cols = matrix[0].size();

	vector<vector<double>> result(cols, vector<double>(rows));

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			result[j][i] = matrix[i][j];
		}
	}

	return result;
}

vector<vector<double>> createIdentityMatrix(int n) {
	vector<vector<double>> identity(n, vector<double>(n, 0.0));
	for (int i = 0; i < n; i++) {
		identity[i][i] = 1.0;
	}
	return identity;
}

vector<double> solveLinearSystem(vector<vector<double>> A, vector<double> b) {
	int n = A.size();

	for (int i = 0; i < n; i++) {
		int pivotRow = i;
		for (int j = i + 1; j < n; j++) {
			if (fabs(A[j][i]) > fabs(A[pivotRow][i])) {
				pivotRow = j;
			}
		}

		if (pivotRow != i) {
			swap(A[i], A[pivotRow]);
			swap(b[i], b[pivotRow]);
		}

		if (isZero(A[i][i])) {
			continue;
		}

		double pivot = A[i][i];
		for (int j = i; j < n; j++) {
			A[i][j] /= pivot;
		}
		b[i] /= pivot;

		for (int j = 0; j < n; j++) {
			if (j != i && !isZero(A[j][i])) {
				double factor = A[j][i];
				for (int k = i; k < n; k++) {
					A[j][k] -= factor * A[i][k];
				}
				b[j] -= factor * b[i];
			}
		}
	}

	return b;
}

void analyzeAndSolve(const vector<vector<double>>& A, const vector<double>& b) {
	int m = A.size();
	if (m == 0) {
		cout << "系数矩阵A为空！" << endl;
		return;
	}
	int n = A[0].size();

	cout << "原始线性方程组 A * x = b:" << endl;
	cout << "A:" << endl;
	printMatrix(A);
	cout << "b: ";
	printVector(b);

	int rankA = matrixRank(A);
	cout << "系数矩阵 A 的秩: " << rankA << endl;

	vector<vector<double>> augmented = A;
	for (int i = 0; i < m; i++) {
		augmented[i].push_back(b[i]);
	}

	int rankAB = matrixRank(augmented);
	cout << "增广矩阵 [A|b] 的秩: " << rankAB << endl << endl;

	if (rankA == rankAB) {
		if (rankA == n) {
			cout << "方程组有唯一解。" << endl;
		} else {
			cout << "方程组有无穷多解。" << endl;
		}

		cout << "由于方程组有解，不需要计算最小二乘解。" << endl;
	} else {
		cout << "方程组无解。" << endl;

		if (rankA == n) {
			cout << "系数矩阵 A 列满秩，可以计算最小二乘解。" << endl;

			vector<vector<double>> A_T = matrixTranspose(A);
			vector<vector<double>> ATA = matrixMultiply(A_T, A);
			vector<double> ATb = matrixMultiply(A_T, { b })[0];

			cout << "正规方程: (A^T * A) * x_hat = A^T * b" << endl;
			cout << "A^T * A:" << endl;
			printMatrix(ATA, "ATA");
			cout << "A^T * b: ";
			printVector(ATb, "ATb");

			vector<double> x_hat = solveLinearSystem(ATA, ATb);

			cout << "最小二乘解 x_hat: ";
			printVector(x_hat);

			vector<double> b_hat(m);
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					b_hat[i] += A[i][j] * x_hat[j];
				}
			}

			cout << "投影后的 b_hat = A * x_hat: ";
			printVector(b_hat);

			double error = 0.0;
			for (int i = 0; i < m; i++) {
				error += (b[i] - b_hat[i]) * (b[i] - b_hat[i]);
			}
			cout << "残差平方和: " << error << endl;

		} else {
			cout << "系数矩阵 A 非列满秩，无法求得唯一的最小二乘解。" << endl;
			cout << "可以考虑使用奇异值分解(SVD)等方法求广义逆解。" << endl;
		}
	}
}

void interactiveInput() {
	int m, n;

	cout << "================================" << endl;
	cout << "线性方程组求解器（支持最小二乘解）" << endl;
	cout << "求解 A * x = b" << endl;
	cout << "================================" << endl;

	cout << "请输入方程个数 m (矩阵A的行数): ";
	cin >> m;
	cout << "请输入未知数个数 n (矩阵A的列数): ";
	cin >> n;

	vector<vector<double>> A(m, vector<double>(n));
	cout << "\n请输入系数矩阵 A (" << m << "行 × " << n << "列):" << endl;
	for (int i = 0; i < m; i++) {
		cout << "输入第 " << i + 1 << " 行的 " << n << " 个元素: ";
		for (int j = 0; j < n; j++) {
			cin >> A[i][j];
		}
	}

	vector<double> b(m);
	cout << "\n请输入常数向量 b (" << m << "个元素):" << endl;
	for (int i = 0; i < m; i++) {
		cout << "输入 b[" << i + 1 << "]: ";
		cin >> b[i];
	}

	cout << "\n" << string(50, '=') << endl;
	cout << "求解结果:" << endl;
	cout << string(50, '=') << endl;

	analyzeAndSolve(A, b);
}

int main() {
	interactiveInput();

	cout << "\n程序运行结束。" << endl;
	return 0;
}