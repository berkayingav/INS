#pragma once
#ifndef matrixP_H
#define matrixP_H

#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>

template <typename T>

class Matrix {
	private:
		std::vector<std::vector<T>> data;
		size_t rows, cols;
	public:

		Matrix() : rows(0), cols(0) {}

		Matrix(size_t row,size_t col) : rows(row),cols(col),data(row,std::vector<T>(col,0)) {}

		Matrix(const std::vector<std::vector<T>>& input) : data (input){  //initializer list
			rows = input.size();
			cols = (rows > 0) ? input[0].size() : 0;
		}

		size_t getRows() const { return rows; }
		size_t getCols() const { return cols; }


		// (operator overloading)
		std::vector<T>& operator[](size_t i) { return data[i]; }
		const std::vector<T>& operator[](size_t i) const { return data[i]; }

		bool operator==(const Matrix<T>& other) const {
			if (rows != other.rows || cols != other.cols) return false; // Boyutlar farklýysa eþit olamaz

			for (size_t i = 0; i < rows; i++) {
				for (size_t j = 0; j < cols; j++) {
					if (data[i][j] != other.data[i][j]) return false;
				}
			}
			return true;
		}

		Matrix<T> operator*(double scalar) const {
			Matrix<T> result(rows, cols);
			for (size_t i = 0; i < rows; ++i) {
				for (size_t j = 0; j < cols; ++j) {
					result[i][j] = data[i][j] * scalar;
				}
			}
			return result;
		}

		friend Matrix<T> operator*(double scalar, const Matrix<T>& mat) {
			return mat * scalar; // Üye fonksiyonu çaðýrýyoruz
		}

		Matrix<T> operator+(const Matrix<T>& other) const {
			if (rows != other.rows || cols != other.cols) {
				throw std::invalid_argument("Matrix dimensions must match for addition");
			}

			Matrix<T> result(rows, cols);
			for (size_t i = 0; i < rows; ++i) {
				for (size_t j = 0; j < cols; ++j) {
					result[i][j] = data[i][j] + other[i][j];
				}
			}
			return result;
		}

		Matrix<T> operator-(const Matrix<T>& other) const {
			if (rows != other.rows || cols != other.cols) {
				throw std::invalid_argument("Matrix dimensions must match for addition");
			}

			Matrix<T> result(rows, cols);
			for (size_t i = 0; i < rows; ++i) {
				for (size_t j = 0; j < cols; ++j) {
					result[i][j] = data[i][j] - other[i][j];
				}
			}
			return result;
		}

		Matrix<T> operator*(const Matrix<T>& other) const {
			if (cols != other.rows) {
				throw std::invalid_argument("Matrix dimensions not suitable for multiplication!");
			}

			Matrix<T> result(rows, other.cols);

			for (size_t i = 0; i < rows; ++i) {
				for (size_t j = 0; j < other.cols; ++j) {
					result[i][j] = 0;  // Baþlangýç deðeri
					for (size_t k = 0; k < cols; ++k) {
						result[i][j] += data[i][k] * other[k][j];
					}
				}
			}
			return result;
		}


		static Matrix<T> add(const Matrix<T>& A, const Matrix<T>& B) {
			if (A.rows != B.rows || A.cols != B.cols)
				throw std::invalid_argument(" Matrix Sizes dont match! Check rows and cols.");

			Matrix<T> result(A.rows, A.cols);
			for (size_t i = 0; i < A.rows; ++i)
				for (size_t j = 0; j < A.cols; ++j)
					result[i][j] = A[i][j] + B[i][j];

			return result;
		}

		static Matrix<T> subtract(const Matrix<T>& A, const Matrix<T>& B) {
			if (A.rows != B.rows || A.cols != B.cols)
				throw std::invalid_argument(" Matrix dimension dont same! Check rows and cols.");

			Matrix<T> result(A.rows, A.cols);
			for (size_t i = 0; i < A.rows; ++i)
				for (size_t j = 0; j < A.cols; ++j)
					result[i][j] = A[i][j] - B[i][j];

			return result;
		}

		static Matrix<T> multiply(const Matrix<T>& A, const Matrix<T>& B) {
			if (A.cols != B.rows) {
				throw std::invalid_argument("Matrix dimension not suitable for multiplication!");
			}
			
			Matrix<T> result(A.rows, B.cols);
			for (size_t i = 0; i < A.rows; i++) {
				for (size_t j = 0; j < B.cols; j++) {
					result[i][j] = 0;  // Ýlk deðerleri sýfýra ayarla

					for (size_t k = 0; k < A.cols; k++) {
						result[i][j] += A[i][k] * B[k][j];
					}
				}
			}
			return result;
		}

		static Matrix<T> transpose(const Matrix<T>& A) {
			Matrix<T> result(A.cols, A.rows);
			for (size_t i = 0; i < A.cols; i++)
				for (size_t j = 0; j < A.rows; j++)
					result[i][j] = A[j][i];

			return result;
		}

		static Matrix<T> SkewSymetric(const Matrix<T>& A) {
			if (!((A.getRows() == 3 && A.getCols() == 1) || (A.getRows() == 1 && A.getCols() == 3))) {
				throw std::invalid_argument("Input must be 3x1 or 1x3 matrix.");
			}

			Matrix<T> result(3, 3);
			if (A.getRows() == 3) {
				result[0][0] = 0;
				result[0][1] = -A[2][0];
				result[0][2] = A[1][0];

				result[1][0] = A[2][0];
				result[1][1] = 0;
				result[1][2] = -A[0][0];

				result[2][0] = -A[1][0];
				result[2][1] = A[0][0];
				result[2][2] = 0;
			}

			else {
				result[0][0] = 0;
				result[0][1] = -A[0][2];
				result[0][2] = A[0][1];

				result[1][0] = A[0][2];
				result[1][1] = 0;
				result[1][2] = -A[0][0];

				result[2][0] = -A[0][1];
				result[2][1] = A[0][0];
				result[2][2] = 0;
			}


			return result;
		}


};








#endif