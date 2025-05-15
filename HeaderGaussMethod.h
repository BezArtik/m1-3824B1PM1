#pragma once

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <cstdlib>
#include <random>
#include <iomanip>
#include <limits>

template <typename T>
class Vector {
protected:
    T* data{ nullptr };
    size_t size{ 0 };
public:
    Vector() = default;
    Vector(size_t n) {
        try {
            data = new T[n]();
            size = n;
        }
        catch (...) {
            throw std::runtime_error("Failed to allocate memory.");
        }
    }
    ~Vector() { delete[] data; }

    Vector(const Vector& other) : data{ nullptr }, size{ 0 } {
        T* new_data = new T[other.size];
        try {
            for (size_t i = 0; i < other.size; ++i) {
                new_data[i] = other.data[i];
            }
            data = new_data;
            size = other.size;
        }
        catch (...) {
            delete[] new_data;
            throw std::runtime_error("Vector copy failed");
        }
        
    }

    size_t GetSize() const noexcept {
        return size;
    }

    Vector& operator=(const Vector& other) {
        if (this == &other) {
            return *this;
        }
        Vector temp(other);
        this->Swap(temp);
        return *this;
    }

    T& operator[](size_t i) {
        if (i >= size) {
            throw std::out_of_range("Vector index out of range");
        }
        return data[i];
    }

    const T& operator[](size_t i) const {
        if (i >= size) {
            throw std::out_of_range("Vector index out of range");
        }
        return data[i];
    }

    void Swap(Vector& other) noexcept {
        std::swap(data, other.data);
        std::swap(size, other.size);
    }
};

template <typename T>
class Matrix : public Vector<Vector<T>> {
    size_t rows{ 0 };
    size_t cols{ 0 };
public:
    Matrix() = default;
    Matrix(size_t r, size_t c) : Vector<Vector<T>>(r), rows(r), cols(c) {
        try {
            for (size_t i = 0; i < rows; ++i) {
                new (&this->data[i]) Vector<T>(cols);
            }
        }
        catch (...) {
            if (this->data) {
                for (size_t i = 0; i < rows; ++i) {
                    this->data[i].~Vector<T>();
                }
            }
            if (this->data) {
                delete[] this->data;
                this->data = nullptr;
            }
            this->size = 0;
            throw std::runtime_error("Failed to construct Matrix");
        }
    }
    ~Matrix() = default;
    Matrix(const Matrix& other) : Vector<Vector<T>>(other), rows(other.rows), cols(other.cols) {}

    Matrix& operator=(const Matrix& other) {
        if (this == &other) {
            return *this;
        }
        Matrix temp(other);
        this->Swap(temp);
        return *this;
    }

    size_t GetRows() const noexcept {
        return rows;
    }
    size_t GetCols() const noexcept {
        return cols;
    }

    Vector<T>& operator[](size_t row) {
        if (row >= rows) {
            throw std::out_of_range("Matrix row index out of range");
        }
        return Vector<Vector<T>>::operator[](row);
    }

    const Vector<T>& operator[](size_t row) const {
        if (row >= rows) {
            throw std::out_of_range("Matrix row index out of range");
        }
        return Vector<Vector<T>>::operator[](row);
    }

    void Swap(Matrix& other) noexcept {
        Vector<Vector<T>>::Swap(other);
        std::swap(rows, other.rows);
        std::swap(cols, other.cols);
    }

    void SwapRows(size_t i, size_t j) {
        if (i >= rows || j >= rows) {
            throw std::out_of_range("Row indices out of range");
        }
        if (i != j) {
            (*this)[i].Swap((*this)[j]);
        }
    }
};

template <typename T>
class Solution {
public:
    bool is_consistent{ false };
    bool has_unique{ false };
    Vector<T> particular;
    Vector<T>* basis{ nullptr };
    size_t basis_size{ 0 };
    size_t* free_vars{ nullptr };
    size_t rank{ 0 };

    Solution() = default;
    ~Solution() {
        delete[] basis;
        delete[] free_vars;
    }
};

template <typename T>
Solution<T> GaussMethod(const Matrix<T>& A_, const Vector<T>& b_) {
    Matrix<T> A = A_;
    Vector<T> b = b_;
    const size_t m = A.GetRows();
    const size_t n = A.GetCols();
    const T eps = std::numeric_limits<T>::epsilon() * 1000;
    Solution<T> sol;
    sol.is_consistent = true;
    size_t rank = 0;
    size_t* pivot_cols = new size_t[m];

    for (size_t i = 0; i < m; ++i) {
        pivot_cols[i] = n;
    }

    for (size_t col = 0; col < n && rank < m; ++col) {
        size_t max_row = rank;
        for (size_t i = rank; i < m; ++i) {
            if (std::abs(A[i][col]) > std::abs(A[max_row][col])) {
                max_row = i;
            }
        }

        if (std::abs(A[max_row][col]) < eps) continue;

        pivot_cols[rank] = col;

        if (max_row != rank) {
            A.SwapRows(rank, max_row);
            std::swap(b[rank], b[max_row]);
        }

        T pivot = A[rank][col];
        for (size_t j = col; j < n; ++j) {
            A[rank][j] /= pivot;
        }
        b[rank] /= pivot;

        for (size_t i = 0; i < m; ++i) {
            if (i != rank && std::abs(A[i][col]) > eps) {
                T factor = A[i][col];
                for (size_t j = col; j < n; ++j) {
                    A[i][j] -= factor * A[rank][j];
                }
                b[i] -= factor * b[rank];
            }
        }
        ++rank;
    }

    for (size_t i = rank; i < m; ++i) {
        if (std::abs(b[i]) > eps) {
            sol.is_consistent = false;
            delete[] pivot_cols;
            return sol;
        }
    }

    sol.rank = rank;
    sol.has_unique = (rank == n);
    sol.particular = Vector<T>(n);

    bool* is_pivot = new bool[n]();
    for (size_t i = 0; i < rank; ++i) {
        if (pivot_cols[i] < n) {
            is_pivot[pivot_cols[i]] = true;
        }
    }

    for (ptrdiff_t i = rank - 1; i >= 0; --i) {
        size_t col = pivot_cols[i];
        if (col >= n) continue;

        sol.particular[col] = b[i];
        for (size_t j = col + 1; j < n; ++j) {
            if (!is_pivot[j]) {
                sol.particular[col] -= A[i][j] * sol.particular[j];
            }
        }
    }

    sol.basis_size = 0;
    for (size_t j = 0; j < n; ++j) {
        if (!is_pivot[j]) ++sol.basis_size;
    }

    if (sol.basis_size > 0) {
        sol.basis = new Vector<T>[sol.basis_size];
        sol.free_vars = new size_t[sol.basis_size];
        size_t basis_idx = 0;

        for (size_t free_col = 0; free_col < n; ++free_col) {
            if (!is_pivot[free_col]) {
                sol.free_vars[basis_idx] = free_col;
                sol.basis[basis_idx] = Vector<T>(n);
                sol.basis[basis_idx][free_col] = 1;

                for (ptrdiff_t i = rank - 1; i >= 0; --i) {
                    size_t pivot_col = pivot_cols[i];
                    if (pivot_col >= n) continue;

                    sol.basis[basis_idx][pivot_col] = 0;
                    for (size_t j = pivot_col; j < n; ++j) {
                        sol.basis[basis_idx][pivot_col] -= A[i][j] * sol.basis[basis_idx][j];
                    }
                }
                ++basis_idx;
            }
        }
    }

    delete[] is_pivot;
    delete[] pivot_cols;
    return sol;
}

template <typename T>
void CheckSolution(const Matrix<T>& A, const Vector<T>& b, const Solution<T>& solution) {
    const T eps = std::numeric_limits<T>::epsilon() * static_cast<T>(A.GetRows()) * 1000;
    if (!solution.is_consistent) {
        std::cout << "Solution is correct!\n";
        return;
    }
    for (size_t i = 0; i < A.GetRows(); ++i) {
        T sum = 0;
        for (size_t j = 0; j < A.GetCols(); ++j) {
            sum += A[i][j] * solution.particular[j];
        }
        T error = std::abs(sum - b[i]);
        T max_val = std::max(T(1), std::max(std::abs(sum), std::abs(b[i])));
        if (error > eps * max_val) {
            throw std::runtime_error("Solution is incorrect.\n");
        }
    }
    std::cout << "Solution is correct!\n";
}

template <typename T>
void PrintSolution(const Solution<T>& sol) {
    const T eps = std::numeric_limits<T>::epsilon() * 100;

    if (!sol.is_consistent) {
        std::cout << "System is inconsistent (no solutions)\n";
        return;
    }

    std::cout << "Solution:\n";

    if (sol.has_unique) {
        for (size_t i = 0; i < sol.particular.GetSize(); ++i) {
            std::cout << "x[" << i + 1 << "] = " << sol.particular[i] << "\n";
        }
    }
    else {
        for (size_t i = 0; i < sol.particular.GetSize(); ++i) {
            std::cout << "x[" << i + 1 << "] = " << sol.particular[i];

            for (size_t k = 0; k < sol.basis_size; ++k) {
                if (std::abs(sol.basis[k][i]) > eps) {
                    std::cout << " + " << sol.basis[k][i] << "*t" << k + 1;
                }
            }
            std::cout << "\n";
        }

        for (size_t k = 0; k < sol.basis_size; ++k) {
            std::cout << "t" << k + 1 << " is free (x[" << sol.free_vars[k] + 1 << "])\n";
        }
    }
}

enum Choice {
    RandomMatrix = 1,
    ManualInput = 2,
    Exit = 3
};

template <typename T>
void ProcessRandomMatrix() {
    size_t n;
    std::cout << "Enter matrix size (n x n): ";
    std::cin >> n;

    if (n == 0) {
        throw std::invalid_argument("Matrix size must be positive");
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    T max_val = static_cast<T>(1000.0 / std::sqrt(static_cast<double>(n)));
    std::uniform_real_distribution<T> dist(-max_val, max_val);

    Matrix<T> A(n, n);
    Vector<T> b(n);

    for (size_t i = 0; i < n; ++i) {
        T row_sum = 0;
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                A[i][j] = dist(gen);
                row_sum += std::abs(A[i][j]);
            }
        }
        A[i][i] = row_sum * (1 + dist(gen) / static_cast<T>(1000.0));
        b[i] = dist(gen);
    }

    const auto start{ std::chrono::steady_clock::now() };

    Solution<T> sol = GaussMethod(A, b);

    const auto finish{ std::chrono::steady_clock::now() };
    const std::chrono::duration<double> elapsed_seconds{ finish - start };

    CheckSolution(A, b, sol);
    std::cout << "Time of work: ";
    std::cout << elapsed_seconds.count() << "s\n";
}

template <typename T>
void ProcessManualInput() {
    size_t m, n;
    std::cout << "Enter number of equations (m): ";
    std::cin >> m;
    std::cout << "Enter number of variables (n): ";
    std::cin >> n;

    if (m == 0 || n == 0) {
        throw std::invalid_argument("Matrix size must be positive");
    }

    Matrix<T> A(m, n);
    Vector<T> b(m);

    std::cout << "Enter augmented matrix of coefficients:\n";
    for (size_t i = 0; i < m; ++i) {
        std::cout << "Row " << i + 1 << ": ";
        for (size_t j = 0; j < n; ++j) {
            if (!(std::cin >> A[i][j])) {
                throw std::invalid_argument("Invalid input");
            }
        }
        if (!(std::cin >> b[i])) {
            throw std::invalid_argument("Invalid input");
        }
    }
    Solution<T> sol = GaussMethod(A, b);

    CheckSolution(A, b, sol);
    PrintSolution(sol);
}

template <typename T>
void Run() {
    
    std::cout << "1. Generate random matrix\n";
    std::cout << "2. Enter matrix manually\n";
    std::cout << "3. Exit\n";
    std::cout << "Your choice: ";

    size_t choice;
    std::cin >> choice;

    switch (choice) {
    case RandomMatrix:
        ProcessRandomMatrix<T>();
        break;
    case ManualInput:
        ProcessManualInput<T>();
        break;
    case Exit:
        std::cout << "Exiting program...\n";
        return;
    default:
        throw std::invalid_argument("Invalid choice");
    }
    std::cout << std::endl;
}