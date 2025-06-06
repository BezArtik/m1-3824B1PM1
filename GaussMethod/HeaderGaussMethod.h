#pragma once

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <random>
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

    Vector(const Vector& other): data(new T[other.size]), size(other.size) {
        try {
            for (size_t i = 0; i < size; ++i) {
                data[i] = other.data[i];  
            }
        }
        catch (...) {
            delete[] data;  
            throw;  
        }
    }

    size_t GetSize() const noexcept { return size; }

    void Swap(Vector& other) noexcept; 

    Vector& operator=(const Vector& other);

    T& operator[](size_t i); 

    const T& operator[](size_t i) const; 

    Vector& operator/=(T scalar); 

    Vector& operator-=(const Vector& rhs); 
    
    Vector& operator+=(const Vector& rhs); 

    Vector operator*(T scalar) const;

private:
    void CheckSize(const Vector& other) const;
};

template <typename T>
void Vector<T>::CheckSize(const Vector& other) const {
    if (size != other.size) {
        throw std::invalid_argument("Vector sizes don't match");
    }
}

template <typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& other) {
    if (this == &other) {
        return *this;
    }
    Vector temp(other);
    this->Swap(temp);
    return *this;
}

template <typename T>
T& Vector<T>::operator[](size_t i) {
    if (i >= size) {
        throw std::out_of_range("Vector index out of range");
    }
    return data[i];
}

template <typename T>
const T& Vector<T>::operator[](size_t i) const {
    if (i >= size) {
        throw std::out_of_range("Vector index out of range");
    }
    return data[i];
}

template <typename T>
Vector<T>& Vector<T>::operator/=(T scalar) {
    if (scalar == T(0)) {
        throw std::invalid_argument("Division by zero");
    }
    for (size_t i = 0; i < size; ++i) {
        data[i] /= scalar;
    }
    return *this;
}

template <typename T>
Vector<T>& Vector<T>::operator-=(const Vector& rhs) {
    CheckSize(rhs);
    for (size_t i = 0; i < size; ++i) {
        data[i] -= rhs[i];
    }
    return *this;
}

template <typename T>
Vector<T>& Vector<T>::operator+=(const Vector& rhs) {
    CheckSize(rhs);
    for (size_t i = 0; i < size; ++i) {
        data[i] += rhs[i];
    }
    return *this;
}

template <typename T>
Vector<T> Vector<T>::operator*(T scalar) const {
    Vector res(size);
    for (size_t i = 0; i < size; ++i) {
        res[i] = data[i] * scalar;
    }
    return res;
}

template <typename T>
void Vector<T>::Swap(Vector& other) noexcept {
    std::swap(data, other.data);
    std::swap(size, other.size);
}

template <typename T>
class Matrix : public Vector<Vector<T>> {
    size_t rows{ 0 };
    size_t cols{ 0 };
public:
    Matrix() = default;
    Matrix(size_t r, size_t c) : Vector<Vector<T>>(r), rows(r), cols(c) {
        for (size_t i = 0; i < rows; ++i) {
            (*this)[i] = Vector<T>(cols);
        }
    }
    ~Matrix() = default;
    Matrix(const Matrix& other) : Vector<Vector<T>>(other), rows(other.rows), cols(other.cols) {}

    size_t GetRows() const noexcept { return rows; }
    size_t GetCols() const noexcept { return cols; }    

    Matrix& operator=(const Matrix& other); 

    Vector<T>& operator[](size_t row);

    const Vector<T>& operator[](size_t row) const; 

    Matrix operator*(T scalar) const; 

    void Swap(Matrix& other) noexcept;

    void SwapRows(size_t i, size_t j); 

};

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix& other) {
    if (this == &other) {
        return *this;
    }
    Matrix temp(other);
    this->Swap(temp);
    return *this;
}

template <typename T>
Vector<T>& Matrix<T>::operator[](size_t row) {
    if (row >= rows) {
        throw std::out_of_range("Matrix row index out of range");
    }
    return Vector<Vector<T>>::operator[](row);
}

template <typename T>
const Vector<T>& Matrix<T>::operator[](size_t row) const {
    if (row >= rows) {
        throw std::out_of_range("Matrix row index out of range");
    }
    return Vector<Vector<T>>::operator[](row);
}

template <typename T>
Matrix<T> Matrix<T>::operator*(T scalar) const {
    Matrix res(*this);
    for (size_t i = 0; i < res.GetRows(); ++i) {
        res[i] *= scalar;
    }
    return res;
}

template <typename T>
void Matrix<T>::Swap(Matrix& other) noexcept {
    Vector<Vector<T>>::Swap(other);
    std::swap(rows, other.rows);
    std::swap(cols, other.cols);
}

template <typename T>
void Matrix<T>::SwapRows(size_t i, size_t j) {
    if (i >= rows || j >= rows) {
        throw std::out_of_range("Row indices out of range");
    }
    if (i != j) {
        (*this)[i].Swap((*this)[j]);
    }
}

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
Solution<T> GaussMethod(const Matrix<T>& _A, const Vector<T>& _b) {
    Matrix<T> A = _A;
    Vector<T> b = _b;
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

        A.SwapRows(rank, max_row);
        std::swap(b[rank], b[max_row]);

        T pivot = A[rank][col];
        A[rank] /= pivot;
        b[rank] /= pivot;

        for (size_t i = rank + 1; i < m; ++i) {
            T factor = A[i][col];    
            if (std::abs(A[i][col]) > eps) {
                A[i] -= A[rank] * factor;
                b[i] -= b[rank] * factor;
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
            sol.particular[col] -= A[i][j] * sol.particular[j];
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
T MaxAbs(const Vector<T>& v) {
    if (v.GetSize() == 0) {
        throw std::runtime_error("Cannot find max of empty vector!");
    }
    T max = std::abs(v[0]);
    for (size_t i = 1; i < v.GetSize(); ++i) {
        T curr = std::abs(v[i]);
        if (curr > max) {
            max = curr;
        }
    }
    if (max <= 0) {
        return 0;
    }
    return max;
}

template <typename T>
void CheckSolution(const Matrix<T>& A, const Vector<T>& b, const Solution<T>& sol) {
    const T eps = std::numeric_limits<T>::epsilon() * 1000;
    Vector<T> res(b.GetSize());
    
    if (!sol.is_consistent) {
        std::cout << "Solution is correct!\n";
        return;
    }

    for (size_t i = 0; i < A.GetRows(); ++i) {
        T sum = 0;
        for (size_t j = 0; j < A.GetCols(); ++j) {
            sum += A[i][j] * sol.particular[j];
        }
        res[i] = std::abs(sum - b[i]);
    }

    T max_error = MaxAbs(res);
    T max_b = MaxAbs(b); 

    std::cout << "\n=== Solution verification ===\n";
    std::cout << "Absolute error: " << max_error << "\n";
    std::cout << "Relative error: " << max_error / max_b << "\n";
    std::cout << "Epsilon scale: " << eps << "\n";

    if (max_error > eps * std::max(T(1), max_b)) {
        throw std::runtime_error("Solution is numerically unstable!");
    }
    std::cout << "Solution is correct!\n";
}

template <typename T>
void PrintSolution(const Solution<T>& sol) {
    const T eps = std::numeric_limits<T>::epsilon() * 100;

    if (!sol.is_consistent) {
        std::cout << "System is inconsistent (no solution exists).\n";
        return;
    }

    std::cout << "=== Solution ===\n";

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

template <typename T>
void ProcessRandomMatrix() {
    size_t n;
    std::cout << "Enter matrix size (n x n): ";
    std::cin >> n;

    if (n == 0) {
        throw std::invalid_argument("Matrix size must be non-zero!");
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    T max_val = static_cast<T>(10000.0 / std::sqrt(n));
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
        A[i][i] = row_sum;
        b[i] = dist(gen);
    }

    const auto start{ std::chrono::steady_clock::now() };

    Solution<T> sol = GaussMethod(A, b);

    const auto finish{ std::chrono::steady_clock::now() };
    const std::chrono::duration<double> elapsed_seconds{ finish - start };

    CheckSolution(A, b, sol);
    std::cout << "\n=== Time of work ===\n";
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
        throw std::invalid_argument("Matrix size must be non-zero!");
    }

    Matrix<T> A(m, n);
    Vector<T> b(m);

    std::cout << "Enter augmented matrix of coefficients:\n";
    for (size_t i = 0; i < m; ++i) {
        std::cout << "Row " << i + 1 << ": ";
        for (size_t j = 0; j < n; ++j) {
            if (!(std::cin >> A[i][j])) {
                throw std::invalid_argument("Invalid input!");
            }
        }
        if (!(std::cin >> b[i])) {
            throw std::invalid_argument("Invalid input!");
        }
    }
    Solution<T> sol = GaussMethod(A, b);

    CheckSolution(A, b, sol);
    PrintSolution(sol);
}

enum Choice {
    RandomMatrix = 1,
    ManualInput = 2,
    Exit = 3
};

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
        throw std::invalid_argument("Invalid choice!");
    }
    std::cout << std::endl;
}

enum TypeChoice {
    FloatType = 1,
    DoubleType = 2
};


void FinalRun() {
    std::cout << "Select numeric type:\n";
    std::cout << "1. float (single precision)\n";
    std::cout << "2. double (double precision)\n";
    std::cout << "3. Exit\n";
    std::cout << "Your choice: ";
    size_t type_choice;
    std::cin >> type_choice;

    switch (type_choice) {
    case FloatType:
        Run<float>();
        break;
    case DoubleType:
        Run<double>();
        break;
    case Exit:
        std::cout << "Exiting program...\n";
        break;
    default:
        throw std::invalid_argument("Invalid choice!");
    }
}