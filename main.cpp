#include <iostream>
#include <vector>

template<typename T>
class Matrix {
 private:
    std::vector<std::vector<T>> matrix;
    std::vector<T> elements;
 public:
    explicit Matrix(const std::vector<std::vector<T>> &v) {
        matrix = v;
        for (const auto &vec : v) {
            for (auto el : vec) {
                elements.push_back(el);
            }
        }
    }

    friend std::ostream &operator<<(std::ostream &out, const Matrix<T> &m) {
        for (size_t row = 0; row < m.size().first - 1; ++row) {
            for (size_t col = 0; col < m.size().second; ++col) {
                if (col != m.size().second - 1) {
                    out << m.matrix[row][col] << "\t";
                } else {
                    out << m.matrix[row][col];
                }
            }
            out << "\n";
        }
        for (size_t col = 0; col < m.size().second; ++col) {
            if (col != m.size().second - 1) {
                out << m.matrix[m.size().first - 1][col] << "\t";
            } else {
                out << m.matrix[m.size().first - 1][col];
            }
        }

        return out;
    }

    Matrix &operator+=(const Matrix &other) {
        for (size_t row = 0; row < size().first; ++row) {
            for (size_t col = 0; col < size().second; ++col) {
                matrix[row][col] += other.matrix[row][col];
            }
        }

        return *this;
    }

    Matrix operator+(const Matrix &other) const {
        Matrix result(matrix);
        for (size_t row = 0; row < size().first; ++row) {
            for (size_t col = 0; col < size().second; ++col) {
                result.matrix[row][col] = matrix[row][col] + other.matrix[row][col];
            }
        }
        return result;
    }

    template<typename U>
    Matrix &operator*=(const U &number) {
        for (size_t row = 0; row < size().first; ++row) {
            for (size_t col = 0; col < size().second; ++col) {
                matrix[row][col] *= number;
            }
        }
        return *this;
    }

    template<typename U>
    Matrix operator*(const U &number) const {
        Matrix result(matrix);
        for (size_t row = 0; row < size().first; ++row) {
            for (size_t col = 0; col < size().second; ++col) {
                result.matrix[row][col] *= number;
            }
        }
        return result;
    }

    Matrix &transpose() {
        std::vector<std::vector<T>> v(size().second, std::vector<T>(size().first, T()));
        Matrix transposed(v);
        for (size_t row = 0; row < transposed.size().first; ++row) {
            for (size_t col = 0; col < transposed.size().second; ++col) {
                transposed.matrix[row][col] = matrix[col][row];
            }
        }
        matrix = transposed.matrix;
        return *this;
    }

    [[nodiscard]] Matrix transposed() const {
        std::vector<std::vector<T>> v(size().second, std::vector<T>(size().first, T()));
        Matrix transposed(v);
        for (size_t row = 0; row < transposed.size().first; ++row) {
            for (size_t col = 0; col < transposed.size().second; ++col) {
                transposed.matrix[row][col] = matrix[col][row];
            }
        }
        return transposed;
    }

    Matrix &operator*=(const Matrix &other) {
        assert(other.size().first == size().second);
        std::vector<std::vector<T>> v(size().first, std::vector<T>(other.size().second, T()));
        auto first = matrix;
        matrix = v;
        for (size_t row = 0; row < size().first; ++row) {
            for (size_t col = 0; col < size().second; ++col) {
                auto row_to_work_with = first[row];
                for (size_t column = 0; column < row_to_work_with.size(); ++column) {
                    matrix[row][col] += (row_to_work_with[column] * other.matrix[column][col]);
                }
            }
        }
        return *this;
    }

    Matrix operator*(const Matrix &other) const {
        assert(other.size().first == size().second);
        std::vector<std::vector<T>> v(size().first, std::vector<T>(other.size().second, T()));
        Matrix<T> new_one(v);
        for (size_t row = 0; row < new_one.size().first; ++row) {
            for (size_t col = 0; col < new_one.size().second; ++col) {
                auto row_to_work_with = matrix[row];
                for (size_t column = 0; column < row_to_work_with.size(); ++column) {
                    new_one.matrix[row][col] +=
                        (row_to_work_with[column] * other.matrix[column][col]);
                }
            }
        }
        return new_one;
    }

    using iter = typename std::vector<T>::iterator;
    using const_iter = typename std::vector<T>::const_iterator;

    iter begin() {
        return elements.begin();
    }

    iter end() {
        return elements.end();
    }

    [[nodiscard]] const_iter begin() const {
        return elements.begin();
    }

    [[nodiscard]] const_iter end() const {
        return elements.end();
    }

    [[nodiscard]] std::pair<size_t, size_t> size() const {
        return std::make_pair(matrix.size(), matrix[0].size());
    }
};
