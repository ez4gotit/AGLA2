#include<iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstring>


using namespace std;
class SquareMatrix;
class IdentityMatrix;
class EliminationMatrix;
class PermutationMatrix;

class Matrix {
public:

    Matrix() {
        this->rows = 0;
        this->columns = 0;
    }

    Matrix(int rows, int columns, vector<vector<double>> matrix) {
        this->rows = rows;
        this->columns = columns;

        vector<vector<double>> new_matrix(rows, vector<double>());
        for (int x = 0; x < rows; x++) {
            new_matrix[x] = vector<double>(columns, 0);
            for (int y = 0; y < columns; y++) {
                new_matrix[x][y] = matrix[x][y];
            }
        }

        this->matrix = new_matrix;
    }

    Matrix(int rows, int columns) {
        this->rows = rows;
        this->columns = columns;

        vector<vector<double>> new_matrix(rows, vector<double>());
        for (int x = 0; x < rows; x++) {
            new_matrix[x] = vector<double>(columns, 0);
            for (int y = 0; y < columns; y++) {
                new_matrix[x][y] = 0;
            }
        }
        this->matrix = new_matrix;
    }

    friend ostream& operator << (ostream& out, const Matrix& m) {
        vector<vector<double>> matrix = m.matrix;
        for (int x = 0; x < m.rows; x++) {
            out << matrix[x][0];
            for (int y = 1; y < m.columns; y++) {
                out << " " << matrix[x][y];
            }
            out << endl;
        }
        return out;
    }

    friend istream& operator >> (istream& in, Matrix& m) {
        in >> m.rows >> m.columns;
        vector<vector<double>> new_matrix(m.rows, vector<double>());
        for (int x = 0; x < m.rows; x++) {
            new_matrix[x] = vector<double>(m.columns, 0);
            for (int y = 0; y < m.columns; y++) {
                in >> new_matrix[x][y];
            }
        }
        m.matrix = new_matrix;
        return in;
    }

    Matrix operator+ (const Matrix& b) {

        if (rows == b.rows && columns == b.columns) {
            Matrix temp(rows, columns, matrix);
            for (int x = 0; x < rows; x++) {
                for (int y = 0; y < columns; y++) {
                    temp.matrix[x][y] += b.matrix[x][y];
                }
            }
            return temp;
        }
        Matrix temp;
        cout << "Error: the dimensional problem occurred" << endl;
        return temp;
    }

    Matrix operator- (const Matrix& b) {

        if (rows == b.rows && columns == b.columns) {
            Matrix temp(rows, columns, matrix);
            for (int x = 0; x < rows; x++) {
                for (int y = 0; y < columns; y++) {
                    temp.matrix[x][y] -= b.matrix[x][y];
                }
            }
            return temp;
        }
        Matrix temp;
        cout << "Error: the dimensional problem occurred" << endl;
        return temp;
    }

    Matrix operator* (const Matrix& b) {

        if (columns == b.rows) {
            Matrix temp(rows, b.columns);
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < b.columns; j++) {
                    for (int k = 0; k < columns; k++) {
                        temp.matrix[i][j] += matrix[i][k] * b.matrix[k][j];
                    }
                }
            }
            return temp;
        }
        Matrix temp;
        cout << "Error: the dimensional problem occurred" << endl;
        return temp;
    }

    /**
     * Calculating transpose matrix
     *
     * @return Resulting matrix of the operation
     */
    Matrix transpose()
    {
        Matrix temp(columns, rows);
        for (int i = 0; i < columns; i++)
            for (int j = 0; j < rows; j++)
                temp.matrix[i][j] = matrix[j][i];
        return temp;
    }

    Matrix& operator=(const Matrix& other) = default;
public:
    int rows;
    int columns;
    vector<vector<double>> matrix;
};


class SquareMatrix :public Matrix {
public:

    SquareMatrix() : Matrix() {
    }


    SquareMatrix(int size, vector<vector<double>> matrix) : Matrix(size, size, matrix) {
    }

    SquareMatrix(int size) : Matrix() {
        this->rows = size;
        this->columns = size;

        vector<vector<double>> new_matrix(size, vector<double>());
        for (int x = 0; x < size; x++) {
            new_matrix[x] = vector<double>(size, 0);
            for (int y = 0; y < size; y++) {
                new_matrix[x][y] = 0;
            }
        }
        this->matrix = new_matrix;
    }


    friend ostream& operator << (ostream& out, const SquareMatrix& m) {
        vector<vector<double>> matrix = m.matrix;
        for (int x = 0; x < m.rows; x++) {
            out << matrix[x][0];
            for (int y = 1; y < m.rows; y++) {
                out << " " << matrix[x][y];
            }
            out << endl;
        }
        return out;
    }


    friend istream& operator >> (istream& in, SquareMatrix& m) {
        in >> m.rows;
        m.columns = m.rows;
        vector<vector<double>> new_matrix(m.rows, vector<double>());
        for (int x = 0; x < m.rows; x++) {
            new_matrix[x] = vector<double>(m.rows, 0);
            for (int y = 0; y < m.rows; y++) {
                in >> new_matrix[x][y];
            }
        }
        m.matrix = new_matrix;
        return in;
    }



    SquareMatrix operator+ (const SquareMatrix& b) {
        Matrix newMatrix = (Matrix)*this + (Matrix)b;
        return *(SquareMatrix*)(&newMatrix);
    }


    SquareMatrix operator- (const SquareMatrix& b) {
        Matrix newMatrix = (Matrix)*this - (Matrix)b;
        return *(SquareMatrix*)(&newMatrix);
    }

    SquareMatrix operator* (const SquareMatrix& b) {
        Matrix newMatrix = (Matrix)*this * (Matrix)b;
        return *(SquareMatrix*)(&newMatrix);
    }


    SquareMatrix transpose()
    {
        Matrix newMatrix = (Matrix)*this;
        newMatrix = newMatrix.transpose();
        return *(SquareMatrix*)(&newMatrix);
    }

};



class IdentityMatrix : public SquareMatrix {
public:


    IdentityMatrix(int size) : SquareMatrix(size) {
        for (int x = 0; x < size; x++) {
            this->matrix[x][x] = 1;
        }
    }

};


class EliminationMatrix : public IdentityMatrix {
public:


    EliminationMatrix(int x, int y, SquareMatrix matrix) : IdentityMatrix(matrix.rows) {
        double value;
        double a = matrix.matrix[x - 1][y - 1];
        double k = matrix.matrix[y - 1][y - 1];

        value = -1 * a / k;

        this->matrix[x - 1][y - 1] = value;
    }
};


class PermutationMatrix : public IdentityMatrix {
public:

    PermutationMatrix(int x, int y, SquareMatrix matrix) : IdentityMatrix(matrix.rows) {
        vector<double> temp;

        temp = this->matrix[x - 1];
        this->matrix[x - 1] = this->matrix[y - 1];
        this->matrix[y - 1] = temp;
    }
};




void swap(int x, int y, SquareMatrix A)
{
    vector<double> temp;
    temp = A.matrix[x];
    A.matrix[x] = A.matrix[y];
    A.matrix[y] = temp;
}


double determinant(SquareMatrix A) {
    double answer = 1;
    for (int i = 0; i < A.rows; i++) answer *= A.matrix[i][i];
    return answer;
}


SquareMatrix inverseMatrix(SquareMatrix A) {
    SquareMatrix I = IdentityMatrix(A.rows);

    int count = 1;
    for (int k = 0; k < A.rows; k++)
    {
        int i_max = k;
        double v_max = abs(A.matrix[i_max][k]);

        for (int i = A.rows - 1; i > k; i--)
            if (abs(A.matrix[i][k]) > v_max)
                v_max = A.matrix[i][k], i_max = i;


        if (i_max != k) {
            swap(k, i_max, A);
            swap(k, i_max, I);

            count++;
        }


        for (int i = k + 1; i < A.rows; i++)
        {

            if (A.matrix[i][k] != 0) {
                EliminationMatrix E(i + 1, k + 1, A);
                A = E * A;
                I = E * I;

                count++;
            }

        }
    }

    for (int k = A.rows - 1; k > 0; k--)
    {
        for (int i = A.rows - 2 - (A.rows - 1 - k); i >= 0; i--)
        {

            if (A.matrix[i][k] != 0) {
                EliminationMatrix E(i + 1, k + 1, A);
                A = E * A;
                I = E * I;
                count++;
            }

        }
    }

    for (int i = 0; i < A.rows; i++) {
        double temp = 1 / A.matrix[i][i];
        for (int k = 0; k < I.rows; k++) {
            I.matrix[i][k] *= temp;
        }
        A.matrix[i][i] = 1;
    }

    return I;

}



Matrix leastSquareApproximation(Matrix A, Matrix b) {

    cout << "A:\n" << A;
    cout << "A_T*A:\n" << (A.transpose() * A);
    Matrix temp = A.transpose() * A;
    Matrix a_new = inverseMatrix(*(SquareMatrix*)(&temp));
    cout << "(A_T*A)^-1:\n" << a_new;
    Matrix b_new = A.transpose() * b;
    cout << "A_T*b:\n" << b_new;
    cout << "x~:\n";

    return a_new * b_new;

}
#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif
int main() {
#ifdef WIN32
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE* pipe = popen(GNUPLOT_NAME, "w");
#endif
    cout << fixed << setprecision(4);
    int size;
    cin >> size;

    Matrix b(size, 1);

    Matrix t(size, 1);


    for (int i = 0; i < size; i++) {
        cin >> t.matrix[i][0];
        cin >> b.matrix[i][0];
    }

    int degree;
    cin >> degree;
    Matrix A(size, degree + 1);

    for (int i = 0; i < size; i++) {
        for (int k = 0; k <= degree; k++) {
            A.matrix[i][k] = pow(t.matrix[i][0], k);
        }
    }
    Matrix X = leastSquareApproximation(A, b);
    cout << X;
    fprintf(pipe, "plot [-30 : 30] [-30 : 30] %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n", X.matrix[3][0], X.matrix[2][0], X.matrix[1][0], X.matrix[0][0]);
    for (int i = 0; i < A.rows; i++)
    {
        
        fprintf(pipe, "%f\t%f\n", t.matrix[i][0], b.matrix[i][0]);
    }
    fprintf(pipe, "e\n");
    fflush(pipe);
#ifdef WIN32
    _pclose(pipe);
#else
    pclose(pipe);
#endif
    return 0;
}
