#include <iostream>
#include <cassert>
#include <cstdint>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <sstream>

#define MATRIX_FILE "C:/Users/Vova/CLionProjects/gauss/input.txt"

using namespace std;

struct RationalComplex {
    RationalComplex(RationalComplex const &complex) {}

    long long int a_num, b_num, a_div, b_div;

    RationalComplex() : a_num(0), b_num(0), a_div(1), b_div(1) {};

    RationalComplex(long long int a_num, long long int b_num, long long int a_div, long long int b_div):
                    a_num(a_num), b_num(b_num), a_div(a_div), b_div(b_div) {};

    RationalComplex(RationalComplex &var) = default; //TODO


    void operator+=(RationalComplex var) { //TODO: вынести из класса
        a_num = a_num * var.a_div + a_div * var.a_num;
        a_div = a_div * var.a_div;
        b_num = b_num * var.b_div + b_div * var.b_num;
        b_div = b_div * var.b_div;
    }

    void operator-=(RationalComplex var) {
        a_num = a_num * var.a_div - a_div * var.a_num;
        a_div = a_div * var.a_div;
        b_num = b_num * var.b_div - b_div * var.b_num;
        b_div = b_div * var.b_div;
    }

    void operator*=(RationalComplex var) {
        long long int a = a_num;
        long long int b = b_num;
        long long int c = a_div;
        long long int d = b_div;

        a_num = a * var.a_num * d * var.b_div - b * var.b_num * c * var.a_div;
        a_div = c * var.a_div * d * var.b_div;
        b_num = a * var.b_num * d * var.a_div + b * var.a_num * c * var.b_div;
        b_div = c * var.a_div * d * var.b_div;
        reduce();
    }

    void operator/=(RationalComplex var) { //TODO: division by 0, reduce every dividing
        long long int a = a_num;
        long long int b = b_num;
        long long int c = a_div;
        long long int d = b_div;

        long long int div = var.a_num * var.b_div * var.a_num * var.b_div + var.b_num * var.a_div * var.b_num * var.a_div;
        long long int denum = var.a_div * var.b_div * var.a_div * var.b_div;

        a_num = (a * var.a_num * d * var.b_div + b * var.b_num * c * var.a_div) * denum;
        b_num = (-a * var.b_num * d * var.a_div + b * var.a_num * c * var.b_div) * denum;
        a_div = (c * var.a_div * d * var.b_div) * div;
        b_div = (c * var.a_div * d * var.b_div) * div;
        reduce();
    }

    bool operator==(RationalComplex var) {
        return (a_num * var.a_div == a_div * var.a_num && b_num * var.b_div == b_div * var.b_num);
    }

    void reduce() {
        long long int t;
        long long int a1 = abs(a_num);
        long long int a2 = abs(a_div);
        long long int b1 = abs(b_num);
        long long int b2 = abs(b_div);
        //std::cout << "@@@";
        //std::cout << "[" << a1 << "!" << a2 << "]";
        while (a2 != 0) {
            //std::cout << "[" << a1 << "!" << a2 << "]";
            t = a2;
            a2 = a1 % a2;
            a1 = t;
        }
        //std::cout << "#";
        a_num /= a1;
        a_div /= a1;
        //std::cout << "$$$";
        while (b2 != 0) {
            //std::cout << "[" << b1 << "!" << b2 << "]";
            t = b2;
            b2 = b1 % b2;
            b1 = t;
        }
        b_num /= b1;
        b_div /= b1;
        //std::cout << std::endl;
    }

    void print() { //TODO: make const
        //std::cout << "*(" << a_num << "/" << a_div << ")" << "+i" << "(" << b_num << "/" << b_div << ")* -> ";
        reduce();
        if (a_num != 0) {
            if ((a_num < 0 && a_div > 0) || (a_num > 0 && a_div < 0)) {
                std::cout << "-";
            }
            if (abs(a_div) == 1) {
                std::cout << abs(a_num);
            } else{
                std::cout << "(" << abs(a_num) << "/" << abs(a_div) << ")";
            }
        }
        if (b_num != 0) {
            if ((b_num > 0 && b_div > 0) || (b_num < 0 && b_div < 0)) {
                std::cout << "+";
            } else {
                std::cout << "-";
            }
            if (abs(b_div) == 1) {
                std::cout << abs(b_num) << "i";
            } else{
                std::cout << "(" << abs(b_num) << "/" << abs(b_div) << ")i";
            }
        }
        if (is_zero()) {
            std::cout << "0";
        }
    }

    void check() { //TODO: перенести в конструктор
        if (a_div == 0) {
            a_div = 1;
            std::cout << "You can't divide by 0!!!! I turned the divider into 1" << std::endl;
        }
        if (b_div == 0) {
            b_div = 1;
            std::cout << "You can't divide by 0!!!! I turned the divider into 1" << std::endl;
        }
    }

    void write() { //TODO: delete
        //std::cout << "a[" << i + 1 << "][" << j + 1 << "]=";
        std::cin >> this->a_num;
        std::cout << "/";
        std::cin >> this->a_div;
        std::cout << "+i";
        std::cin >> this->b_num;
        std::cout << "/";
        std::cin >> this->b_div;

    }

    bool is_zero() const {
        return a_num == 0 && b_num == 0;
    }
};


template<typename T, size_t dim = 2>        //начальная размерность равна 2
class Grid final { //TODO: final в других местах и другие правила из файла
    using SGrid = Grid<T, dim - 1>;         //тип подсетки на 1 размерность меньше
    SGrid *data;                            //массив подсеток
    size_t ax_size;                         //длина высшей строчки
public:
    Grid() : data(nullptr), ax_size(0) {}   //пустая сетка, дефолтный конструктор

    ~Grid() {
        delete[] data;                      // деструктор
    }

    template<class ... sizes>
    Grid(const T &val, size_t ax_size, sizes ... args) : ax_size(ax_size) {         //pack is intended to be at the end
        data = new SGrid[ax_size];                  //создание массива подсеток
        for (int i = 0; i < ax_size; i++)           //для каждой подсетки
            data[i] = SGrid(val, args...);          //передаём оствашиеся размеры в такой же конструктор
    }

    Grid(const Grid<T, dim> &cpy) : ax_size(cpy.ax_size) {  //копирование объекта
        data = new SGrid[ax_size];
        for (int i = 0; i < ax_size; i++)
            data[i] = cpy.data[i];
    }

    Grid(Grid<T, dim> &&cpy) : ax_size(cpy.ax_size) { //перемещение объекта
        data = cpy.data;
        cpy.data = nullptr;
        cpy.ax_size = 0;
    }

    Grid<T, dim> &operator=(const Grid<T, dim> &rha) {   //  отождествление сеток
        if (&rha == this)
            return *this;
        if (ax_size != rha.ax_size) {
            ax_size = rha.ax_size;
            delete[] data;
            data = new SGrid[ax_size];
        }
        for (int i = 0; i < ax_size; i++)
            data[i] = rha.data[i];
        return *this;
    }

    Grid<T, dim> &
    operator=(Grid<T, dim> &&rha) {  //перенос данных с другой сетки с последующим удалением исходных сеток
        if (&rha == this)
            return *this;
        ax_size = rha.ax_size;
        delete[] data;
        data = rha.data;
        rha.data = nullptr;
        rha.ax_size = 0;
        return *this;
    }

    const SGrid &operator[](size_t i) const {
        return data[i];
    }

    SGrid &operator[](size_t i) {
        return data[i];
    }

    template<class ... indexes>
    const T &operator()(size_t i, indexes ... args) const {
        return data[i](args...);            //выбираем подсетку и предаём остальные индексы дальше оператору ()
    }

    template<class ... indexes>
    T &operator()(size_t i, indexes ... args) {
        return data[i](args...);
    }
};


template<typename T>
class Grid<T, 2> final {
private:
    size_t y_size, x_size;
    T *data;
public:
    Grid() : data(nullptr), x_size(0), y_size(0) {}  // конструктор по умолчанию

    ~Grid() {
        delete[] data;          // деструктор
    }

    Grid(T *data, size_t y_size, size_t x_size) :    // сетка с передаваемыми данными
            data(data), y_size(y_size), x_size(x_size) {}

    Grid(const Grid<T> &rha) : y_size(rha.y_size), x_size(rha.x_size) {      // копирование объекта
        data = new T[x_size * y_size];
        for (int i = 0; i < x_size * y_size; i++)
            data[i] = rha.data[i];
    }

    Grid(Grid<T> &&rha) : y_size(rha.y_size), x_size(rha.x_size) {   //перемещение объекта
        data = rha.data;
        rha.data = nullptr;
        rha.x_size = 0;
        rha.y_size = 0;
    }

    Grid(const T &t) : x_size(1), y_size(1) {        //сетка 1Х1 с переданным элементом
        data = new T[1];
        *data = t;
    }

    Grid(size_t y_size, size_t x_size) : y_size(y_size),
                                         x_size(x_size) {            //сетка с переданными размерами и дефолтными значеиями
        data = new T[y_size * x_size]();
    }

    Grid(const T &t, size_t y_size, size_t x_size) : y_size(y_size),
                                                     x_size(x_size) {    //сетка с переданными размерами и значеиями
        data = new T[y_size * x_size];
        for (int i = 0; i < x_size * y_size; i++)
            data[i] = t;
    }

    Grid<T> &operator=(const Grid<T> &rha) {   //  отождествление сеток
        if (&rha == this)
            return *this;
        delete data;
        x_size = rha.x_size;
        y_size = rha.y_size;
        data = new T[x_size * y_size];
        for (int i = 0; i < x_size * y_size; i++)
            data[i] = rha.data[i];
        return *this;
    }

    Grid<T> &operator=(Grid<T> &&rha) {    //перенос данных с другой сетки с последующим удалением исходных сеток
        if (&rha == this)
            return *this;
        delete data;
        x_size = rha.x_size;
        y_size = rha.y_size;
        data = rha.data;
        rha.data = nullptr;
        rha.x_size = 0;
        rha.y_size = 0;
        return *this;
    }

    T operator()(size_t y, size_t x) const {
        return data[y * x_size + x];
    }

    T &operator()(size_t y, size_t x) {
        return data[y * x_size + x];
    }

    Grid<T> &operator=(const T &t) {   //замена объектов сетки на переданные
        for (auto it = data, end = data + x_size * y_size; it != end; it++)
            *it = t;
        return *this;
    }

    T *operator[](size_t y) {
        return data + x_size * y;           //возвращает T* в начале y-й строки
    }

    const T *operator[](size_t y) const {
        return data + x_size * y;
    }

    size_t get_y_size() const { return y_size; }

    size_t get_x_size() const { return x_size; }


    Grid<RationalComplex> write(std::string filename, int rows, int columns) {
        Grid<RationalComplex> matrix = Grid<RationalComplex>(rows, columns);
        std::ifstream input;
        input.open(filename);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                long long int a, b, c, d;
                input >> a;
                cout << a;
                input >> b;
                cout << b;
                input >> c;
                cout << c;
                input >> d;
                cout << d;
                matrix[i][j].a_num = a;
                matrix[i][j].a_div = b;
                matrix[i][j].b_num = c;
                matrix[i][j].b_div = d;
            }
        }
        input.close();
        return matrix;
    }

    void print_matrix() {
        for (int i = 0; i < y_size; i++) {
            for (int j = 0; j < x_size; j++) {
                (*this)[i][j].print();
                cout << ",";
            }
            cout << endl;
        }
    }

};

class Matrix final {
public:
    Matrix() {};

    void
    division(Grid<RationalComplex> &matrix, int num) { //деление строчки, чтобы первое ненулевое число стало единичкой
        for (int j = 0; j < matrix.get_x_size(); j++) {
            if (not matrix[num][j].is_zero()) {  //TODO
                for (int k = j + 1; k < matrix.get_x_size(); k++) {
                    matrix[num][k] /= matrix[num][j];
                    std::cout << "Dividing" << std::endl;
                    matrix.print_matrix();
                }
                matrix[num][j] /= matrix[num][j];
                std::cout << "Dividing" << std::endl;
                matrix.print_matrix();
                break;
            }
        }
    }

    void
    subtracting(Grid<RationalComplex> &matrix, int num1, int num2) { //вычитание из второй строки первой, домноженной
        // на такой коэффициент, чтобы первый ненулевой элемент вычитаемой строчки обнулился
        for (int j = 0; j < matrix.get_x_size(); j++) {
            if (not matrix[num1][j].is_zero()) {
                for (int k = j + 1; k < matrix.get_x_size(); k++) {
                    RationalComplex help = matrix[num2][j];
                    help /= matrix[num1][j]; //TODO mb optimize
                    help *= matrix[num1][k];
                    matrix[num2][k] -= help;
                    std::cout << "Subtracting" << std::endl;
                    //matrix.print_matrix();
                }
                matrix[num2][j] -= matrix[num2][j];
                std::cout << "Subtracting" << std::endl;
                //matrix.print_matrix();
                break;
            }
        }
    }

    void reverse(Grid<RationalComplex> &matrix, int num1, int num2) {
        //Меняет две строчки
        for (int j = 0; j < matrix.get_x_size(); j++) {
            RationalComplex help = matrix[num1][j];
            matrix[num1][j] = matrix[num2][j];
            matrix[num2][j] = help;
            std::cout << "Reversing" << std::endl;
        }
    }

    int solve(Grid<RationalComplex> &matrix) {
        //Приводит матрицу к диагональному виду, выводит ранг
        int rank = 0;
        for (int step = 0; step < matrix.get_x_size(); step++) { //прямой ход
            if (matrix[rank][step].is_zero()) {
                for (int i = step; i < matrix.get_y_size(); i++) { //поиск ненулевого и замена
                    if (not matrix[i][step].is_zero()) {
                        reverse(matrix, step, i);
                        break;
                    }

                }
                continue; //к следующему столбцу, если первый нулевой
            }
            division(matrix, step);
            for (int j = step + 1; j < matrix.get_y_size(); j++) {
                subtracting(matrix, step, j);
            }
            rank += 1;
        }
        for (int step = matrix.get_y_size() - 1; step > 0; step--) { //обратный ход
            for (int substep = step - 1; substep >= 0; substep--) {
                subtracting(matrix, step, substep);
            }
        }
        for (int key = 0; key < matrix.get_y_size(); key++) {
            division(matrix, key);
        }
        return rank;
    }

    Grid<RationalComplex> final_solver(Grid<RationalComplex> &matrix) { //TODO: dim=0
        //Выдает матрицу - фундаментальную систему решений однородного уравнения
        int rank = solve(matrix);
        Grid<RationalComplex> answer = Grid<RationalComplex>(matrix.get_x_size(), matrix.get_x_size() - rank);
        assert(matrix.get_x_size() >= rank);
        for (int i = 0; i < rank; i++) {
            for (int j = 0; j < matrix.get_x_size() - rank; j++) {
                RationalComplex help;
                help.a_num = 0;
                help.b_num = 0;
                help.a_div = 1;
                help.b_div = 1;
                help -= matrix[i][j + rank];
                answer[i][j] = help;
            }
        }
        for (int i = rank; i < matrix.get_x_size(); i++) {
            RationalComplex help;
            help.a_num = 1;
            help.b_num = 0;
            help.a_div = 1;
            help.b_div = 1;
            answer[i][i - rank] = help;
        }
        return answer;

    }

    Grid<RationalComplex> eigenvectors(Grid<RationalComplex> matrix, Grid<RationalComplex> eigenvalues) {
        //Нахождение собственных векторов по известным собственным значениям
        assert(matrix.get_x_size() == matrix.get_y_size());
        assert(eigenvalues.get_y_size() <= matrix.get_y_size());
        Grid<RationalComplex> answer = Grid<RationalComplex>(matrix.get_y_size(), matrix.get_y_size());
        Grid<RationalComplex> help = matrix;
        int x = 0;
        for (int i = 0; i < eigenvalues.get_y_size(); i++) {
            matrix = help;
            for (int j = 0; j < matrix.get_y_size(); j++) {
                matrix[j][j] -= eigenvalues[0][i];
            }
            Grid<RationalComplex> eigen = final_solver(matrix);

            for (int j = 0; j < eigen.get_y_size(); j++) {
                for (int k = 0; k < eigen.get_x_size(); k++) {
                    answer[j][k + x] = eigen[j][k];
                }
            }
            x += eigen.get_x_size();
        }
        return answer;
    }

    void JNF(Grid<RationalComplex> matrix, Grid<RationalComplex> eigenvalues) {
        //Функция, позволяющая понять, как выглядит ЖНФ для заданной матрицы 2x2
        assert(matrix.get_y_size() == 2 && matrix.get_x_size() == 2);
        if (eigenvalues.get_y_size() == 2) {
            eigenvalues[0][0].print();
            cout << " 0";
            cout << endl;
            cout << "0 ";
            eigenvalues[0][1].print();
        } else if (eigenvectors(matrix, eigenvalues).get_x_size() == 2) {
            eigenvalues[0][0].print();
            cout << " 0";
            cout << endl;
            cout << "0 ";
            eigenvalues[0][0].print();
        } else {
            eigenvalues[0][0].print();
            cout << " 1";
            cout << endl;
            cout << "0 ";
            eigenvalues[0][0].print();
        }
    }

    Grid<RationalComplex> transposition(Grid<RationalComplex> matrix) {
        //Транспонирование матрицы
        Grid<RationalComplex> answer = Grid<RationalComplex>(matrix.get_x_size(), matrix.get_y_size());
        for (int i = 0; i < answer.get_y_size(); i++) {
            for (int j = 0; j < answer.get_x_size(); j++) {
                answer[i][j] = matrix[j][i];
            }
        }
        return answer;
    }

    Grid<RationalComplex> duality(Grid<RationalComplex> matrix) {
        //Нахождение матрицы по ФСР - использование "дуальности"
        Grid<RationalComplex> help = transposition(matrix);
        Grid<RationalComplex> answer = final_solver(help);
        return transposition(answer);
    }

    RationalComplex trace(Grid<RationalComplex> matrix) {
        //Нахождение следа матрицы
        assert(matrix.get_y_size() == matrix.get_x_size());
        RationalComplex trace = matrix[0][0];
        for (int i = 1; i < matrix.get_y_size(); i++) {
            trace += matrix[i][i];
        }
        return trace;
    }

    Grid<RationalComplex> sum(Grid<RationalComplex> matrix1, Grid<RationalComplex> matrix2) {
        // Сумма двух матриц (с одинаковыми размерами)
        assert(matrix1.get_y_size() == matrix2.get_y_size() && matrix1.get_x_size() == matrix2.get_x_size());
        Grid<RationalComplex> matrix = Grid<RationalComplex>(matrix1.get_y_size(), matrix1.get_x_size());
        for (int i = 0; i < matrix.get_y_size(); i++) {
            for (int j = 0; j < matrix.get_x_size(); j++) {
                matrix[i][j] = matrix1[i][j];
                matrix[i][j] += matrix2[i][j];
            }
        }
        return matrix;

    }

    Grid<RationalComplex> product(Grid<RationalComplex> matrix1, Grid<RationalComplex> matrix2) {
        // Произведение двух матриц (с подходящими размерами)
        assert(matrix1.get_x_size() == matrix2.get_y_size());
        Grid<RationalComplex> matrix = Grid<RationalComplex>(matrix1.get_y_size(), matrix2.get_x_size());
        for (int i = 0; i < matrix.get_y_size(); i++) {
            for (int j = 0; j < matrix.get_x_size(); j++) {
                matrix[i][j].a_num = 0;
                matrix[i][j].b_num = 0;
                matrix[i][j].a_div = 1;
                matrix[i][j].b_div = 1;
                for (int k = 0; k < matrix.get_x_size(); k++) {
                    RationalComplex help = matrix1[i][k];
                    help *= matrix2[k][j];
                    matrix[i][j] += help;
                }
            }
        }
        return matrix;
    }


};

int main() { //TODO create other code files

    int rows;
    int columns;
    cout << "Enter amount of rows:";
    cin >> rows;
    cout << "Enter amount of columns:";
    cin >> columns;

    Grid<RationalComplex> matrix = Grid<RationalComplex>(rows, columns);
    matrix = matrix.write(MATRIX_FILE, rows, columns); //TODO excess arguments
    Grid<RationalComplex> &matr = matrix;

    matrix.print_matrix();

    /*   for (int i = 0; i < rows; i++) {
           for (int j = 0; j < columns; j++) {
               cout << "a[" << i + 1 << "," << j + 1 << "]=";
               matrix[i][j].write();
               matrix[i][j].check();
           }
       }*/

    Matrix solver = Matrix();

    Grid<RationalComplex> answer = solver.final_solver(matr);
    matrix.print_matrix();

    answer.print_matrix();

    return 0;
}
