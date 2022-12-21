#include <iostream>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <sstream>

class Grid final {
    size_t y_size, x_size;
    std::vector < std::vector <int> > v;
public:
    Grid() : x_size(0), y_size(0){
        v.resize(0);
    }

    explicit Grid(const int &t): x_size(1), y_size(1) {
        v.resize(1);
        v[0].resize(1);
        v[0][0] = t;
    }

    Grid(size_t y_size, size_t x_size) : y_size(y_size), x_size(x_size) {
        v.resize(y_size);
        for (int i = 0; i < y_size; i++) {
            v[i].resize(x_size);
        }
    }

    int operator()(size_t x, size_t y) const {
        return v[y][x];
    }
    //Grid(std::vector < std::vector <int> > v, size_t y_size, size_t x_size) :    // сетка с передаваемыми данными
    //        v(v), y_size(y_size), x_size(x_size) {}

    void write(std::string filename) {
        std::ifstream input;
        input.open(filename);
        for (int i = 0; i < y_size; i++) {
            for (int j = 0; j < x_size; j++) {
                long long int a, b, c, d;
                input >> a;
                input >> b;
                input >> c;
                input >> d;
                v[i][j].a_num = a;
                v[i][j].a_div = b;
                v[i][j].b_num = c;
                v[i][j].b_div = d;
            }
        }
        input.close();
    }

    void print_matrix() {
        for (int i = 0; i < y_size; i++) {
            for (int j = 0; j < x_size; j++) {
                v[i][j].print();
                cout << ",";
            }
            cout << endl;
        }
    }

};


int main() {
    Grid m(3, 3);


    return 0;
}
