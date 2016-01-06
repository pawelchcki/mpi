#include <mpi.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <random>
#include <vector>

#define MAT_LENGTH 10
#define MAT_WIDTH 10

bool cool(char x){
    if (x == '1'){
        return true;
    } else {
        return false;
    }
}
bool isWithinRange(boost::numeric::ublas::matrix<char> &mat, int x, int y){
    return x >= 0 && y >= 0 && x < mat.size1() && y < mat.size2();
}

int yVector(int sx, int sy, int x, int y){
    if (y - sy == 0){
        return 0;
    } else if (y - sy > 0) {
        return 1;
    } else {
        return -1;
    }
}

void dupa(boost::numeric::ublas::matrix<char> &mat, int sx, int sy, int x, int y){
    if (isWithinRange(mat, x, y) && cool(mat(x, y))){
        int next_x = x+1, next_y = x + yVector(sx, sy, x, y);

        if (isWithinRange(mat, next_x, next_y)){
            dupa(mat, sx, sy, next_x, next_y);
        }
    }
}

void handle(boost::numeric::ublas::matrix<char> &mat, int x, int y){
    if (isWithinRange(mat, x, y) && cool(mat(x,y))){
        dupa(mat, x, y, x+1, y-1);
        dupa(mat, x, y, x+1, y);
        dupa(mat, x, y, x+1, y+1);
    }
}

int main(/*int argc, char **argv*/)
{
    std::random_device rd;
    std::mt19937 e2(rd());

    std::uniform_int_distribution<char> dist(48,49);
    boost::numeric::ublas::matrix<char> mat(MAT_LENGTH, MAT_WIDTH);

    e2.seed(10);

    for(auto i = 0; i < mat.size1(); i++){
        for(auto j = 0; j < mat.size2(); j++){
            mat(i,j) = dist(e2);
        }
    }

    for(auto i = 0; i < mat.size1(); i++){
        for(auto j = 0; j < mat.size2(); j++){
            handle(mat, i, j);
        }
    }



    std::cout << mat;

    std::cout << "\n";

    auto dupa = "asf";
    std::cout <<dupa;

    return 0;
}
