#include <mpi.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>

#include <iostream>
#include <random>
#include <vector>

#define MAT_LENGTH 10000
#define MAT_WIDTH 1000

class substr {

    int sx, sy, x,y;
public:
    substr(int sx, int sy, int x, int y): sx(sx), sy(sy), x(x), y(y){}
    boost::basic_format<char> toString() const {
        return boost::format("[%1%, %2%, %3%, %4% | %5%, %6% ]") % sx % sy % x % y % (x - sx) % (y - sy);
    }
};

std::ostream &operator<<(std::ostream &os, substr const &m) {
    return os << m.toString();
}

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
    if ((y - sy) == 0){
        return 0;
    } else if (y - sy > 0) {
        return 1;
    } else {
        return -1;
    }
}
void processResult(std::vector<substr> &result, int sx, int sy, int x, int y){
    int length = x - sx;
    if (length >= 2){
        result.push_back(substr(sx, sy, x, y));
    }
}

bool dupa(boost::numeric::ublas::matrix<char> &mat, std::vector<substr> &result, int sx, int sy, int x, int y){
//    std::cout << boost::format("--| %1% %2% %3% %4% \n") % sx % sy % x % y ;
    bool finished = false;
    if (isWithinRange(mat, x, y)) {
        if (cool(mat(x, y))){
            int next_x = x+1, next_y = y + yVector(sx, sy, x, y);
            if (dupa(mat, result, sx, sy, next_x, next_y)){
                processResult(result, sx, sy, x, y);
            }
        } else {
            finished = true;

        }
    } else {
        finished = true;
    }
    return finished;
}

void handle(boost::numeric::ublas::matrix<char> &mat, std::vector<substr> &result, int x, int y){
    if (isWithinRange(mat, x, y) && cool(mat(x,y))){
        dupa(mat, result, x, y, x+1, y-1);
        dupa(mat, result, x, y, x+1, y);
        dupa(mat, result, x, y, x+1, y+1);
    }
}

int main(/*int argc, char **argv*/)
{
    std::random_device rd;
    std::mt19937 e2(rd());

    std::uniform_int_distribution<char> dist(48,49);
    boost::numeric::ublas::matrix<char> mat(MAT_LENGTH, MAT_WIDTH);

    e2.seed(10);

    std::vector<substr> result;

    for(auto i = 0; i < mat.size1(); i++){
        for(auto j = 0; j < mat.size2(); j++){
            mat(i,j) = dist(e2);
        }
    }

    for(auto i = 0; i < mat.size1(); i++){
        for(auto j = 0; j < mat.size2(); j++){
            handle(mat, result, i, j);
        }
    }

//    std::cout << mat;
//    std::cout << "\n";
//    for (const auto i: result){
//      std::cout << i << '\n';
//    }
//    std::cout << "\n";



    return 0;
}
