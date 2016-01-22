#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/serialization/vector.hpp>


#include <iostream>
#include <random>
#include <vector>

#include <fstream>

#include <boost/archive/binary_oarchive.hpp>
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>

#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <cassert>

#define MAT_LENGTH 100
#define MAT_WIDTH 100

namespace mpi = boost::mpi;
namespace ublas = boost::numeric::ublas;

class OpVector{
    friend class boost::serialization::access;

    inline void generatevals(){
        length = abs(posx - startx);
        vx = (posx-startx)/length;
        vy = (posy-startx)/length;
        assert(vx == 0 || vx == 1 || vx == -1);
        assert(vy == 0 || vy == 1 || vy == -1);
    }

    template<class Archive>  void serialize(Archive & ar, const unsigned int version)
    {
       ar & startx;
       ar & starty;
       ar & posx;
       ar & posy;

    }
public:
    int length;
    char vx,vy; // vector to move
    int startx, starty;
    int posx, posy;

    OpVector(){}
    OpVector(int startx, int starty, int posx, int posy): startx(startx), starty(starty), posx(posx), posy(posy){
        generatevals();
    }
    boost::basic_format<char> toString() const {
        return boost::format("[%1%, %2%, %3%, %4% | %5%, %6%, %7% ]") % startx % starty % posx % posy % length % vx % vy;
    }
};

std::ostream &operator<<(std::ostream &os, OpVector const &m) {
    return os << m.toString();
}

bool isCorrect(char x){
    if (x == '1'){
        return true;
    } else {
        return false;
    }
}

inline bool isWithinRange(boost::numeric::ublas::matrix<char> &mat, int x, int y){
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
void processResult(std::vector<OpVector> &result, int sx, int sy, int x, int y){
    OpVector data(sx,sy,x,y);
    if (data.length >= 2){
        result.push_back(data);
    }
}

bool recurse(boost::numeric::ublas::matrix<char> &mat, std::vector<OpVector> &result, int sx, int sy, int x, int y){
    bool finished = false;
    if (isWithinRange(mat, x, y)) {
        if (isCorrect(mat(x, y))){
            int next_x = x+1, next_y = y + yVector(sx, sy, x, y);
            if (recurse(mat, result, sx, sy, next_x, next_y)){
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

class PositionHelper {
    int cellWidth = 2;
    int cellHeight = 2;
    int width;
    int cellsInWidth;
public:
    PositionHelper(int cellWidth, int cellHeight, int width) : cellWidth(cellWidth), cellHeight(cellHeight), width(width) {
        cellsInWidth = cellsInWidthCalc();
    }

    void processCell(ublas::matrix<char> &mat, std::vector<OpVector> &result, int cellNo) {
        int x = firstXInCell(cellNo);
        int y = firstYInCell(cellNo);

        for(int i=x; i <= (x + cellWidth); i++){
            for(int j=y; j <= (y + cellHeight); j++){
                if (isCorrect(mat(x,y))){
                    recurse(mat, result, x, y, x+1, y-1);
                    recurse(mat, result, x, y, x+1, y);
                    recurse(mat, result, x, y, x+1, y+1);
                }
            }
        }
    }

    inline int cellsInWidthCalc(){
        int result = width / cellWidth;
        if (width % cellWidth > 0){
            result += 1;
        }
        return result;
    }

    int cellNumber(int x, int y){
        int cellFromLeft = x / cellWidth;
        int cellFromTop = y / cellHeight;
        int rv = cellFromLeft + cellFromTop * cellsInWidth;
//        std::cout << rv << cellFromLeft << cellFromTop << std::endl;
        return rv;
    }

    int cellRank(int x, int y, int maxRank){
        return cellNumber(x, y) % maxRank;
    }

    int firstXInCell(int cellNo){
        return (cellNo % cellsInWidth) * cellWidth;
    }

    inline int firstYInCell(int cellNo){
        return (cellNo / cellsInWidth) * cellHeight;
    }
};

int test() {
    PositionHelper h3(2, 2, 3);

    assert(h3.cellNumber(0, 0) == 0);
    assert(h3.cellNumber(1, 0) == 0);
    assert(h3.cellNumber(2, 0) == 1);
    assert(h3.cellNumber(0, 1) == 0);
    assert(h3.cellNumber(0, 2) == 2);
    assert(h3.cellNumber(2, 2) == 3);
    assert(h3.firstXInCell(0) == 0);
    assert(h3.firstXInCell(1) == 2);
    assert(h3.firstXInCell(3) == 2);
    assert(h3.firstXInCell(2) == 0);

    assert(h3.firstYInCell(2) == 2);
    assert(h3.firstYInCell(0) == 0);
    assert(h3.firstYInCell(1) == 0);
    assert(h3.firstYInCell(3) == 2);
}

int main(int argc, char **argv){
    test();
    std::random_device rd;
    std::mt19937 e2(rd());

    std::uniform_int_distribution<char> dist(48,49);
    boost::numeric::ublas::matrix<char> mat(MAT_LENGTH, MAT_WIDTH);

    e2.seed(10);

    std::vector<OpVector> result;

    mpi::environment env;
    mpi::communicator world;
    std::cout << "I am process " << world.rank() << " of " << world.size()
              << "." << std::endl;

    int size,rank;

    for(auto i = 0; i < mat.size1(); i++){
        for(auto j = 0; j < mat.size2(); j++){
            mat(i,j) = dist(e2);
        }
    }
    rank = 0;
    if (rank == 0){
        for(auto i = 0; i < mat.size1(); i++){
            for(auto j = 0; j < mat.size2(); j++){
//                handle(mat, result, i, j);
            }
        }
    }


    std::string serial_str;
    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string> > s(inserter);
    boost::archive::binary_oarchive oa(s);

    oa << result;

    // don't forget to flush the stream to finish writing into the buffer
    s.flush();

//    std::cout << serial_str
//                 ;
    std::cout << "\n";

    return 0;
}
