#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/serialization/vector.hpp>


#include <iostream>
#include <random>
#include <algorithm>
#include <vector>

#include <fstream>

#include <boost/archive/binary_oarchive.hpp>

#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <cassert>

#define MAT_HEIGHT 10000
#define MAT_WIDTH 10000

namespace mpi = boost::mpi;
namespace ublas = boost::numeric::ublas;

bool isCorrect(char x){
    if (x == 'a' || x == 'u' || x == 'i' || x == 'e' || x == 'o' ||
        x == 'A' || x == 'U' || x == 'I' || x == 'E' || x == 'O'){
        return true;
    } else {
        return false;
    }
}

class OpVector{
    friend class boost::serialization::access;

    inline void generatevals(){
        length = std::max(abs(posx - startx), abs(posy-starty));
        vx = (posx-startx)/length;
        vy = (posy-starty)/length;
        assert(vx == 0 || vx == 1 || vx == -1);
        assert(vy == 0 || vy == 1 || vy == -1);
    }

    template<class Archive>  void serialize(Archive & ar, const unsigned int version)
    {
       ar & startx;
       ar & starty;
       ar & posx;
       ar & posy;
       generatevals();
    }
public:
    int length;
    int vx,vy; // vector to move
    int startx, starty;
    int posx, posy;

    OpVector(){}
    OpVector(int startx, int starty, int posx, int posy): startx(startx), starty(starty), posx(posx), posy(posy){
        generatevals();
    }

    bool canAdvance(ublas::matrix<char> &mat){
        return mat.size1() > nextX() && mat.size2() > nextY();
    }

    inline int nextX(){
        return posx + vx;
    }

    inline int nextY(){
        return posy + vy;
    }


    inline bool isCurrentCorrect(ublas::matrix<char> &mat){
        return  posx >= 0 &&
                posy >= 0 &&
                posx < mat.size1() &&
                posy < mat.size2() && isCorrect(mat(posx, posy));
    }

    inline bool isNextCorrect(ublas::matrix<char> &mat){
        return canAdvance(mat) && isCorrect(mat(nextX(), nextY()));
    }

    void advance() {
        posx += vx;
        posy += vy;
        length += 1;
    }

    boost::basic_format<char> toString() const {
        return boost::format("[%1%, %2%, %3%, %4% | %5%, %6%, %7% ]") % startx % starty % posx % posy % length % vx % vy;
    }
};

std::ostream &operator<<(std::ostream &os, OpVector const &m) {
    return os << m.toString();
}

bool recurse(boost::numeric::ublas::matrix<char> &mat, std::vector<OpVector> &result, OpVector &data){
    if (data.isCurrentCorrect(mat)){
        return false;
    }
    while (data.isNextCorrect(mat)) {
        data.advance();
    }

    if (data.length >= 5){
        result.push_back(data);
        // clean accepted to reduce duplicates
        int x = data.startx;
        int y = data.starty;
        for (int i=0; i <= data.length; i++){
            mat(x,y) = -1;
            x += data.vx;
            y += data.vy;
        }

    }

    return true;
}

class PositionHelper {
    int cellWidth = 2;
    int cellHeight = 2;
    int cellsInWidth;
    ublas::matrix<char> &mat;

public:
    PositionHelper(int cellWidth, int cellHeight, ublas::matrix<char> &mat) : cellWidth(cellWidth), cellHeight(cellHeight), mat(mat) {
        cellsInWidth = cellsInWidthCalc();
    }

    void processCell(ublas::matrix<char> &mat, std::vector<OpVector> &result, int cellNo) {
        int x = firstXInCell(cellNo);
        int y = firstYInCell(cellNo);
        std::cout << cellNo << " " << x << " " << y << std::endl;

        for(int i=x; i < std::min(x + cellWidth, (int)mat.size1()) ; i++){
            for(int j=y; j < std::min(y + cellHeight, (int)mat.size2()); j++){
                if (isCorrect(mat(i,j))){
                    OpVector opv1(i, j, i+1, j);
                    OpVector opv2(i, j, i+1, j+1);
                    OpVector opv3(i, j, i,   j+1);
                    OpVector opv4(i, j, i-1, j+1);

                    recurse(mat, result, opv1);
                    recurse(mat, result, opv2);
                    recurse(mat, result, opv3);
                    recurse(mat, result, opv4);
                }
            }
        }
    }

    inline int cellsInWidthCalc(){
        int result = mat.size1() / cellWidth;
        if (mat.size1() % cellWidth > 0){
            result += 1;
        }
        return result;
    }

    int numCells(){
        int wCells = cellsInWidthCalc();
        int hCells = mat.size2() / cellHeight;
        if (mat.size2() % cellHeight > 0){
            hCells +=1;
        }
        return wCells * hCells;
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
    ublas::matrix<char> mat(3,2);
    PositionHelper h3(2, 2, mat);

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

    std::uniform_int_distribution<char> dist(32,125);

    e2.seed(10);

    std::vector<OpVector> result;

    mpi::environment env(argc, argv);
    mpi::communicator world;
    std::cout << "I am process " << world.rank() << " of " << world.size()
              << "." << std::endl;

    boost::numeric::ublas::matrix<char> mat(MAT_WIDTH, MAT_HEIGHT);

    PositionHelper ph(10000, 1000, mat);
    int ws = world.size();
    int numCells = ph.numCells();
    assert(ws > 1);

    if (world.rank() == 0){
        std::vector<OpVector> results;
        for (int i =0; i < numCells; i++){
            int target = (i % (ws - 1)) + 1;
            world.send(target, 1, i);
        }

        for (int i = 1; i < ws; i++){
            world.send(i, 1, -1);
        }

//        for (int i =1 ; i< ws; i++){
//            std::vector<OpVector> results_tmp;
//            world.recv(i, 2, results_tmp);
//            results.reserve(results.size() + results_tmp.size());
//            results.insert(results.end(), results_tmp.begin(), results_tmp.end());
//        }
        std::cout << "num results: " << results.size() << "\n";
    } else {

        for(auto i = 0; i < mat.size1(); i++){
            for(auto j = 0; j < mat.size2(); j++){
                mat(i,j) = dist(e2);
            }
        }

        std::vector<OpVector> results;
        while (true){
            int msg;
            world.recv(0, 1, msg);
            if (msg == -1){
                break;
            }
            ph.processCell(mat, results, msg);
        }
        std::cout << "Sending n: " << results.size() << "\n";
//        world.send(0, 2, results);

    }
    std::cout << "\n";

    return 0;
}
