#include <boost/mpi/environment.hpp>``
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
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>

#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <cassert>

#define MAT_HEIGHT 5
#define MAT_WIDTH 5

namespace mpi = boost::mpi;
namespace ublas = boost::numeric::ublas;

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

    void advance() {
        posx += vx;
        posy += vy;
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


inline void processResult(std::vector<OpVector> &result, OpVector &data){

    if (data.length >= 2){
        result.push_back(data);
    }
}

bool recurse(boost::numeric::ublas::matrix<char> &mat, std::vector<OpVector> &result, OpVector &data){
    while (data.canAdvance(mat)) {
        if (isCorrect(mat(data.nextX(), data.nextY()))){
           data.advance();
        } else {
           std::cout << data << "\n";
           break;
        }
    }
    if (data.length >= 2){
        result.push_back(data);
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

    std::uniform_int_distribution<char> dist(48,49);
    boost::numeric::ublas::matrix<char> mat(MAT_WIDTH, MAT_HEIGHT);

    e2.seed(10);

    std::vector<OpVector> result;

    mpi::environment env(argc, argv);
    mpi::communicator world;
    std::cout << "I am process " << world.rank() << " of " << world.size()
              << "." << std::endl;

    for(auto i = 0; i < mat.size1(); i++){
        for(auto j = 0; j < mat.size2(); j++){
            mat(i,j) = dist(e2);
        }
    }

    PositionHelper ph(10, 10, mat);
    int ws = world.size();
    int numCells = ph.numCells();
//    assert(ws > 1);
    std::vector<OpVector> results;
    if (world.rank() == 0){
        for (int i =0; i < numCells; i++){
//            int target = (i % (ws - 1)) + 1;
//            int target =1;
//            world.send(target, 1, i);
            ph.processCell(mat, results, i);
            std::cout << "aaa " << " " << i << "\n";
        }
        for (int i = 1; i < ws; i++){
            world.send(i, 1, -1);
        }

    } else {
        std::cout << "AsfasfAFSAFS\n";
        std::vector<OpVector> results;
        while (true){
            int msg;
            world.recv(0, 1, msg);
            ph.processCell(mat, results, msg);
//            std::cout << "dupa " << msg << std::endl;
            if (msg == -1){
                break;
            }

        }
//        for(int i=0; i< results.size(); i++){
//            std::cout << results[i] << std::endl;
//        }
    }
    std::cout << mat;
    for(int i=0; i< results.size(); i++){
        std::cout << results[i] << std::endl;
    }


//    std::string serial_str;
//    boost::iostreams::back_insert_device<std::string> inserter(serial_str);
//    boost::iostreams::stream<boost::iostreams::back_insert_device<std::string> > s(inserter);
//    boost::archive::binary_oarchive oa(s);

//    oa << result;

//    // don't forget to flush the stream to finish writing into the buffer
//    s.flush();

////    std::cout << serial_str
//                 ;
    std::cout << "\n";

    return 0;
}
