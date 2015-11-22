#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <boost/serialization/vector.hpp>

#include <iostream>
#include <vector>
#include <fstream>

#include <boost/iostreams/stream.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace mpi = boost::mpi;

class SerExample{
    friend class boost::serialization::access;

    template<class Archive>  void serialize(Archive & ar, const unsigned int /*version*/)
    {
        ar & val;
    }

public:
    std::string val;
    SerExample() {}
    SerExample(std::string val): val(val) {}
};

std::ostream &operator<<(std::ostream &os, SerExample const &m) {
    return os << "[" << m.val << "]";
}


int main(int argc, char **argv){
    mpi::environment env(argc, argv);
    mpi::communicator world;
    std::cout << "I am process " << world.rank() << " of " << world.size()
              << "." << std::endl;

    if (world.rank() == 0){
        for (int i =1; i < world.size(); i++){
            world.send(i, 1, std::string("ehlo ") + std::to_string(i));
        }

        for (int i = 1; i< world.size(); i++){
            std::vector<SerExample> tmp;
            world.recv(i, 1, tmp);

            std::cout << "received vector: " << std::endl;
            for (auto const & value: tmp){
                std::cout << value << value;
            }
        }

    } else {
        std::string tmp;
        world.recv(0, 1, tmp);

        std::cout << "Received " << tmp << "\n";
        std::vector<SerExample> response;
        response.push_back(SerExample("test"));
        response.push_back(SerExample(tmp));
        world.send(0, 1, response);
    }    
    return 0;
}
