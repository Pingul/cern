#include "particles.hh"
#include "hamiltonian.hh"

#include <fstream>
#include <sstream>
#include <iomanip>

namespace stron {


//template <typename T>
//void ParticleCollection<T>::write(std::string filePath) const
//{
    //std::ofstream file(filePath.c_str(), std::ios::app);
    //if (!file.is_open())
        //throw std::runtime_error("Could not write particles");
    
    //for (size_t i = 0; i < size(); ++i) {
        //std::stringstream ss;
        //ss << std::setprecision(16) << 
            //mEnergy[i] << "," << 
            //mPhase[i] << "," << 

            // What to do about this one?
            //hamiltonian(mAcc, mEnergy[i], mPhase[i]) << std::endl;
        //file << ss.str();
    //}
//}


}


template class stron::ParticleCollection<double>;
