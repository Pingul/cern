#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include "common.hh"
#include "timer.hh"

#include "accelerator.hh"
#include "settings.hh"
#include "hamiltonian.hh"
#include "particles.hh"

#include "ramp_program.hh"

#ifdef IO_TO_SAME_DIR
#define RESOURCE_DIR "."
#define OUTPUT_DIR "."
#else
#define RESOURCE_DIR "resources"
#define OUTPUT_DIR "calc"
#endif

namespace stron {



template <typename T>
class SimpleSynchrotron 
{
    struct FileNotFound : public std::runtime_error {
        FileNotFound(const std::string& fn) : std::runtime_error(fn) {}
    };

public:
    using Acc = Accelerator<T>;
    using Particles = ParticleCollection<T>;
    using ParticlesPtr = typename Particles::Ptr;
    using Collimat = typename Acc::Collimat;
    using ProgramPtr = typename Program<Acc>::Ptr;

    SimpleSynchrotron(const Acc acc) 
        : mAcc(acc), mParticles(nullptr), mProgram(nullptr) {}

    void addParticles(ParticlesPtr particles)
    {
        mParticles = particles;
        mCollHits.assign(mParticles->size(), -1);
        writeSingleDistribution(STARTDIST_FILE);
    }

    void createTurnFileHeader(std::string filePath, int turns) const
    {
        std::ofstream file(filePath.c_str());
        if (file.is_open()) {
            std::cout << "Saving turns data in '" << filePath << "'" << std::endl;
            file << mParticles->size() << "," << turns << std::endl; 
        } else 
            std::cerr << "could not write headers" << std::endl;
        file.close();
    }

    void writeMetaData(const std::string& filePath) const 
    {
        std::ofstream f(filePath.c_str());
        if (!f.is_open())
            throw FileNotFound(filePath);
    }

    void writeDistribution(std::string filePath) const 
    {
        std::ofstream file(filePath.c_str(), std::ios::app);

        if (!file.is_open()) 
            throw FileNotFound(filePath);

        // FILE:
        //      x, px, ∆energy,phase,h        p1
        //      ...                           p2
        for (size_t i = 0; i < mParticles->size(); ++i) {
            std::stringstream ss;
            const T de = mParticles->momentum[i];
            const T phase = mParticles->phase[i];
            ss << std::setprecision(16) << 
                de << "," << 
                phase << "," << 
                hamiltonian(mAcc, de, phase) << std::endl;
            file << ss.str();
        }
        file.close();
    }

    void writeSingleDistribution(std::string filePath) const
    {
        createTurnFileHeader(filePath, 1);
        writeDistribution(filePath);
    }

    void writeCollHits(std::string filePath) const
    {
        std::cout << "Saving coll hits to '" << filePath << "'" << std::endl;
        std::ofstream file(filePath.c_str());
        if (!file.is_open())
            throw std::runtime_error("could not open file");
        else if (mCollHits.empty())
            throw std::runtime_error("could not save coll hits -- no data was recorded");

        // FILE:
        //     id / turn_lost / phase_lost / energy_lost
        file << mParticles->size() << std::endl;
        //file << mAcc.coll_top << ", " << mAcc.coll_bot << std::endl;
        file << 0.0 << ", " << 0.0 << std::endl;
        for (size_t i = 0; i < mParticles->size(); ++i) {
            if (mCollHits[i] == -1) continue;
            std::stringstream ss;
            file << i << ", " << 
                mCollHits[i] << ", " << std::setprecision(16) << 
                mParticles->phase[i] << ", " << 
                mParticles->momentum[i] << std::endl;
        }
    }

    //
    // SIMULATION
    //
    struct ParticleStats { T emax, emin, phmax, phmin; int pleft; };
    ParticleStats getStats() const
    {
        T emax, emin, phmax, phmin;
        emin = phmin = std::numeric_limits<T>::max();
        emax = phmax = -emin;
        int pleft = 0;
        for (size_t i = 0; i < mParticles->size(); ++i) {
            if (mCollHits[i] < 0) {
                ++pleft;
                emax = mParticles->momentum[i] > emax ? mParticles->momentum[i] : emax;
                emin = mParticles->momentum[i] < emin ? mParticles->momentum[i] : emin;
                phmax = mParticles->phase[i] > phmax ? mParticles->phase[i] : phmax;
                phmin = mParticles->phase[i] < phmin ? mParticles->phase[i] : phmin;
            }
        }
        return ParticleStats{emax, emin, phmax, phmin, pleft};
    }

    void simulateTurn(int stepID) 
    {
        CalcOp op(stepID, mAcc, *mParticles, mCollHits);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, mParticles->size()), op);
    }
    
    //void simulateTurns(int n, std::string filePath = "", int saveFreq = 1)
    void simulateTurns(ProgramPtr program, std::string filePath = "", int saveFreq = 1)
    {
        mTurns = program->steps();
        program->setup();

        if (filePath.empty())
            std::cout << "Will not save particle path data" << std::endl;
        else  {
            createTurnFileHeader(filePath, 1 + mTurns/saveFreq); // +1 as we are saving the first starting configuration as well
            writeDistribution(filePath);
        }

        std::cout << "Tracking " << mParticles->size() << " particles for " << mTurns << " turns" << std::endl;
        std::cout << "Starting simulation..." << std::endl;

        common::SilentTimer timer;

        timer.start();
        for (int i = 1; i < mTurns; ++i) {
            simulateTurn(i);
            program->step();
            mAcc.recalc();

            // reduce the memory footprint a little if it's not necessary
            if (!filePath.empty() && (i + 1) % saveFreq == 0)
                writeDistribution(filePath);

            const int d = mTurns/10;
            if (d > 0 && i % d == 0) {
                int percent = 10*i/d;
                ParticleStats stats = getStats();
                std::cout << "\t" << std::setw(9) << i << " of " << mTurns << " turns (" << percent << "%).";

                int pleft = (100*stats.pleft)/mParticles->size();
                std::cout << " #=" << pleft << "%. φ=[" << std::setprecision(3) << stats.phmin << ", " << stats.phmax << "]" << std::endl;
            }
        }

        double ms = timer.stop();

        auto dur = common::MillisecondsToElapsedTime(unsigned(ms));
        std::cout << "Finished in ";
        if (dur.d > 0) std::cout << dur.d << "d ";
        if (dur.h > 0) std::cout << dur.h << "h ";
        if (dur.m > 0) std::cout << dur.m << "m ";
        if (dur.s > 0) std::cout << dur.s << "s ";
        std::cout << dur.ms << "ms" << std::endl;

        writeSingleDistribution(ENDDIST_FILE);
    }

    void runLossmap(ProgramPtr program)
    {
    
        simulateTurns(program, stron::PATH_FILE, 11245); 
        writeCollHits(COLL_FILE);
    
        int ilast = 0;
        for (size_t i = 1; i < mCollHits.size(); ++i) {
            if (mCollHits[i] > mCollHits[ilast]) 
                ilast = i;
        }
    
        int tlast = mCollHits[ilast];
        if (tlast == -1)
            std::cout << "No losses" << std::endl;
        else {
            std::cout 
                << "Latest hit:\n\tparticle " << ilast << ", turn " << tlast 
                << "(approx. after " << std::setprecision(3) << (double(tlast)/mAcc.f_rev) << " s)\n";
        }
    
    }

    Acc& getAcc() { return mAcc; }

private:
    struct CalcOp
    {
        CalcOp(int stepID, const Acc& acc, Particles& particles, std::vector<int>& collHits)
            : mStepID(stepID), mAcc(acc), mPart(particles), mCollHits(collHits)
        {}

        bool particleInactive(size_t index) const
        {
            return !mCollHits.empty() && mCollHits[index] != -1;
        }

        bool particleCollided(size_t index) const
        {
            bool collided = false;
            for (auto& coll : mAcc.collimators) {
                switch (coll.type) {
                    case Collimat::Type::TCP_IR3: {
                        const T& momentum = mPart.momentum[index];
                        const T dispersion = -2.07e3;
                        const T cut = (std::abs(coll.left) + std::abs(coll.right))/(2.0*dispersion);
                        const T mcut = mAcc.E()*cut;
                        collided |= momentum < mcut || momentum > -mcut;
                        break;
                    } case Collimat::Type::TCPc_IR7: {
                        break;
                    }
                }
            }
            return collided;
        }

        bool outsideBucket(size_t index) const 
        {
            return mPart.phase[index] < -0.1;
        }

        void longitudinalStep(size_t index) const
        {
            using std::sin;
            using cnst::pi;

            T deltaRef = mAcc.E() - mAcc.E_prev();

            int Ns = 100;
            // Particles outside of the bucket does not need to be tracked as carefully
            if (outsideBucket(index)) Ns = 1;

            T& momentum = mPart.momentum[index];
            T& phase = mPart.phase[index];
            for (int i = 0; i < Ns; ++i) {
                momentum += (mAcc.V_rf*(sin(phase)) - deltaRef)/T(Ns);
                auto p = mAcc.calcParticleProp(momentum, 0.0);
                phase -= T(2)*pi*mAcc.h_rf*p.eta/(T(Ns)*p.b2*mAcc.E())*momentum;
            }
        }

        void transverseStep(size_t index) const
        {
        }

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t n = range.begin(), N = range.end(); n < N; ++n) {
                if (particleInactive(n)) 
                    continue;

                longitudinalStep(n);
                transverseStep(n);

                if (particleCollided(n)) mCollHits[n] = mStepID;
            }
        }
    private:
        int mStepID;
        const Acc& mAcc;
        Particles& mPart;
        std::vector<int>& mCollHits;
    }; // CalcOp

    int mTurns;
    Acc mAcc;
    ParticlesPtr mParticles;
    std::vector<int> mCollHits; // Turn hitting collimator
    ProgramPtr mProgram;
};

} // namespace stron
