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
        : mAcc(acc), mParticles(nullptr), mProgram(nullptr) 
    {
        // Erasing previous files
        { std::ofstream f(PATH_FILE); }
        { std::ofstream f(STARTDIST_FILE); }
        { std::ofstream f(ENDDIST_FILE); }
    }

    void addParticles(ParticlesPtr particles)
    {
        mParticles = particles;
        for (const auto& c : mAcc.collimators)
            mCollHits.emplace_back(CollimatorHits(c, *mParticles));
        writeDistribution(STARTDIST_FILE);
    }

    void writeDistribution(std::string filePath) const 
    {
        std::ofstream file(filePath.c_str(), std::ios::app);

        if (!file.is_open()) 
            throw FileNotFound(filePath);

        
        for (const Collimat& c : mAcc.collimators) {
            if (c.type == Collimat::TCP_IR3) {
                for (size_t i = 0; i < mParticles->size(); ++i) {
                    auto p = mAcc.calcParticleProp(mParticles->momentum[i], 0.0);
                    auto xc = mParticles->xBeta(i, c.alpha, c.beta, p.b, p.g);
                    file << i << ","
                         << xc.x << ","
                         << xc.px << ","
                         //<< mParticles->x[i] << ","
                         //<< mParticles->px[i] << ","
                         << mParticles->momentum[i] << ","
                         << mParticles->phase[i] << ","
                         << hamiltonian(mAcc, mParticles->momentum[i], mParticles->phase[i]) 
                         << std::endl;
                }
            }
        }
    }

    void writeCollHits(std::string filePath) const
    {
        std::cout << "** Collimator stats **" << std::endl;
        for (const auto& ch : mCollHits) {
            if (ch.collimat.type == Collimat::Type::TCP_IR3) {
                ch.printStats();
                ch.write(filePath);
            }
        }
        std::cout << "**" << std::endl;
    }

    //
    // SIMULATION
    //
    void simulateTurn(int turn) 
    {
        CalcOp op(turn, mAcc, *mParticles, mCollHits);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, mParticles->size()), op);
    }
    
    //void simulateTurns(int n, std::string filePath = "", int saveFreq = 1)
    void simulateTurns(ProgramPtr program, std::string filePath = "", int saveFreq = 1)
    {
        MetaData meta;
        meta.nbrp = mParticles->size();
        meta.turns = program->steps();
        meta.saveFreq = saveFreq;

        if (filePath.empty())
            std::cout << "Will not save particle path data" << std::endl;
        else  {
            writeDistribution(filePath);
            ++meta.savedTurns;
        }

        std::cout << "Tracking " << meta.nbrp << " particles for " << meta.turns << " turns" << std::endl;
        std::cout << "Starting simulation..." << std::endl;

        common::SilentTimer timer;
        timer.start();

        program->setup();
        for (int i = 0; i < meta.turns - 1; ++i) {
            simulateTurn(i);
            program->step();
            mAcc.recalc();

            if (!filePath.empty() && (i + 1) % saveFreq == 0) {
                writeDistribution(filePath);
                ++meta.savedTurns;
            }

            const int d = meta.turns/10;
            if (d > 0 && i % d == 0) {
                int percent = 10*i/d;
                SimStats stats = getStats();
                std::cout << "\t" << std::setw(9) << i << " of " << meta.turns << " turns (" << percent << "%).";
                if (percent < 10) std::cout << " ";

                int pleft = (100*stats.pleft)/mParticles->size();
                std::cout << " #=" << pleft << "%. Eref=" << mAcc.E() << ", φ=[" << std::setprecision(3) << stats.phmin << ", " << stats.phmax << "]" << std::endl;
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

        writeDistribution(ENDDIST_FILE);
        writeCollHits(COLL_FILE);
        meta.write(META_FILE);
    }

    void runLossmap(ProgramPtr program)
    {
        simulateTurns(program, stron::PATH_FILE, 11245); 
        std::cout << "** Collimator stats **" << std::endl;
        for (const auto& ch : mCollHits) 
            ch.printStats();
        std::cout << "**" << std::endl;
    }

    Acc& getAcc() { return mAcc; }

private:
    struct CollimatorHits 
    {
        using Collimat = typename Acc::Collimat;
        CollimatorHits(const Collimat& c, const Particles& p) 
            : collimat(c), mParticles(p), mTurns(p.size(), -1), mX(p.size()) {}

        bool isHit(int pId) { return mTurns.at(pId) == -1; }
        void registerHit(int pId, int turn, T xpos) { mTurns.at(pId) = turn; mX.at(pId) = xpos; }

        void write(const std::string& filePath) const
        {        
            std::cout << "Saving collimator (" << collimat.type << ") hits to '"  << filePath << "'" << std::endl;
            std::ofstream file(filePath.c_str());
            if (!file.is_open())
                throw FileNotFound(filePath);

            file << mParticles.size() << std::endl;
            //file << mAcc.coll_top << ", " << mAcc.coll_bot << std::endl;
            file << 0.0 << ", " << 0.0 << std::endl;
            for (size_t i = 0; i < mParticles.size(); ++i) {
                if (mTurns[i] == -1) continue;
                file << i << ", " 
                    << mTurns[i] << "," << std::setprecision(16) 
                    << mParticles.phase[i] << "," 
                    << mParticles.momentum[i] << ","
                    << mX[i] << std::endl;
            }
        }

        void printStats() const 
        {
            int n = 0;
            int last = 0;
            T maxX = 0.0;
            for (size_t i = 0; i < mTurns.size(); ++i) {
                if (mTurns[i] == -1) continue;
                last = (mTurns[i] > last) ? mTurns[i] : last;
                maxX = (mX[i] > maxX) ? mX[i] : maxX;
                n++;
            }

            std::cout << "For collimator " << collimat.type << ":" << std::endl;
            std::cout << "\tHits : " << n << std::endl;
            std::cout << "\tLast : " << last << " (" << std::setprecision(4) << (last/cnst::s_to_turn) << " s)" << std::endl;
            std::cout << "\tMax x: " << maxX << std::endl;
        }

        const Collimat& collimat;
    private:
        const Particles& mParticles;
        std::vector<int> mTurns;
        std::vector<T> mX;
    }; // CollimatorHits

    struct CalcOp
    {
        CalcOp(int turn, const Acc& acc, Particles& particles, std::vector<CollimatorHits>& chits)
            : mTurn(turn), mAcc(acc), mPart(particles), mCHits(chits)
        {}

        bool particleCollided(size_t index) const
        {
            bool collided = false;
            for (auto& ch : mCHits) {
                auto& coll = ch.collimat;
                switch (coll.type) {
                    case Collimat::Type::TCP_IR3: {
                        const T& dp = mPart.momentum[index];
                        auto p = mAcc.calcParticleProp(dp, 0.0);
                        const T xd = dp/mAcc.E()*coll.dispersion;
                        const T xb = mPart.xBeta(index, coll.alpha, coll.beta, p.b, p.g).x;
                        const T x = xb + xd;
                        const T xcut = (std::abs(coll.left) + std::abs(coll.right))/2.0*1e-3;
                        bool cond = x < -xcut || x > xcut;
                        if (cond) ch.registerHit(index, mTurn, x);
                        collided |= cond;
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
            return mPart.phase[index] < -0.1 || mPart.phase[index] > 2.1*cnst::pi;
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
            using std::cos;
            using std::sin;
            using cnst::Qx;

            T& x = mPart.x[index];
            const T xc = x;
            T& px = mPart.px[index];
            x = cos(Qx)*xc + sin(Qx)*px;
            px = -sin(Qx)*xc + cos(Qx)*px;
        }

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t n = range.begin(), N = range.end(); n < N; ++n) {
                if (!mPart.isActive(n)) 
                    continue;

                longitudinalStep(n);
                transverseStep(n);

                if (particleCollided(n))
                    mPart.setActive(n, false);
            }
        }
    private:
        int mTurn;
        const Acc& mAcc;
        Particles& mPart;
        std::vector<CollimatorHits>& mCHits;
    }; // CalcOp

    struct MetaData
    {
        int nbrp{0};
        int turns{0};
        int saveFreq{0};
        int savedTurns{0};

        void write(const std::string& filePath) const 
        {
            std::ofstream f(filePath.c_str());
            if (!f.is_open())
                throw FileNotFound(filePath);

            f << "{" << std::endl;
            f << "\t\"nbrp\": " << nbrp << "," << std::endl;
            f << "\t\"turns\": " << turns << "," << std::endl;;
            f << "\t\"saveFreq\": " << saveFreq << "," << std::endl;;
            f << "\t\"savedTurns\": " << savedTurns << std::endl;;
            f << "}" << std::endl;
        }
    };    

    struct SimStats { T emax, emin, phmax, phmin; int pleft; };
    SimStats getStats() const
    {
        T emax, emin, phmax, phmin;
        emin = phmin = std::numeric_limits<T>::max();
        emax = phmax = -emin;
        int pleft = 0;
        for (size_t i = 0; i < mParticles->size(); ++i) {
            if (!mParticles->isActive(i)) 
                continue;
            ++pleft;
            emax = mParticles->momentum[i] > emax ? mParticles->momentum[i] : emax;
            emin = mParticles->momentum[i] < emin ? mParticles->momentum[i] : emin;
            phmax = mParticles->phase[i] > phmax ? mParticles->phase[i] : phmax;
            phmin = mParticles->phase[i] < phmin ? mParticles->phase[i] : phmin;
            
        }
        return SimStats{emax, emin, phmax, phmin, pleft};
    }

    Acc mAcc;
    ParticlesPtr mParticles;
    std::vector<CollimatorHits> mCollHits;
    ProgramPtr mProgram;

};

} // namespace stron
