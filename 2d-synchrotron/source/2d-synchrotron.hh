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
#include "common.hh"
#include "timer.hh"

#include "accelerator.hh"
#include "settings.hh"
#include "hamiltonian.hh"
#include "particles.hh"

#ifdef IO_TO_SAME_DIR
#define RESOURCE_DIR "."
#define OUTPUT_DIR "."
#else
#define RESOURCE_DIR "resources"
#define OUTPUT_DIR "calc"
#endif

namespace stron {

enum RAMP_TYPE
{
    NO_RAMP,
    LHC_RAMP,
    EXTENDED_LHC_RAMP,
    SEMI_AGGRESSIVE_RAMP,
    AGGRESSIVE_RAMP,
    EXTERNAL_RAMP,
};


template <typename T>
inline void readRampFile(int steps, std::string filePath, std::vector<T>& E_ramp)
{
    using common::skip;

    E_ramp.reserve(steps);
    T data;
    std::cout << "Reading '" << filePath << "'..." << std::endl;
    std::ifstream ramp_file(filePath);
    for (int i = 0; i < steps; ++i) {
        ramp_file >> skip >> data;
        E_ramp.push_back(data*1e6);
    }
}

template <typename T>
inline void readRamp(int steps, std::vector<T>& E_ramp, RAMP_TYPE type)
{
    E_ramp.reserve(steps);
    switch (type) {
        default:
        case LHC_RAMP: {
            readRampFile(steps, LHC_RAMP_FILE, E_ramp);
            break;
        }
        case EXTENDED_LHC_RAMP: {
            int extra_steps = 150; // A particle close to the separatrix takes ~600 turns to go around
            readRamp(steps - extra_steps, E_ramp, LHC_RAMP);
            std::vector<T> extension(extra_steps, 450e9); 
            E_ramp.insert(E_ramp.begin(), extension.begin(), extension.end());

            std::cout << "--- VERIFY EXTENSION DATA ----" << std::endl;
            for (int i = 0; i < 50; ++i) {
                for (int j = 0; j < 8; ++j) {
                    int index = i*8 + j;
                    std::cout << std::setw(20) << std::setprecision(16) << E_ramp[index];
                }
                std::cout << std::endl;
            }
            std::cout << "------------------------------" << std::endl;
            break;
        }
        case EXTERNAL_RAMP: 
            readRampFile(steps, EXTERNAL_RAMP_FILE, E_ramp);
            break;
        case AGGRESSIVE_RAMP:
            for (int i = 0; i < steps; ++i) E_ramp.push_back(450e9 + i*1e7);
            break;
        case SEMI_AGGRESSIVE_RAMP:
            for (int i = 0; i < steps; ++i) E_ramp.push_back(450e9 + i*3e6);
            break;
    }
}

template <typename T>
inline void readCollimators(int steps, std::vector<std::pair<T, T>>& collimators, const std::vector<T>& E_ramp)
{
    using common::skip;

    std::cout << "Reading collimators from '" << COLL_MOTOR_FILE << "'" << std::endl;
    std::ifstream coll_file(COLL_MOTOR_FILE);
    if (!coll_file.is_open())
        throw std::runtime_error("could not open file");

    for (int i = 0; i < steps; ++i) {
        T bot_coll, top_coll;
        coll_file >> skip >> bot_coll >> top_coll;
        collimators.push_back(std::make_pair(bot_coll*E_ramp[i], top_coll*E_ramp[i]));
    }
}

template <typename T>
struct SimpleSynchrotron 
{
    using Acc = Accelerator<T>;
    using PColl = ParticleCollection<T>;

    SimpleSynchrotron(size_t n, RAMP_TYPE type)
        : mAcc(Acc::getLHC()), mParticles(n), mType(type)
    {
        initStorage(n);

        std::random_device rdev;
        std::mt19937 generator(rdev());

        std::normal_distribution<> e_dist(0, 0.20e9);
        std::normal_distribution<> ph_dist(cnst::pi, cnst::pi/4);
        for (size_t i = 0; i < n; ++i) {
            const T deltaE = e_dist(generator);
            const T phase = ph_dist(generator);
            mParticles.momentum[i] = deltaE;
            mParticles.phase[i] = phase;
        }
        writeSingleDistribution(STARTDIST_FILE);
    }


    // For parameter passing in the next constructors
    // These are all different distributions to choose from
    struct AroundSeparatrix {};
    struct ActionValues {};
    struct LinearDecay {};
    struct ExponentialDecay {};
    
    struct SixTrackTest {};

    SimpleSynchrotron(int n, RAMP_TYPE type, AroundSeparatrix)
        : mAcc(Acc::getLHC()), mParticles(n), mType(type)
    {
        initStorage(n);

        std::random_device rdev;
        std::mt19937 generator(rdev());

        std::uniform_real_distribution<> e_dist(-0.5e9, 0.5e9);
        std::uniform_real_distribution<> ph_dist(0, 2*cnst::pi);
        int count = 0;
        const T sep = separatrix<T>(mAcc);
        while (count < n) {
            const T deltaE = e_dist(generator);
            const T phase = ph_dist(generator);
            const T H = hamiltonian(mAcc, deltaE, phase);
            if ((sep - 1e6) < H && H < (sep + 1e6)) {
                mParticles.momentum[count] = deltaE;
                mParticles.phase[count] = phase;
                count++;
            }
        }
        writeSingleDistribution(STARTDIST_FILE);
    }

    SimpleSynchrotron(int N, RAMP_TYPE type, ActionValues)
        : mAcc(Acc::getLHC()), mParticles(N), mType(type)
    {
        const T sep = separatrix(mAcc);
        //std::vector<T> d_actions = {-1e4, -9e3, -8e3, -7e3, -6e3, -5e3, -4e3, -3e3, -2e3, -1e3, -100};
        std::vector<T> d_actions;
        
        int outside = 80;
        double maxdE = 1.6e9;
        double de = maxdE/double(outside);
        for (int i = 1; i <= outside; ++i) {
            d_actions.push_back(hamiltonian(mAcc, de*double(i), cnst::pi));
        }
        
        int inside = 20;
        for (int i = 1; i <= inside; ++i) {
            d_actions.push_back(T(-8000 + 500*i));
        }

        const int n = N/(2*d_actions.size());
        initStorage(2*n*d_actions.size());

        int count = 0;
        std::random_device rdev;
        std::mt19937 generator(rdev());
        for (T action : d_actions) {
            action += sep;
            std::uniform_real_distribution<> dist(0.0, 2*cnst::pi);
            for (T sign : std::vector<T>({-1.0, 1.0})) {
                for (int i = 0; i < n; ++i) {
                    const T phase = dist(generator);
                    const T energy = levelCurve(mAcc, phase, action, sign);
                    if (std::isnan(energy)) { --i; continue; }
                    mParticles.momentum[count] = energy;
                    mParticles.phase[count] = phase;
                    ++count;
                }
            }
        }
        std::cout << "Tried to create " << N << " particles, only initialized " << size() << std::endl;
        writeSingleDistribution(STARTDIST_FILE);
    }

    SimpleSynchrotron(int n, RAMP_TYPE type, LinearDecay)
        : mAcc(Acc::getLHC()), mParticles(n), mType(type)
    {
        initStorage(n);
        const T sep = separatrix(mAcc);
        std::random_device rdev;
        std::mt19937 generator(rdev());

        // Distribution as (relative to separatrix)
        //  -15k        -10k         +5k
        //    | constant  | linear dec.|
        std::vector<T> Hs{sep - 15e3, sep - 10e3, sep + 5e3, sep + 2.7e7};
        std::vector<T> prob{100.0, 100.0, 0.1, 0.1};
        std::piecewise_linear_distribution<> H_dist(Hs.begin(), Hs.end(), prob.begin());

        std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);

        for (int i = 0; i < n; ++i) {
            const T phase = uni_dist(generator);
            const T sign = uni_dist(generator) < cnst::pi ? 1.0 : -1.0;
            const T action = H_dist(generator);
            const T energy = levelCurve(mAcc, phase, action, sign);
            if (std::isnan(energy)) { --i; continue; }
            mParticles.momentum[i] = energy;
            mParticles.phase[i] = phase;
        }
        std::cout << "Initialized " << size() << " particles" << std::endl;
        writeSingleDistribution(STARTDIST_FILE);
    }

    SimpleSynchrotron(int n, RAMP_TYPE type, ExponentialDecay)
        : mAcc(Acc::getLHC()), mParticles(n), mType(type)
    {
        initStorage(n);

        std::random_device rdev;
        std::mt19937 generator(rdev());

        std::exponential_distribution<> H_dist(0.0003);
        std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);

        const T H_low = separatrix(mAcc) - 15000;
        for (int i = 0; i < n; ++i) {
            const T phase = uni_dist(generator);
            const T sign = uni_dist(generator) < cnst::pi ? 1.0 : -1.0;
            const T action = H_low + H_dist(generator);
            const T energy = levelCurve(mAcc, phase, action, sign);
            if (std::isnan(energy)) { --i; continue; }
            mParticles.momentum[i] = energy;
            mParticles.phase[i] = phase;
        }
        std::cout << "Initialized " << size() << " particles" << std::endl;
        writeSingleDistribution(STARTDIST_FILE);
    }

    SimpleSynchrotron(SixTrackTest, T energy)
        : mAcc(Acc::getLHC()), mParticles(1), mType(LHC_RAMP)
    {
        mParticles.momentum[0] = energy;
        mParticles.phase[0] = cnst::pi;
    }

    void initStorage(int n)
    {
        mCollHits.assign(n, -1);
    }


    //
    // FILE IO
    //

    void readDistribution(std::string filePath)
    {
        std::cout << "Reading distribution from '" << filePath << "'" << std::endl;
        std::ifstream file(filePath.c_str());
        if (!file.is_open())
            throw std::runtime_error("could not open file");

        using common::skip;

        int n;
        file >> n >> skip;
        initStorage(n);
    
        for (int i = 0; i < n; ++i) {
            T de, phase;
            file >> de >> skip<char> >> phase >> skip; // we need to consume end of line
            mParticles.momentum[i] = de;
            mParticles.phase[i] = phase;
        }
        std::cout << "Read " << size() << " particles" << std::endl;
    }

    void createTurnFileHeader(std::string filePath, int turns) const
    {
        std::ofstream file(filePath.c_str());
        if (file.is_open()) {
            std::cout << "Saving turns data in '" << filePath << "'" << std::endl;
            file << size() << "," << turns << std::endl; 
        } else 
            std::cerr << "could not write headers" << std::endl;
        file.close();
    }

    void writeDistribution(std::string filePath) const 
    {
        std::ofstream file(filePath.c_str(), std::ios::app);

        if (!file.is_open()) 
            throw std::runtime_error("could not save particle coordinates");

        // FILE:
        //      ∆energy,phase,h        p1
        //      ...                    p2
        for (size_t i = 0; i < size(); ++i) {
            std::stringstream ss;
            const T de = mParticles.momentum[i];
            const T phase = mParticles.phase[i];
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
        file << size() << std::endl;
        file << mAcc.coll_top << ", " << mAcc.coll_bot << std::endl;
        for (size_t i = 0; i < size(); ++i) {
            if (mCollHits[i] == -1) continue;
            std::stringstream ss;
            file << i << ", " << 
                mCollHits[i] << ", " << std::setprecision(16) << 
                mParticles.phase[i] << ", " << 
                mParticles.momentum[i] << std::endl;
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
        for (size_t i = 0; i < size(); ++i) {
            if (mCollHits.size() == size() && mCollHits[i] < 0) {
                ++pleft;
                emax = mParticles.momentum[i] > emax ? mParticles.momentum[i] : emax;
                emin = mParticles.momentum[i] < emin ? mParticles.momentum[i] : emin;
                phmax = mParticles.phase[i] > phmax ? mParticles.phase[i] : phmax;
                phmin = mParticles.phase[i] < phmin ? mParticles.phase[i] : phmin;
            }
        }
        return ParticleStats{emax, emin, phmax, phmin, pleft};
    }

    void simulateTurn(int stepID) 
    {
        CalcOp op(stepID, mAcc, mParticles, mCollHits);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size()), op);
    }
    
    void simulateTurns(int n, std::string filePath = "", int saveFreq = 1)
    {
        if (filePath.empty())
            std::cout << "Will not save particle path data" << std::endl;
        else  {
            createTurnFileHeader(filePath, 1 + n/saveFreq); // +1 as we are saving the first starting configuration as well
            writeDistribution(filePath);
        }

        std::vector<T> E_ramp;
        std::vector<std::pair<T, T>> ext_collimators;
        if (mType > NO_RAMP) {
            readRamp<T>(n, E_ramp, mType);
            readCollimators<T>(n, ext_collimators, E_ramp);
            mAcc.setE(E_ramp[0], true);
            mAcc.coll_bot = ext_collimators[0].first;
            mAcc.coll_top = ext_collimators[0].second;
        }

        std::cout << "Tracking " << size() << " particles for " << n << " turns" << std::endl;
        std::cout << "Starting simulation..." << std::endl;

        common::SilentTimer timer;

        timer.start();
        for (int i = 0; i < n; ++i) {
            if (mType > NO_RAMP) {
                mAcc.setE(E_ramp[i]);

                // Caluclated from LHC_ramp.dat
                const T k = 2.9491187074838457087e-07;
                mAcc.V_rf = (6 + k*i)*1e6;

                // Comment out to remove motor motion
                if (!ext_collimators.empty()) {
                    mAcc.coll_bot = ext_collimators[i].first;
                    mAcc.coll_top = ext_collimators[i].second;
                }
            }

            simulateTurn(i);

            // reduce the memory footprint a little if it's not necessary
            if (!filePath.empty() && (i + 1) % saveFreq == 0)
                writeDistribution(filePath);

            const int d = n/10;
            if (d > 0 && i % d == 0) {
                int percent = 10*i/d;
                ParticleStats stats = getStats();
                std::cout << "\t" << std::setw(9) << i << " of " << n << " turns (" << percent << "%).";

                int pleft = (100*stats.pleft)/size();
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

    void runLossmap(int seconds)
    {
        T freq = mAcc.f_rev;
        int turns = seconds*freq;
    
        //simulateTurns(turns); 
        simulateTurns(turns, stron::PATH_FILE, 11245); 
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
                << "(approx. after " << std::setprecision(3) << (double(tlast)/freq) << " s)\n";
        }
    
    }

    size_t size() const { return mParticles.size(); }
    std::vector<int> getCollHits() const { return mCollHits; }

private:
    struct CalcOp
    {
        CalcOp(int stepID, const Acc& acc, PColl& particles, std::vector<int>& collHits)
            : mStepID(stepID), mAcc(acc), mParticles(particles), mCollHits(collHits)
        {}

        bool particleInactive(size_t index) const
        {
            return !mCollHits.empty() && mCollHits[index] != -1;
        }

        bool particleCollided(size_t index) const
        {
            return !mCollHits.empty() // not populated means we don't look for collisions
                    && mCollHits[index] == -1 // not previously hit
                    && (mParticles.momentum[index] >= mAcc.coll_top || mParticles.momentum[index] <= mAcc.coll_bot);
        }

        bool outsideBucket(size_t index) const 
        {
            return mParticles.phase[index] < -0.1;
        }

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t n = range.begin(), N = range.end(); n < N; ++n) {
                if (particleInactive(n)) 
                    continue;
                
                using std::sin;
                using cnst::pi;

                T deltaRef = mAcc.E() - mAcc.E_prev();

                int Ns = 100;
                // Particles outside of the bucket does not need to be tracked as carefully
                if (outsideBucket(n)) Ns = 1;

                T& momentum = mParticles.momentum[n];
                T& phase = mParticles.phase[n];
                for (int i = 0; i < Ns; ++i) {
                    momentum += (mAcc.V_rf*(sin(phase)) - deltaRef)/T(Ns);
                    auto p = mAcc.calcParticleProp(momentum, 0.0);
                    phase -= T(2)*pi*mAcc.h_rf*p.eta/(T(Ns)*p.b2*mAcc.E())*momentum;
                }

                if (particleCollided(n)) mCollHits[n] = mStepID;
            }
        }
    private:
        int mStepID;
        const Acc& mAcc;
        PColl& mParticles;
        std::vector<int>& mCollHits;
    };

    int mStepID;
    Acc mAcc;
    // Particle properties
    PColl mParticles;
    std::vector<int> mCollHits; // Turn hitting collimator

    RAMP_TYPE mType;
};

} // namespace stron
