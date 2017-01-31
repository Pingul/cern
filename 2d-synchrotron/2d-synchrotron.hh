#include <iostream>
#include <vector>
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

namespace CONST {

constexpr double pi = 3.14159265359;
constexpr double c = 299792458.0; // m/s
constexpr double m_proton = 938.2796e6; // eV

}; // namespace CONST

namespace jwc {

// Same as used in the python visualisation code
static constexpr double FRAME_X_LOW = -2*CONST::pi;
static constexpr double FRAME_X_HIGH = 4*CONST::pi;
static constexpr double FRAME_Y_LOW = -2e9;
static constexpr double FRAME_Y_HIGH = 2e9;

static constexpr const char* PATH_FILE = "calc/particles.dat";
static constexpr const char* LINE_FILE = "calc/lines.dat";
static constexpr const char* COLL_FILE = "calc/coll.dat";
static constexpr const char* STARTDIST_FILE = "calc/startdist.dat";
static constexpr const char* ENDDIST_FILE = "calc/enddist.dat";
static constexpr const char* SIXTRACK_TEST_FILE = "calc/toymodel_track.dat";
// static constexpr const char* RAMP_FILE = "resources/ramp.txt";
static constexpr const char* RAMP_FILE = "resources/LHC_ramp.dat";
static constexpr const char* COLL_MOTOR_FILE = "resources/motor_tcp.txt";

template <typename T>
struct Accelerator
{
    using ValType = T;
    using Acc = Accelerator<T>;

private:
    // Use accessor methods for these instead
    T mE_ref;
    T mE_pref;
public:
    T rf_voltage;
    T rf_freq;
    T harmonic;
    T m_compaction;
    T coll_top;
    T coll_bot;
    T revolution_freq;
    T w_revolution_freq;

    static Acc getLHC() 
    {
        // Parameters for LHC can be found at http://atlas.physics.arizona.edu/~kjohns/downloads/linac/UTclass2lhc.pdf
        Acc acc;
        acc.mE_ref = T(450e9); // eV
        acc.mE_pref = acc.mE_ref;
        acc.rf_voltage = T(6e6); // V
        acc.rf_freq = T(398765412.66);
        acc.harmonic = T(35640);
        acc.m_compaction = T(0.0003225);
        acc.coll_top = T(0.5e9); // ∆eV
        acc.coll_bot = T(-0.5e9);
        acc.revolution_freq = T(acc.rf_freq/acc.harmonic); // Hz
        acc.w_revolution_freq = 2*CONST::pi*acc.revolution_freq; // Hz
        return acc;
    }

    static Acc getLHC_NOCOLL()
    {
        Acc a = getLHC();
        a.coll_top = std::numeric_limits<T>::max();
        a.coll_bot = std::numeric_limits<T>::min();;
        return a;
    }

    void setE(T v, bool reset = false) { mE_pref = mE_ref; if (reset) mE_pref = v; mE_ref = v; }
    T E() const { return mE_ref; }
    T E_prev() const { return mE_pref; }
};

template <typename T>
inline T hamiltonian(const Accelerator<T>& acc, T deltaE, T phase)
{
    const T rev = acc.w_revolution_freq;
    const T gamma = (acc.E() + deltaE)/CONST::m_proton;
    const T gamma_2 = T(1)/(gamma*gamma);
    const T eta = gamma_2 - acc.m_compaction;
    const T beta2 = T(1) - gamma_2;
    const T beta = std::sqrt(beta2);
    const T k = acc.harmonic*rev/(beta*CONST::c);
    const T Omega2 = rev*rev*acc.harmonic*eta*acc.rf_voltage/(T(2)*CONST::pi*beta*acc.E());

    const T H = T(1)/2*beta2*pow(CONST::c*k*eta*deltaE/acc.E(), 2) - Omega2*std::cos(phase);
    return H;
}

template <typename T>
struct ToyModel 
{
    using Accelerator = Accelerator<T>;

    enum RAMP_TYPE
    {
        NO_RAMP,
        LHC_RAMP,
        SEMI_AGGRESSIVE_RAMP,
        AGGRESSIVE_RAMP,
    };

    ToyModel(size_t n, RAMP_TYPE type)
        : mAcc(Accelerator::getLHC()), mType(type)
    {
        mEnergy.reserve(n);
        mPhase.reserve(n);
        mCollHits.assign(n, -1);

        std::random_device rdev;
        std::mt19937 generator(rdev());

        // Good distribution for the full frame
        // std::normal_distribution<> e_dist(0, 0.2e9);
        // std::normal_distribution<> ph_dist(CONST::pi, CONST::pi);

        std::normal_distribution<> e_dist(0, 0.15e9);
        std::normal_distribution<> ph_dist(CONST::pi, CONST::pi/5);
        for (size_t i = 0; i < n; ++i) {
            const T deltaE = e_dist(generator);
            const T phase = ph_dist(generator);
            mEnergy.push_back(deltaE); 
            mPhase.push_back(phase); 
        }

        writeSingleDistribution(STARTDIST_FILE);
    }

    // For parameter passing in the next constructors
    struct LineSim {};
    struct LossAnalysis {};
    struct SixTrackTest {};

    ToyModel(RAMP_TYPE type, LineSim)
        : mAcc(Accelerator::getLHC_NOCOLL()), mType(type)
    {
        for (T y = FRAME_Y_LOW; y < FRAME_Y_HIGH; y += FRAME_Y_HIGH/9.0) {
            mEnergy.push_back(y);
            mPhase.push_back(FRAME_X_LOW);
        }

        switch (mType) {
            default:
            case NO_RAMP:
            case LHC_RAMP:
                for (T x = FRAME_X_LOW; x < FRAME_X_HIGH; x += CONST::pi/4.0) {
                    mEnergy.push_back(0.0f);
                    mPhase.push_back(x);
                }
                for (T y = FRAME_Y_LOW; y < FRAME_Y_HIGH; y += FRAME_Y_HIGH/9.0) {
                    mEnergy.push_back(y);
                    mPhase.push_back(FRAME_X_HIGH);
                }
                break;
            case AGGRESSIVE_RAMP:
                break;
            case SEMI_AGGRESSIVE_RAMP:
                for (T x = -CONST::pi; x < FRAME_X_HIGH; x += 2*CONST::pi) {
                    T d = CONST::pi/2;
                    for (T delta = -d; delta < d; delta += CONST::pi/5) {
                        mEnergy.push_back(0.0f);
                        mPhase.push_back(x + delta);
                    }
                }
                break;
        }
    }

    ToyModel(RAMP_TYPE type, LossAnalysis) 
        : mAcc(Accelerator::getLHC()), mType(type)
    {
        const T Emax = 7e8;
        const T Emin = -Emax;
        const T Estep = (Emax - Emin)/2000;
        const T PHmax = 2.01*CONST::pi;
        const T PHmin = 0;
        const T PHstep = (PHmax - PHmin)/2000;

        for (T de = Emin; de < Emax; de += Estep) {
            for (T ph = PHmin; ph < PHmax; ph += PHstep) {
                const T H = hamiltonian(mAcc, de, ph);
                if (H > 1.17e5 && H < 1.199e5) {
                    mEnergy.push_back(de);
                    mPhase.push_back(ph);
                }
            }
        }
        mCollHits.assign(size(), -1);

        writeSingleDistribution(STARTDIST_FILE);
    }

    ToyModel(int n, RAMP_TYPE type, LossAnalysis)
        : mAcc(Accelerator::getLHC()), mType(type)
    {    
        mEnergy.reserve(n);
        mPhase.reserve(n);
        mCollHits.assign(n, -1);

        std::random_device rdev;
        std::mt19937 generator(rdev());

        std::uniform_real_distribution<> e_dist(-0.5e9, 0.5e9);
        std::uniform_real_distribution<> ph_dist(0, 2*CONST::pi);
        int count = 0;
        while (count < n) {
            const T deltaE = e_dist(generator);
            const T phase = ph_dist(generator);
            const T H = hamiltonian(mAcc, deltaE, phase);
            if (H > 1.193e5 && H < 1.195e5) {
                mEnergy.push_back(deltaE); 
                mPhase.push_back(phase); 
                count++;
            }
        }
        writeSingleDistribution(STARTDIST_FILE);
    }

    ToyModel(SixTrackTest)
        : mAcc(Accelerator::getLHC()), mType(LHC_RAMP)
    {
        mEnergy.push_back(T(0.41646726612503554e6));
        mPhase.push_back(CONST::pi);
        // mPhase.push_back(0);
    }

    size_t size() const { return mEnergy.size(); }

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
        //     ∆energy / phase 
        for (size_t i = 0; i < size(); ++i) {
            std::stringstream ss;
            ss << std::setprecision(16) << mEnergy[i] << "," << mPhase[i] << std::endl;
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
            file << i << ", " << mCollHits[i] << ", " << std::setprecision(16) << mPhase[i] << ", " << mEnergy[i] << std::endl;
        }
    }

    void takeTimestep(int stepID) 
    {
        CalcOp op(stepID, mAcc, mEnergy, mPhase, mCollHits);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size()), op);
    }

    void loadRamp(int steps, std::vector<T>& E_ramp, std::vector<std::pair<T, T>>& collimators) 
    {
        E_ramp.reserve(steps);
        collimators.reserve(steps);

        T data;
        switch (mType) {
            default:
            case LHC_RAMP: {
                std::ifstream ramp_file(RAMP_FILE);
                std::cout << "Reading '" << RAMP_FILE << "'..." << std::endl;
                for (int i = 0; i < steps; ++i) {
                    ramp_file >> skip >> data;
                    E_ramp.push_back(data*1e6);
                }
                ramp_file.close();

                std::ifstream coll_file(COLL_MOTOR_FILE);
                std::cout << "Reading '" << COLL_MOTOR_FILE << "'..." << std::endl;
                for (int i = 0; i < steps; ++i) {
                    T bot_coll, top_coll;
                    coll_file >> skip >> bot_coll >> top_coll;
                    collimators.push_back(std::make_pair(bot_coll*E_ramp[i], top_coll*E_ramp[i]));
                }
                break;
            }
            case AGGRESSIVE_RAMP:
                for (int i = 0; i < steps; ++i) E_ramp.push_back(450e9 + i*1e7);
                break;
            case SEMI_AGGRESSIVE_RAMP:
                for (int i = 0; i < steps; ++i) E_ramp.push_back(450e9 + i*3e6);
                break;
        }

    }

    void takeTimesteps(int n, std::string filePath = "", int saveFreq = 1)
    {
        if (filePath.empty())
            std::cout << "Will not save particle path data" << std::endl;
        else  {
            createTurnFileHeader(filePath, n + 1); // +1 as we are saving the first starting configuration as well
            writeDistribution(filePath);
        }

        std::vector<T> E_ramp;
        std::vector<std::pair<T, T>> ext_collimators;
        if (mType > NO_RAMP) {
            loadRamp(n, E_ramp, ext_collimators);
            mAcc.setE(E_ramp[0], true);
        }

        std::cout << "Tracking " << size() << " particles for " << n << " turns" << std::endl;
        std::cout << "Starting simulation..." << std::endl;

        SilentTimer timer;

        timer.start();
        for (int i = 0; i < n; ++i) {
            if (mType > NO_RAMP) {
                // T deltaE = E_ramp[i] - mAcc.E();
                mAcc.setE(E_ramp[i]);

                // // Adjust all particle energies
                // for (T& e : mEnergy) 
                //     e -= deltaE;

                // Caluclated from LHC_ramp.dat
                const T k = 2.9491187074838457087e-07;
                mAcc.rf_voltage = (6 + k*n)*1e6;

                if (!ext_collimators.empty()) {
                    mAcc.coll_bot = ext_collimators[i].first;
                    mAcc.coll_top = ext_collimators[i].second;
                }
            }

            takeTimestep(i);

            // reduce the memory footprint a little if it's not necessary
            if (!filePath.empty() && (i + 1) % saveFreq == 0) writeDistribution(filePath);

            const int d = n/10;
            if (i % d == 0) {
                int percent = 10*i/d;
                std::cout << "\t" << std::setw(9) << i << " of " << n << " turns (" << percent << "%)" << std::endl;
            }
        }

        double ms = timer.stop();

        auto dur = MillisecondsToElapsedTime(unsigned(ms));
        std::cout << "Finished in ";
        if (dur.d > 0) std::cout << dur.d << "d ";
        if (dur.h > 0) std::cout << dur.h << "h ";
        if (dur.m > 0) std::cout << dur.m << "m ";
        if (dur.s > 0) std::cout << dur.s << "s ";
        std::cout << dur.ms << "ms" << std::endl;

        writeSingleDistribution(ENDDIST_FILE);
    }

    const std::vector<T>& getEnergy() const { return mEnergy; }
    const std::vector<T>& getPhase() const { return mPhase; }
    std::vector<int> getCollHits() const { return mCollHits; }

private:
    struct CalcOp
    {
        CalcOp(int stepID, const Accelerator& acc, std::vector<T>& energy, 
               std::vector<T>& phase, std::vector<int>& collHits)
            : mStepID(stepID), mAcc(acc), mEnergy(energy), mPhase(phase), mCollHits(collHits)
        { }

        bool particleInactive(size_t index) const
        {
            return !mCollHits.empty() && mCollHits[index] != -1;
        }

        bool particleCollided(size_t index) const
        {
            return !mCollHits.empty() // not populated means we don't look for collisions
                    && mCollHits[index] == -1 // not previously hit
                    && (mEnergy[index] >= mAcc.coll_top || mEnergy[index] <= mAcc.coll_bot);
        }

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            const T e = 1.0;
            // const T m = 938e6;
            const T B2_s = 1 - std::pow(CONST::m_proton/mAcc.E(), 2);
            for (size_t n = range.begin(), N = range.end(); n < N; ++n) {
                if (particleInactive(n)) 
                    continue;

                // We changed the reference energy before, update the ∆E accordingly
                T deltaRef = mAcc.E() - mAcc.E_prev();
                mEnergy[n] -= deltaRef;

                // Give the particle a kick
                // mEnergy[n] -= mAcc.rf_voltage*std::sin(deltaRef/mAcc.rf_voltage);

                mEnergy[n] += e*mAcc.rf_voltage*std::sin(mPhase[n]);
                T B2 = T(1) - std::pow(CONST::m_proton/(mAcc.E() + mEnergy[n]), 2);
                T eta = T(1) - B2 - mAcc.m_compaction;
                mPhase[n] -= T(2)*CONST::pi*mAcc.harmonic*eta/(B2_s*mAcc.E())*mEnergy[n];
                if (particleCollided(n)) mCollHits[n] = mStepID;
            }
        }
    private:
        int mStepID;

        const Accelerator& mAcc;

        std::vector<T>& mEnergy;
        std::vector<T>& mPhase;

        std::vector<int>& mCollHits;
    };

    int mStepID;

    Accelerator mAcc;
    // Particle properties
    std::vector<T> mEnergy; // eV
    std::vector<T> mPhase; // radians

    std::vector<int> mCollHits;

    RAMP_TYPE mType;
};

inline double synchrotron_frequency()
{
    const auto acc = ToyModel<double>::Accelerator::getLHC();
    const auto omega = std::sqrt(std::abs(hamiltonian<double>(acc, 0, 0)));
    const auto freq_turns = omega/acc.w_revolution_freq;
    return freq_turns;
}


}; // namespace jwc
