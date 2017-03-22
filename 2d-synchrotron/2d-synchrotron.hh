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

#ifdef IO_TO_SAME_DIR
#define RESOURCE_DIR "."
#define OUTPUT_DIR "."
#else
#define RESOURCE_DIR "resources"
#define OUTPUT_DIR "calc"
#endif

namespace twodsynch {

namespace cnst {
constexpr double pi = 3.14159265359;
constexpr double c = 299792458.0; // m/s
constexpr double m_proton = 938.2796e6; // eV
} // namespace cnst


// Same as used in the python visualisation code
static constexpr double FRAME_X_LOW = -2*cnst::pi;
static constexpr double FRAME_X_HIGH = 4*cnst::pi;
static constexpr double FRAME_Y_LOW = -2e9;
static constexpr double FRAME_Y_HIGH = 2e9;

static constexpr const char* PATH_FILE = OUTPUT_DIR"/particles.dat";
static constexpr const char* LINE_FILE = OUTPUT_DIR"/lines.dat";
static constexpr const char* COLL_FILE = OUTPUT_DIR"/coll.dat";
static constexpr const char* STARTDIST_FILE = OUTPUT_DIR"/startdist.dat";
static constexpr const char* ENDDIST_FILE = OUTPUT_DIR"/enddist.dat";
static constexpr const char* SIXTRACK_TEST_FILE = OUTPUT_DIR"/toymodel_track.dat";
// static constexpr const char* RAMP_FILE = "resources/ramp.txt";
static constexpr const char* RAMP_FILE = RESOURCE_DIR"/LHC_ramp.dat";
static constexpr const char* COLL_MOTOR_FILE = RESOURCE_DIR"/motor_tcp.txt";

enum RAMP_TYPE
{
    NO_RAMP,
    LHC_RAMP,
    EXTENDED_LHC_RAMP,
    SEMI_AGGRESSIVE_RAMP,
    AGGRESSIVE_RAMP,
};

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
    T C;
    T V_rf;
    T f_rf;
    T h_rf;
    T k_rf;
    T m_compaction;
    T coll_top;
    T coll_bot;
    T f_rev;
    T w_rev;

    static Acc getLHC() 
    {
        // Parameters for LHC can be found at http://atlas.physics.arizona.edu/~kjohns/downloads/linac/UTclass2lhc.pdf
        Acc acc;
        acc.C = 26658.8832;
        acc.mE_ref = T(450e9); // eV
        acc.mE_pref = acc.mE_ref;
        acc.V_rf = T(6e6); // V
        //acc.f_rf = T(398765412.66);
        //acc.f_rf = T(400.8e6);
        acc.h_rf = T(35640);
        acc.k_rf = acc.h_rf*T(2)*cnst::pi/acc.C;
        acc.m_compaction = T(0.0003225);
        acc.coll_top = T(0.5e9); // ∆eV
        acc.coll_bot = T(-0.5e9);
        
        auto p = acc.calcParticleProp(0.0, 0.0); // This is ok: we only want b
        acc.f_rev = p.b*cnst::c/acc.C;
        acc.f_rf = acc.f_rev*acc.h_rf;

        //acc.f_rev = T(acc.f_rf/acc.h_rf); // Hz
        acc.w_rev = 2*cnst::pi*acc.f_rev; // Hz
        return acc;
    }

    static Acc getLHC_NOCOLL()
    {
        Acc a = getLHC();
        a.coll_top = std::numeric_limits<T>::max();
        a.coll_bot = std::numeric_limits<T>::min();;
        return a;
    }

    struct ParticleProp 
    { 
        T pc;
        T g;    // gamma
        T g_2;    // 1/gamma^2
        T eta;
        T b2;    // beta^2
        T b;    // beta
        T W2;    // Omega2
    };
    ParticleProp calcParticleProp(T de, T phase) const 
    {
        ParticleProp p;
        //p.g = (E() + de)/cnst::m_proton;

        // We treat ∆E really as ∆pc
        const T pc_E0 = (E() + de)/cnst::m_proton;
        p.g = std::sqrt(T(1) + pc_E0*pc_E0);
        p.g_2 = T(1)/(p.g*p.g);
        p.eta = p.g_2 - m_compaction;
        p.b2 = T(1) - p.g_2;
        p.b = std::sqrt(p.b2);
        p.W2 = w_rev*w_rev*h_rf*p.eta*V_rf/(T(2)*cnst::pi*p.b2*E())*std::cos(lag_phase()); 
        return p;
    }


    void setE(T v, bool reset = false) { mE_pref = mE_ref; if (reset) mE_pref = v; mE_ref = v; }
    T E() const { return mE_ref; }
    T E_prev() const { return mE_pref; }
    T lag_phase() const { return cnst::pi - (std::asin((E() - E_prev())/V_rf)); }
    //T lag_phase() const { return cnst::pi - 0.3; }
};

template <typename T>
inline T hamiltonian(const Accelerator<T>& acc, T de, T ph)
{
    using std::cos;
    using std::sin;
    using cnst::pi;

    auto p = acc.calcParticleProp(de, ph);
    const T ph_s = acc.lag_phase();
    const T A = -T(0.5)*acc.h_rf*p.eta/(p.b2*acc.E());
    const T H = A*de*de + acc.V_rf/(T(2)*pi)*(cos(ph) - cos(ph_s) + (ph - ph_s)*sin(ph_s));
    return H;
}

template <typename T>
inline T levelCurve(const Accelerator<T>& acc, T ph, T H, T sign=1.0) 
{
    using std::cos;
    using std::sin;
    using std::sqrt;
    using cnst::pi;

    sign = sign/std::abs(sign);

    const T ph_s = acc.lag_phase();
    const T threshold = 0.005;
    T de = 0.0;
    do {
        const auto p = acc.calcParticleProp(de, ph);
        const T A = -T(0.5)*acc.h_rf*p.eta/(p.b2*acc.E());
        de = sign*sqrt(T(1.0)/A*(H - acc.V_rf/(T(2)*pi)*(cos(ph) - cos(ph_s) + (ph - ph_s)*sin(ph_s))));
    } while (std::abs(hamiltonian(acc, de, ph) - H) > threshold);
    return de;
}

template <typename T>
inline T separatrix(const Accelerator<T>& acc)
{
    return hamiltonian<T>(acc, 0.0, cnst::pi - acc.lag_phase());
}

template <typename T>
inline void writePhasespaceFrame(const Accelerator<T>& acc, std::string filePath)
{
    std::cout << "Writing frame to '" << filePath << "'" << std::endl;
    std::ofstream file(filePath.c_str());
    if (!file.is_open())
        throw std::runtime_error("could not open file");

    const T ph_step = 0.005;
    const int ph_steps = std::ceil((FRAME_X_HIGH - FRAME_X_LOW)/ph_step);

    const T Hstart = hamiltonian<T>(acc, 0.0, cnst::pi) + 5e4;
    const T Hstep = 2.0e5;
    const T Hdstep = 1.0e5;
    const int Hsteps = 20;
    //const int Hsteps = 0;

    int lines = 2 + 2*Hsteps;

    const T Hseparatrix = separatrix<T>(acc);

    file << lines << "," << ph_steps << std::endl;
    for (T ph = FRAME_X_LOW; ph <= FRAME_X_HIGH; ph += ph_step) {
        T de = levelCurve(acc, ph, Hseparatrix);
        file << std::setprecision(16) << de  << "," << ph << "," << Hseparatrix << std::endl
                                      << -de  << "," << ph << "," << Hseparatrix << std::endl;
        
        int n = 0;
        T H = Hstart;
        while (n++ < Hsteps) {
            de = levelCurve(acc, ph, H);
            file << std::setprecision(16) << de  << "," << ph << "," << H << std::endl
                                          << -de  << "," << ph << "," << H << std::endl;
            H += Hstep + T(n)*Hdstep;
        }
    }
}

template <typename T>
inline void readRamp(int steps, std::vector<T>& E_ramp, RAMP_TYPE type)
{
    using common::skip;

    E_ramp.reserve(steps);

    T data;
    switch (type) {
        default:
        case LHC_RAMP: {
            std::ifstream ramp_file(RAMP_FILE);
            std::cout << "Reading '" << RAMP_FILE << "'..." << std::endl;
            for (int i = 0; i < steps; ++i) {
                ramp_file >> skip >> data;
                E_ramp.push_back(data*1e6);
            }
            ramp_file.close();

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
struct ToyModel 
{
    using Acc = Accelerator<T>;

    ToyModel(size_t n, RAMP_TYPE type)
        : mAcc(Acc::getLHC()), mType(type)
    {
        initStorage(n);

        std::random_device rdev;
        std::mt19937 generator(rdev());

        std::normal_distribution<> e_dist(0, 0.15e9);
        std::normal_distribution<> ph_dist(cnst::pi, cnst::pi/5);
        for (size_t i = 0; i < n; ++i) {
            const T deltaE = e_dist(generator);
            const T phase = ph_dist(generator);
            mEnergy.push_back(deltaE); 
            mPhase.push_back(phase); 
        }
        writeSingleDistribution(STARTDIST_FILE);
    }


    // For parameter passing in the next constructors
    // These are all different distributions to choose from
    struct AroundSeparatrix {};
    struct ActionValues {};
    struct LinearlyDecaying {};
    struct ExponentiallyDecaying {};
    
    struct SixTrackTest {};

    ToyModel(int n, RAMP_TYPE type, AroundSeparatrix)
        : mAcc(Acc::getLHC()), mType(type)
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
            if ((sep - 5000) < H && H < (sep + 2000)) {
                mEnergy.push_back(deltaE); 
                mPhase.push_back(phase); 
                count++;
            }
        }
        writeSingleDistribution(STARTDIST_FILE);
    }

    ToyModel(int N, RAMP_TYPE type, ActionValues)
        : mAcc(Acc::getLHC()), mType(type)
    {
        const T sep = separatrix(mAcc);
        //std::vector<T> d_actions = {-1e4, -9e3, -8e3, -7e3, -6e3, -5e3, -4e3, -3e3, -2e3, -1e3, -100};
        std::vector<T> d_actions;
        for (int i = 1; i <= 80; ++i) {
            d_actions.push_back(T(-8000 + 100*i));
        }
        const int n = N/(2*d_actions.size());
        initStorage(2*n*d_actions.size());

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
                    mEnergy.push_back(energy);
                    mPhase.push_back(phase);
                }
            }
        }
        std::cout << "Tried to create " << N << " particles, only initialized " << size() << std::endl;
        writeSingleDistribution(STARTDIST_FILE);
    }

    ToyModel(int n, RAMP_TYPE type, LinearlyDecaying)
        :    mAcc(Acc::getLHC()), mType(type)
    {
        const T sep = separatrix(mAcc);

        std::random_device rdev;
        std::mt19937 generator(rdev());

        // Distribution as (relative to separatrix)
        //  -15k        -10k          5k
        //    | constant  | linear dec.|
        std::vector<T> Hs{sep - 15e3, sep - 10e3, sep + 5e3};
        std::vector<T> prob{1.0, 1.0, 0.0};
        std::piecewise_linear_distribution<> H_dist(Hs.begin(), Hs.end(), prob.begin());

        std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);

        for (int i = 0; i < n; ++i) {
            const T phase = uni_dist(generator);
            const T sign = uni_dist(generator) < cnst::pi ? 1.0 : -1.0;
            const T action = H_dist(generator);
            const T energy = levelCurve(mAcc, phase, action, sign);
            if (std::isnan(energy)) { --i; continue; }
            mEnergy.push_back(energy);
            mPhase.push_back(phase);
        }
        std::cout << "Initialized " << size() << " particles" << std::endl;
        writeSingleDistribution(STARTDIST_FILE);
    }

    
    ToyModel(const std::string filePath, RAMP_TYPE type)
        : mAcc(Acc::getLHC()), mType(type)
    {
        readDistribution(filePath);
        std::cout << "Initialized " << size() << " particles" << std::endl;
        writeSingleDistribution(STARTDIST_FILE);
    }

    ToyModel(SixTrackTest, T energy)
        : mAcc(Acc::getLHC()), mType(LHC_RAMP)
    {
        mEnergy.push_back(energy);
        mPhase.push_back(cnst::pi);
    }

    void initStorage(int n)
    {
        mEnergy.reserve(n);
        mPhase.reserve(n);
        mCollHits.assign(n, -1);
        mLost.assign(n, -1);
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
            mEnergy.push_back(de);
            mPhase.push_back(phase);
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
            ss << std::setprecision(16) << mEnergy[i] << "," << mPhase[i] << "," << hamiltonian(mAcc, mEnergy[i], mPhase[i]) << std::endl;
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

    void writeLostTurns(std::string filePath) const
    {
        std::cout << "Saving lost turns to '" << filePath << "'" << std::endl;
        std::ofstream file(filePath.c_str());
        if (!file.is_open())
            throw std::runtime_error("could not open file");
        else if (mLost.empty())
            throw std::runtime_error("no data to write");

        for (size_t i = 0; i < size(); ++i) {
            if (mLost[i] == -1) continue;
            std::stringstream ss;
            int delta = mCollHits[i] - mLost[i];
            file << i << "," << mLost[i] << "," << mCollHits[i] << "," << delta << std::endl;
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
                emax = mEnergy[i] > emax ? mEnergy[i] : emax;
                emin = mEnergy[i] < emin ? mEnergy[i] : emin;
                phmax = mPhase[i] > phmax ? mPhase[i] : phmax;
                phmin = mPhase[i] < phmin ? mPhase[i] : phmin;
            }
        }
        return ParticleStats{emax, emin, phmax, phmin, pleft};
    }

    void simulateTurn(int stepID) 
    {
        CalcOp op(stepID, mAcc, mEnergy, mPhase, mCollHits, mLost);
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
        simulateTurns(turns, twodsynch::PATH_FILE, 11245); 
        writeCollHits(COLL_FILE);
    
        int ilast = -1;
        for (size_t i = 0; i < mCollHits.size(); ++i) {
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

    size_t size() const { return mEnergy.size(); }
    const std::vector<T>& getEnergy() const { return mEnergy; }
    const std::vector<T>& getPhase() const { return mPhase; }
    std::vector<int> getCollHits() const { return mCollHits; }

private:
    struct CalcOp
    {
        CalcOp(int stepID, const Acc& acc, std::vector<T>& energy, 
               std::vector<T>& phase, std::vector<int>& collHits, std::vector<int>& lost)
            : mStepID(stepID), mAcc(acc), mEnergy(energy), mPhase(phase), mCollHits(collHits), mLost(lost)
        {}

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

        bool particleLost(size_t index) const 
        {
            // The constant of motion/hamiltionian should be used here, but we've had some small
            // discrepancies in the past, and so it's probably not too reliable. Instead, we 
            // define a lost particle to be one that has some arbitrary (low) phase, as we then
            // know that it is outside of the bucket.

            return !mLost.empty() && 
                    mLost[index] == -1 && 
                    mPhase[index] < -0.1;
            // Below does not work
            //const T mSeparatrix = hamiltonian(mAcc, 0.0, cnst::pi - mAcc.lag_phase());
            //return hamiltonian(mAcc, mPhase[index], mEnergy[index]) > mSeparatrix;
        }

        bool isPreviouslyLost(size_t index) const
        {
            return !mLost.empty() && mLost[index] > -1;
        }

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            //auto prop_s = mAcc.calcParticleProp(0.0, 0.0);
            //const T B2_s = prop_s.b2;
            for (size_t n = range.begin(), N = range.end(); n < N; ++n) {
                if (particleInactive(n)) 
                    continue;
                
                using std::sin;
                using cnst::pi;

                // We changed the reference energy before, update the ∆E accordingly
                // This constitutes the sin(lag_phase) term in actual stepping function
                T deltaRef = mAcc.E() - mAcc.E_prev();
                //mEnergy[n] -= deltaRef;

                int Ns = 100;
                // Particles outside of the bucket does not need to be tracked as carefully
                if (isPreviouslyLost(n)) Ns = 1;

                for (int i = 0; i < Ns; ++i) {
                    mEnergy[n] += (mAcc.V_rf*(sin(mPhase[n])) - deltaRef)/T(Ns);
                    auto p = mAcc.calcParticleProp(mEnergy[n], 0.0);
                    mPhase[n] -= T(2)*pi*mAcc.h_rf*p.eta/(T(Ns)*p.b2*mAcc.E())*mEnergy[n];
                }

                if (particleCollided(n)) mCollHits[n] = mStepID;
                if (particleLost(n)) mLost[n] = mStepID;
            }
        }
    private:
        int mStepID;

        const Acc& mAcc;

        std::vector<T>& mEnergy;
        std::vector<T>& mPhase;

        std::vector<int>& mCollHits;
        std::vector<int>& mLost;

    };

    int mStepID;

    Acc mAcc;
    // Particle properties
    std::vector<T> mEnergy; // eV -- note that this is '∆pc'
    std::vector<T> mPhase; // radians

    std::vector<int> mCollHits; // Turn hitting collimator
    std::vector<int> mLost;        // Turn considered 'lost' from the bunch

    RAMP_TYPE mType;
};

void generatePhasespaceLines(int seconds)
{
    // will generate 1 per second
    std::cout << "Generate phasespace lines" << std::endl;
    auto acc = Accelerator<double>::getLHC();
    int freq = int(acc.f_rev);
    int turns = seconds*freq;

    std::vector<double> E;
    readRamp(turns, E, LHC_RAMP);

    for (int i = 0; i < seconds; ++i) {
        int turn = i*freq;
        acc.setE(E[turn]);
        acc.setE(E[turn + 1]);
        
        // Caluclated from LHC_ramp.dat
        const double k = 2.9491187074838457087e-07;
        acc.V_rf = (6 + k*turn)*1e6;
        
        std::stringstream ss;
        ss << "phasespace/" << i << "lines.dat";
        writePhasespaceFrame(acc, ss.str());
    }
}

} // namespace twodsynch
