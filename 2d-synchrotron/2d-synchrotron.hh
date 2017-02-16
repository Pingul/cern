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

namespace twodsynch {

namespace cnst {
constexpr double pi = 3.14159265359;
constexpr double c = 299792458.0; // m/s
constexpr double m_proton = 938.2796e6; // eV
}; // namespace cnst


// Same as used in the python visualisation code
static constexpr double FRAME_X_LOW = -2*cnst::pi;
static constexpr double FRAME_X_HIGH = 4*cnst::pi;
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

enum RAMP_TYPE
{
    NO_RAMP,
    LHC_RAMP,
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
        acc.f_rf = T(398765412.66);
        acc.h_rf = T(35640);
		acc.k_rf = acc.h_rf*T(2)*cnst::pi/acc.C;
        acc.m_compaction = T(0.0003225);
        acc.coll_top = T(0.5e9); // ∆eV
        acc.coll_bot = T(-0.5e9);
        acc.f_rev = T(acc.f_rf/acc.h_rf); // Hz
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
		T g;	// gamma
		T g_2;	// 1/gamma^2
		T eta;
		T b2;	// beta^2
		T b;	// beta
		T W2;	// Omega2
	};
	ParticleProp calcParticleProp(T de, T phase) const 
	{
		ParticleProp p;
		p.g = (E() + de)/cnst::m_proton;
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
	using cnst::c;
	using cnst::pi;

	auto p = acc.calcParticleProp(0.0, ph); // calculating properties for the synchronous particle
	const T ph_s = acc.lag_phase(); // I do not know why I need to add the pi here to get things to work
	const T phi_dot = -p.b*c*acc.k_rf*p.eta*de/acc.E();
	const T H = T(0.5)*phi_dot*phi_dot - p.W2/cos(ph_s)*(cos(ph) + ph*sin(ph_s));
    return H;
}

template <typename T>
inline T levelCurve(const Accelerator<T>& acc, T ph, T H) 
{
	using std::cos; 
	using std::sin;
	using std::sqrt;
	using cnst::c;
	using cnst::pi;

	auto p = acc.calcParticleProp(0.0, ph);
	const T ph_s = acc.lag_phase();
	const T phi_dot = sqrt(2*H + 2*p.W2/cos(ph_s)*(cos(ph) + ph*sin(ph_s)));
	const T de = phi_dot*acc.E()/(-p.b*c*acc.k_rf*p.eta);
	return de;
}

template <typename T>
inline void phasespaceLevelCurve(const Accelerator<T>& acc, std::string filePath)
{
	std::cout << "Writing level curves" << std::endl;
	std::ofstream file(filePath.c_str());
	if (!file.is_open())
		throw std::runtime_error("could not open file");

	const T ph_step = 0.01;
	const int ph_steps = int((FRAME_X_HIGH - FRAME_X_LOW)/ph_step);

	const T Hstart = hamiltonian(acc, 0.0, cnst::pi) + 1e4;
	const T Hstep = 0.6e5;
	const T Hdstep = 0.2e5;
	const int Hsteps = 20;

	int lines = 2 + 2*Hsteps;

	const T Hseparatrix = hamiltonian(acc, 0.0, cnst::pi - acc.lag_phase());

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
    using Accelerator = Accelerator<T>;

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
        // std::normal_distribution<> ph_dist(cnst::pi, cnst::pi);

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
    struct LossAnalysis {};
    struct LossAnalysis2 {};
    struct SixTrackTest {};

	ToyModel(int n, RAMP_TYPE type, LossAnalysis)
		: mAcc(Accelerator::getLHC()), mType(type)
	{
        mEnergy.reserve(n);
        mPhase.reserve(n);
        mCollHits.assign(n, -1);

        std::random_device rdev;
        std::mt19937 generator(rdev());

        std::uniform_real_distribution<> e_dist(-0.5e9, 0.5e9);
        std::uniform_real_distribution<> ph_dist(0, 2*cnst::pi);
        int count = 0;
        while (count < n) {
            const T deltaE = e_dist(generator);
            const T phase = ph_dist(generator);
            const T H = hamiltonian(mAcc, deltaE, phase);
			//if (-1000 < H && H < 1000) {
			if (H < 0) {
                mEnergy.push_back(deltaE); 
                mPhase.push_back(phase); 
                count++;
			}
        }
        writeSingleDistribution(STARTDIST_FILE);
	}

    ToyModel(int n, RAMP_TYPE type, LossAnalysis2)
        : mAcc(Accelerator::getLHC()), mType(type)
    {    
        mEnergy.reserve(n);
        mPhase.reserve(n);
        mCollHits.assign(n, -1);

        std::random_device rdev;
        std::mt19937 generator(rdev());

        std::normal_distribution<> e_dist(11720647.5520609, 1e4);
        std::normal_distribution<> ph_dist(6.150280636842277, 0.01);
        int count = 0;
        while (count < n) {
            const T deltaE = e_dist(generator);
            const T phase = ph_dist(generator);
            mEnergy.push_back(deltaE); 
            mPhase.push_back(phase); 
            count++;
        }
        writeSingleDistribution(STARTDIST_FILE);
    }


    ToyModel(SixTrackTest)
        : mAcc(Accelerator::getLHC()), mType(LHC_RAMP)
    {
        mEnergy.push_back(T(0.41646726612503554e6));
        mPhase.push_back(cnst::pi);
        // mPhase.push_back(0);
    }


	//
	// FILE IO
	//
    ToyModel(const std::string filePath, RAMP_TYPE type)
        : mAcc(Accelerator::getLHC()), mType(type)
    {
        std::cout << "Reading distribution from '" << filePath << "'" << std::endl;
        std::ifstream file(filePath.c_str());
        if (!file.is_open())
            throw std::runtime_error("could not open file");

		using common::skip;

        int n;
        file >> n >> skip;
        mEnergy.reserve(n);
        mPhase.reserve(n);
        mCollHits.assign(n, -1);
    
        for (int i = 0; i < n; ++i) {
            T de, phase;
            file >> de >> skip<char> >> phase >> skip; // we need to consume end of line
            mEnergy.push_back(de);
            mPhase.push_back(phase);
        }
        
        std::cout << "Initialized " << n << " particles" << std::endl;
        writeSingleDistribution(STARTDIST_FILE);
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
        //      ∆energy,phase,h		p1
        //      ...					p2
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

	//
	// SIMULATION
	//
    void takeTimestep(int stepID) 
    {
        CalcOp op(stepID, mAcc, mEnergy, mPhase, mCollHits);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size()), op);
    }

    struct ParticleStats { T emax, emin, phmax, phmin; int pleft; };
    ParticleStats getStats() const
    {
        T emax, emin, phmax, phmin;
        emin = phmin = std::numeric_limits<double>::max();
        emax = phmax = -emin;
        int pleft = 0;
        for (int i = 0; i < size(); ++i) {
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
    
    void takeTimesteps(int n, std::string filePath = "", int saveFreq = 1)
    {
        takeTimestepsFromTo(0, n, filePath, saveFreq);
    }

    void takeTimestepsFromTo(int from, int to, std::string filePath = "", int saveFreq = 1)
    {
        if (to <= from)
            throw std::runtime_error("the total number of steps must be > 0");
        int n = to - from;
        
        if (filePath.empty())
            std::cout << "Will not save particle path data" << std::endl;
        else  {
            createTurnFileHeader(filePath, 1 + n/saveFreq); // +1 as we are saving the first starting configuration as well
            writeDistribution(filePath);
        }

        std::vector<T> E_ramp;
        std::vector<std::pair<T, T>> ext_collimators;
        if (mType > NO_RAMP) {
			readRamp<T>(to, E_ramp, mType);
			readCollimators<T>(to, ext_collimators, E_ramp);
            //loadRamp(to , E_ramp, ext_collimators);
            mAcc.setE(E_ramp[from], true);
        }

        std::cout << "Tracking " << size() << " particles for " << n << " turns" << std::endl;
        std::cout << "Starting simulation..." << std::endl;

		common::SilentTimer timer;

        timer.start();
        for (int i = from; i < to; ++i) {
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

            takeTimestep(i);

            int itaken = i - from;
            // reduce the memory footprint a little if it's not necessary
            if (!filePath.empty() && (itaken + 1) % saveFreq == 0) writeDistribution(filePath);

            const int d = n/10;
            if (i % d == 0) {
                int percent = 10*itaken/d;
                //std::cout << "\t" << std::setw(9) << itaken << " of " << n << " turns (" << percent << "%)" << std::endl;
                ParticleStats stats = getStats();
                std::cout << "\t" << std::setw(9) << itaken << " of " << n << " turns (" << percent << "%).";

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

    size_t size() const { return mEnergy.size(); }
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
            const T B2_s = 1 - std::pow(cnst::m_proton/mAcc.E(), 2);
            for (size_t n = range.begin(), N = range.end(); n < N; ++n) {
                if (particleInactive(n)) 
                    continue;

                // We changed the reference energy before, update the ∆E accordingly
                T deltaRef = mAcc.E() - mAcc.E_prev();
                mEnergy[n] -= deltaRef;

                // Give the particle a kick
                // mEnergy[n] -= mAcc.V_rf*std::sin(deltaRef/mAcc.V_rf);

                mEnergy[n] += e*mAcc.V_rf*std::sin(mPhase[n]);
                T B2 = T(1) - std::pow(cnst::m_proton/(mAcc.E() + mEnergy[n]), 2);
                T eta = T(1) - B2 - mAcc.m_compaction;
                mPhase[n] -= T(2)*cnst::pi*mAcc.h_rf*eta/(B2_s*mAcc.E())*mEnergy[n];
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



}; // namespace twodsynch
