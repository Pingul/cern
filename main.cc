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

namespace CONST {

constexpr double pi = 3.14159265359;
constexpr double c = 299792458.0;
constexpr double m_proton = 938e6;

}; // namespace CONST

namespace jwc {

static constexpr double FRAME_X_LOW = -2*CONST::pi;
static constexpr double FRAME_X_HIGH = 4*CONST::pi;
static constexpr double FRAME_Y_LOW = -2e9;
static constexpr double FRAME_Y_HIGH = 2e9;

static constexpr const char* PATH_FILE = "data/particles.dat";
static constexpr const char* LINE_FILE = "data/lines.dat";
static constexpr const char* COLL_FILE = "data/coll.dat";
// static constexpr const char* RAMP_FILE = "resources/ramp.txt";
static constexpr const char* RAMP_FILE = "resources/cLHC_momentum_programme6.5TeV.dat";
static constexpr const char* COLL_MOTOR_FILE = "resources/motor_tcp.txt";
static constexpr const char* STARTDIST_FILE = "data/startdist.dat";

template <typename T>
struct Accelerator
{
	using ValType = T;
	using Acc = Accelerator<T>;

	T E_ref;
	T rf_voltage;
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
		acc.E_ref = T(450e9);
		acc.rf_voltage = T(8e6);
		acc.harmonic = T(35640);
		acc.m_compaction = T(0.0003225);
		acc.coll_top = T(0.5e9);
		acc.coll_bot = T(-0.5e9);
		acc.revolution_freq = T(11245);
		acc.w_revolution_freq = 2*CONST::pi*acc.revolution_freq;
		return acc;
	}

	static Acc getLHC_NOCOLL()
	{
		Acc a = getLHC();
		a.coll_top = std::numeric_limits<T>::max();
		a.coll_bot = std::numeric_limits<T>::min();;
		return a;
	}
};

template <typename T>
inline T hamiltonian(const Accelerator<T>& acc, T deltaE, T phase)
{
	const T rev = acc.w_revolution_freq;
	const T gamma = (acc.E_ref + deltaE)/CONST::m_proton;
	const T gamma_2 = T(1)/(gamma*gamma);
	const T eta = gamma_2 - acc.m_compaction;
	const T beta2 = T(1) - gamma_2;
	const T beta = std::sqrt(beta2);
	const T k = acc.harmonic*rev/(beta*CONST::c);
	const T Omega2 = rev*rev*acc.harmonic*eta*acc.rf_voltage/(T(2)*CONST::pi*beta*acc.E_ref);

	const T H = T(1)/2*beta2*pow(CONST::c*k*eta*deltaE/acc.E_ref, 2) - Omega2*std::cos(phase);
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
			// const T H = hamiltonian(mAcc, deltaE, phase);
			// if (H > 0) {
				mEnergy.push_back(deltaE); 
				mPhase.push_back(phase); 
			// }
		}
	}

	// For parameter passing in the next constructors
	struct LineSim {};
	struct LossAnalysis {};

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
		const T Estep = (Emax - Emin)/300;
		const T PHmax = 2.1*CONST::pi;
		const T PHmin = 0;
		const T PHstep = (PHmax - PHmin)/300;

		for (T de = Emin; de < Emax; de += Estep) {
			for (T ph = PHmin; ph < PHmax; ph += PHstep) {
				const T H = hamiltonian(mAcc, de, ph);
				if (H > 1e5 && H < 2.2e5) {
					mEnergy.push_back(de);
					mPhase.push_back(ph);
				}
			}
		}
		mCollHits.assign(size(), -1);

		std::cout << "Generated " << size() << " particles" << std::endl;
	}

	size_t size() const { return mEnergy.size(); }

	void createTurnFile(std::string filePath, int turns) 
	{
		std::ofstream file(filePath.c_str());
		if (file.is_open()) {
			std::cout << "Saving turns data in '" << filePath << "'" << std::endl;
			file << size() << "," << turns << std::endl; 
		} else 
			std::cerr << "could not write headers" << std::endl;
		file.close();
	}

	void saveParticleCoords(std::string filePath) const 
	{
		std::ofstream file(filePath.c_str(), std::ios::app);

		if (!file.is_open()) 
			throw std::runtime_error("could not save particle coordinates");

		// FILE:
		// 	∆energy / phase 
		for (size_t i = 0; i < size(); ++i) {
			std::stringstream ss;
			ss << mEnergy[i] << "," << mPhase[i] << std::endl;
			file << ss.str();
		}
		file.close();
	}

	void saveCollHits(std::string filePath) const
	{
		std::cout << "Saving coll hits to '" << filePath << "'" << std::endl;

		std::ofstream file(filePath.c_str());
		if (!file.is_open())
			throw std::runtime_error("could not open file");
		else if (mCollHits.empty())
			throw std::runtime_error("could not save coll hits -- no data was recorded");

		// FILE:
		// 	id / turn_lost / phase_lost / energy_lost
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
				for (int i = 0; i < steps; ++i) {
					ramp_file >> data;
					ramp_file >> data;
					E_ramp.push_back(data*1e6);
				}
				ramp_file.close();

				std::ifstream coll_file(COLL_MOTOR_FILE);
				std::cout << "Reading '" << COLL_MOTOR_FILE << "'..." << std::endl;
				for (int i = 0; i < steps; ++i) {
					T bot_coll, top_coll;
					// 		 	 line #
					coll_file >> bot_coll >> bot_coll >> top_coll;
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
			createTurnFile(filePath, n + 1); // +1 as we are saving the first starting configuration as well
			saveParticleCoords(filePath);
		}

		std::vector<T> E_ramp;
		std::vector<std::pair<T, T>> ext_collimators;
		if (mType > NO_RAMP) {
			loadRamp(n, E_ramp, ext_collimators);
			mAcc.E_ref = E_ramp[0];
		}

		std::cout << "Starting simulation..." << std::endl;

		for (int i = 0; i < n; ++i) {
			if (mType > NO_RAMP) {
				T deltaE = E_ramp[i] - mAcc.E_ref;
				mAcc.E_ref = E_ramp[i];

				// Adjust all particle energies
				for (T& e : mEnergy) 
					e -= deltaE;

				if (!ext_collimators.empty()) {
					mAcc.coll_bot = ext_collimators[i].first;
					mAcc.coll_top = ext_collimators[i].second;
				}
			}

			takeTimestep(i);

			// reduce the memory footprint a little if it's not necessary
			if (!filePath.empty() && (i + 1) % saveFreq == 0) saveParticleCoords(filePath);
		}
	}

	std::vector<T> getEnergy() const { return mEnergy; }
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
			const T B2_s = 1 - std::pow(CONST::m_proton/mAcc.E_ref, 2);
			for (size_t n = range.begin(), N = range.end(); n < N; ++n) {
				if (particleInactive(n)) 
					continue;

				mEnergy[n] += e*mAcc.rf_voltage*(std::sin(mPhase[n]) - std::sin(0));
				T B2 = T(1) - std::pow(CONST::m_proton/(mAcc.E_ref + mEnergy[n]), 2);
				T eta = T(1) - B2 - mAcc.m_compaction;
				mPhase[n] += - T(2)*CONST::pi*mAcc.harmonic*eta/(B2_s*mAcc.E_ref)*mEnergy[n];
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
	std::vector<T> mEnergy;
	std::vector<T> mPhase;

	std::vector<int> mCollHits;

	RAMP_TYPE mType;
};

template <typename TModel>
void generateLossmap(typename TModel::RAMP_TYPE type)
{
	TModel tm(type, typename TModel::LossAnalysis());
	tm.saveParticleCoords(STARTDIST_FILE);
	auto freq = TModel::Accelerator::getLHC().revolution_freq;
	tm.takeTimesteps(300*freq); // 50 seconds
	tm.saveCollHits(COLL_FILE);

	int ilastHit, lastHit = -1;
	auto hits = tm.getCollHits();
	for (int i = 0; i < hits.size(); ++i) {
		if (hits[i] > lastHit) {
			lastHit = hits[i];
			ilastHit = i;
		}
	}

	auto energy = tm.getEnergy();

	std::cout 
		<< "Latest hit:\n\tparticle " << ilastHit << ", turn " << lastHit 
		<< "(approx. after " << std::setprecision(3) << (double(lastHit)/freq) << " s)\n"
		<< "\tstarting energy " << std::setprecision(5) << energy[ilastHit] << std::endl;
}

}; // namespace jwc

int main(int argc, char* argv[])
{
	std::vector<std::string> args(argv, argv + argc);

	tbb::task_scheduler_init init(4);

	typedef jwc::ToyModel<double> ToyModel;

	// CHANGE FOR DIFFERENT SIMULATIONS
	ToyModel::RAMP_TYPE type = ToyModel::LHC_RAMP;

	if (args.size() < 2) {
		std::cout << "Not enough arguments specified" << std::endl;
	} else if (args[1] == "animate") {
		{
			std::cout << "Creating particle paths" << std::endl;
			ToyModel tm(1000, type);
			// ToyModel tm(type, ToyModel::LossAnalysis());
			tm.takeTimesteps(50000, jwc::PATH_FILE, 10);
			tm.saveCollHits(jwc::COLL_FILE);
		}
		{
			std::cout << "Creating line segments" << std::endl;
			ToyModel tm(type, ToyModel::LineSim());
			tm.takeTimesteps(1700, jwc::LINE_FILE, 5);
		}
	} else if (args[1] == "lossmap") {
		std::cout << "Creating data for lossmap" << std::endl;
		ToyModel tm(10000, type);
		tm.takeTimesteps(10000);
		tm.saveCollHits(jwc::COLL_FILE);
	} else if (args[1] == "energy") {
		std::cout << "Simulating 1 particle" << std::endl;
		ToyModel tm(1, type);
		tm.takeTimesteps(40000, jwc::PATH_FILE);
		tm.saveCollHits(jwc::COLL_FILE);
	} else if (args[1] == "lossmap-analysis") {
		std::cout << "Loss pattern analysis" << std::endl;
		jwc::generateLossmap<ToyModel>(type);
	} else {
		std::cout << "No action with name '" << args[1] << "' found" << std::endl;
	}

	return 0;
}