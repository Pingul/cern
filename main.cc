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


namespace jwc {

static constexpr double PI = 3.14159265359;

static constexpr double FRAME_X_LOW = -2*PI;
static constexpr double FRAME_X_HIGH = 4*PI;
static constexpr double FRAME_Y_LOW = -2e9;
static constexpr double FRAME_Y_HIGH = 2e9;

static constexpr double E_REF = 450e9;

static constexpr int FREQ = 11245;

static constexpr const char* PATH_FILE = "particles.dat";
static constexpr const char* LINE_FILE = "line.dat";
static constexpr const char* COLL_FILE = "coll.dat";
static constexpr const char* RAMP_FILE = "ramp.txt";
static constexpr const char* STARTDIST_FILE = "startdist.dat";

template <typename T>
struct Accelerator
{
	using ValType = T;
	using Acc = Accelerator<T>;

	T rf_voltage;
	T harmonic;
	T m_compaction;
	T coll_top;
	T coll_bot;

	static Acc getLHC() 
	{
		// Parameters for LHC can be found at http://atlas.physics.arizona.edu/~kjohns/downloads/linac/UTclass2lhc.pdf
		return Acc{
			/*rf_voltage=*/T(8e6),
			/*harmonic=*/T(35640),
			/*m_compaction=*/T(0.0003225),
			/*coll_top=*/T(0.5e9),
			/*coll_bot=*/T(-0.5e9)
		};
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
		// std::normal_distribution<> ph_dist(PI, PI);

		std::normal_distribution<> e_dist(0, 0.15e9);
		std::normal_distribution<> ph_dist(PI, PI/5);
		for (size_t i = 0; i < n; ++i) {
			mEnergy.push_back(e_dist(generator)); 
			mPhase.push_back(ph_dist(generator)); 
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
				for (T x = FRAME_X_LOW; x < FRAME_X_HIGH; x += PI/4.0) {
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
				for (T x = -PI; x < FRAME_X_HIGH; x += 2*PI) {
					T d = PI/2;
					for (T delta = -d; delta < d; delta += PI/5) {
						mEnergy.push_back(0.0f);
						mPhase.push_back(x + delta);
					}
				}
				break;
		}
	}

	ToyModel(size_t n, T maxDE, T minDE, RAMP_TYPE type, LossAnalysis) 
		: mAcc(Accelerator::getLHC()), mType(type)
	{
		maxDE = std::abs(maxDE);
		minDE = std::abs(minDE);

		if (! (maxDE > minDE)) 
			throw std::runtime_error("maxDE needs to be bigger than minDE");

		mEnergy.reserve(n);
		mPhase.reserve(n);
		mCollHits.assign(n, -1);

		for (T de = minDE; de < maxDE; de += 2.0*(maxDE - minDE)/n) {
			mEnergy.push_back(de);
			mPhase.push_back(PI);

			mEnergy.push_back(-de);
			mPhase.push_back(PI);
		}
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

	void takeTimestep(int stepID, T E_ref) 
	{
		CalcOp op(stepID, mAcc, E_ref, mEnergy, mPhase, mCollHits);
		tbb::parallel_for(tbb::blocked_range<size_t>(0, size()), op);
	}

	void takeTimesteps(int n, std::string filePath = "", int saveFreq = 1)
	{
		if (filePath.empty())
			std::cout << "Will not save particle path data" << std::endl;
		else  {
			std::cout << "Saving path to '" << filePath << "'" << std::endl;
			createTurnFile(filePath, n + 1); // +1 as we are saving the first starting configuration as well
			saveParticleCoords(filePath);
		}

		T E_ref = E_REF;
		std::vector<T> E_ramp;
		if (mType > NO_RAMP) {
			E_ramp.reserve(n);

			std::ifstream file(RAMP_FILE);
			T data;
			for (int i = 0; i < n; ++i) {
				switch (mType) {
					default:
					case LHC_RAMP:
						file >> data; // the line number
						file >> data; // the actual data
						E_ramp.push_back(data*1e6);
						break;
					case AGGRESSIVE_RAMP:
						// linear ramp
						E_ramp.push_back(E_ref + i*1e7);
						break;
					case SEMI_AGGRESSIVE_RAMP:
						E_ramp.push_back(E_ref + i*3e6);
						break;
				}

			}
			E_ref = E_ramp[0];
		}

		for (int i = 0; i < n; ++i) {
			if (mType > NO_RAMP) {
				T deltaE = E_ramp[i] - E_ref;
				E_ref = E_ramp[i];

				// Adjust all particle energies
				for (T& e : mEnergy) e -= deltaE;
			}

			takeTimestep(i, E_ref);

			// reduce the memory footprint a little if it's not necessary
			if (!filePath.empty() && (i + 1) % saveFreq == 0) saveParticleCoords(filePath);
		}
	}

	std::vector<T> getEnergy() const { return mEnergy; }
	std::vector<int> getCollHits() const { return mCollHits; }

private:
	struct CalcOp
	{
		CalcOp(int stepID, const Accelerator& acc, T E_ref, std::vector<T>& energy, 
			   std::vector<T>& phase, std::vector<int>& collHits)
			: mStepID(stepID), mAcc(acc), mE_ref(E_ref), mEnergy(energy), mPhase(phase), mCollHits(collHits)
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
			const T m = 938e6;
			const T B2_s = 1 - std::pow(m/mE_ref, 2);
			for (size_t n = range.begin(), N = range.end(); n < N; ++n) {
				if (particleInactive(n)) 
					continue;

				mEnergy[n] += e*mAcc.rf_voltage*(std::sin(mPhase[n]) - std::sin(0));
				T B2 = T(1) - std::pow(m/(mE_ref + mEnergy[n]), 2);
				T eta = T(1) - B2 - mAcc.m_compaction;
				mPhase[n] += - T(2)*PI*mAcc.harmonic*eta/(B2_s*mE_ref)*mEnergy[n];
				if (particleCollided(n)) mCollHits[n] = mStepID;
			}
		}
	private:
		int mStepID;

		const Accelerator& mAcc;
		const T mE_ref;

		std::vector<T>& mEnergy;
		std::vector<T>& mPhase;

		std::vector<int>& mCollHits;
	};

	int mStepID;

	const Accelerator mAcc;
	// Particle properties
	std::vector<T> mEnergy;
	std::vector<T> mPhase;

	std::vector<int> mCollHits;

	RAMP_TYPE mType;
};

template <typename TModel>
void generateLossmap(typename TModel::RAMP_TYPE type)
{
	TModel tm(16000, 0.45e9, 0.44e9, type, typename TModel::LossAnalysis());
	tm.saveParticleCoords(STARTDIST_FILE);
	tm.takeTimesteps(50*FREQ); // 50 seconds
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
		<< "(approx. after " << std::setprecision(3) << (double(lastHit)/FREQ) << " s)\n"
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
			tm.takeTimesteps(500, jwc::PATH_FILE);
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