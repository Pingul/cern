#include <iostream>
#include <vector>
#include <string>
#include <tbb/task_scheduler_init.h>
#include "2d-synchrotron.hh"

namespace jwc {

template <typename TModel>
void generateLossmap(typename TModel::RAMP_TYPE type)
{
	TModel tm(type, typename TModel::LossAnalysis());
	tm.createTurnFileHeader(STARTDIST_FILE, 1);
	tm.writeDistribution(STARTDIST_FILE);
	auto freq = TModel::Accelerator::getLHC().revolution_freq;
	int turns = 300*freq;
	std::cout << "Simulating " << turns << " turns" << std::endl;
	tm.takeTimesteps(turns); // 50 seconds
	tm.writeCollHits(COLL_FILE);

	int ilastHit, tlastHit = -1;
	auto hits = tm.getCollHits();
	for (int i = 0; i < hits.size(); ++i) {
		if (hits[i] > tlastHit) {
			tlastHit = hits[i];
			ilastHit = i;
		}
	}

	if (tlastHit == -1)
		std::cout << "No losses" << std::endl;
	else {
		auto energy = tm.getEnergy();
		std::cout 
			<< "Latest hit:\n\tparticle " << ilastHit << ", turn " << tlastHit 
			<< "(approx. after " << std::setprecision(3) << (double(tlastHit)/freq) << " s)\n"
			<< "\tstarting energy " << std::setprecision(5) << energy[ilastHit] << std::endl;
	}

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
			tm.takeTimesteps(5000, jwc::PATH_FILE, 10);
			tm.writeCollHits(jwc::COLL_FILE);
		}
		{
			std::cout << "Creating line segments" << std::endl;
			ToyModel tm(type, ToyModel::LineSim());
			tm.takeTimesteps(1700, jwc::LINE_FILE, 5);
		}
	} else if (args[1] == "energy") {
		std::cout << "Simulating 1 particle" << std::endl;
		ToyModel tm(1, type);
		tm.takeTimesteps(40000, jwc::PATH_FILE);
		tm.writeCollHits(jwc::COLL_FILE);
	} else if (args[1] == "lossmap-analysis" || args[1] == "lossmap") {
		std::cout << "Loss pattern analysis" << std::endl;
		jwc::generateLossmap<ToyModel>(type);
	} else if (args[1] == "sixtrack-comp") {
		std::cout << "Sixtrack comparison" << std::endl;
		ToyModel tm( (ToyModel::SixTrackTest()) );
		tm.takeTimesteps(40000, jwc::SIXTRACK_TEST_FILE);
	} else {
		std::cout << "No action with name '" << args[1] << "' found" << std::endl;
	}

	return 0;
}
