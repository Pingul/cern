#include <iostream>
#include <vector>
#include <string>
#include <tbb/task_scheduler_init.h>
#include "2d-synchrotron.hh"

namespace twodsynch {

template <typename TModel>
void generateLossmap(typename TModel::RAMP_TYPE type)
{
	//TModel tm(2000, type, typename TModel::LossAnalysis2());
	TModel tm("config/highlow.dat", type);
    auto freq = TModel::Accelerator::getLHC().f_rev;
    int turns = 500*freq;

    // hack...
    std::vector<double> startE(tm.getEnergy());
    std::vector<double> startPh(tm.getPhase());
    // ...end of hack

	//tm.takeTimesteps(turns); 
	tm.takeTimesteps(turns, twodsynch::PATH_FILE, freq); 
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
            << "\tstarting H=" << std::setprecision(16) << twodsynch::hamiltonian(TModel::Accelerator::getLHC(),
            startE[ilastHit], startPh[ilastHit]) << std::endl;
    }

}

}; // namespace twodsynch

int main(int argc, char* argv[])
{
    std::vector<std::string> args(argv, argv + argc);

    tbb::task_scheduler_init init;

    typedef twodsynch::ToyModel<double> ToyModel;

    //std::cout << twodsynch::synchrotron_frequency() << std::endl;

    // CHANGE FOR00 DIFFERENT SIMULATIONS
    ToyModel::RAMP_TYPE type = ToyModel::LHC_RAMP;

    if (args.size() < 2) {
        std::cout << "Not enough arguments specified" << std::endl;
    } else if (args[1] == "animate") {
        {
            std::cout << "Creating particle paths" << std::endl;
            //ToyModel tm(1000, type);
            ToyModel tm(500, type, ToyModel::LossAnalysis());
            tm.takeTimesteps(5000, twodsynch::PATH_FILE, 10);
            tm.writeCollHits(twodsynch::COLL_FILE);
        }
        //{
            //std::cout << "Creating line segments" << std::endl;
            //ToyModel tm(type, ToyModel::LineSim());
            //tm.takeTimesteps(1700, twodsynch::LINE_FILE, 5);
        //}
    } else if (args[1] == "energy") {
        std::cout << "Simulating 1 particle" << std::endl;
        ToyModel tm(1, type);
        tm.takeTimesteps(40000, twodsynch::PATH_FILE);
        tm.writeCollHits(twodsynch::COLL_FILE);
    } else if (args[1] == "lossmap-analysis" || args[1] == "lossmap") {
        std::cout << "Loss pattern analysis" << std::endl;
        twodsynch::generateLossmap<ToyModel>(type);
    } else if (args[1] == "sixtrack-comp") {
        std::cout << "Sixtrack comparison" << std::endl;
        ToyModel tm( (ToyModel::SixTrackTest()) );
        tm.takeTimesteps(30000, twodsynch::SIXTRACK_TEST_FILE);
    } else if (args[1] == "startdist") {
        std::cout << "Start distribution" << std::endl;
        ToyModel tm(2000, type, ToyModel::LossAnalysis());
    } else if (args[1] == "phasespace") {
		phasespaceLevelCurve(ToyModel::Accelerator::getLHC(), twodsynch::LINE_FILE);

        //int turn = std::stoi(args[2]);

        //// Just need a small starting distribution
        //{
            //ToyModel tm(1, type, ToyModel::LossAnalysis());
        //}

        //std::cout << "Creating line segments for turn " << turn << std::endl;
        //ToyModel tm(type, ToyModel::LineSim());
        //tm.takeTimestepsFromTo(turn, 1700 + turn, twodsynch::LINE_FILE, 5);
    } else if (args[1] == "restart") {
        if (args.size() < 3)
            std::cout << "Must provide file path" << std::endl;
        else {
            ToyModel tm(args[2], type);
            tm.takeTimesteps(20*11245, twodsynch::PATH_FILE, 4);
        }
	} else if (args[1] == "test") {
		using twodsynch::cnst::pi;

		auto acc = ToyModel::Accelerator::getLHC();
		double ph = pi;
		for (double de = 0; de < 1e8; de += 1e7) {
			double H = hamiltonian(acc, de, ph);
			std::cout << de << " -> " << std::setprecision(32) << H << " -> " << levelCurve(acc, ph, H) << std::endl;
		}
		std::cout << hamiltonian(acc, 0.0, pi) << std::endl;

    } else {
        std::cout << "No action with name '" << args[1] << "' found" << std::endl;
    }

    return 0;
}
