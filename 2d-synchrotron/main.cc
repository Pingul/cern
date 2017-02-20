#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <tbb/task_scheduler_init.h>
#include "2d-synchrotron.hh"

namespace twodsynch {

template <typename TModel>
void generateLossmap(RAMP_TYPE type)
{
	TModel tm(2000, type, typename TModel::LossAnalysis());
	//TModel tm("config/2p.dat", type);
    auto freq = TModel::Accelerator::getLHC().f_rev;
    int turns = 50*freq;

    // hack...
    std::vector<double> startE(tm.getEnergy());
    std::vector<double> startPh(tm.getPhase());
    // ...end of hack

	tm.takeTimesteps(turns); 
	//tm.takeTimesteps(turns, twodsynch::PATH_FILE, 100); 
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

void generatePhasespaceLines()
{
	std::cout << "Generate phasespace lines" << std::endl;
	auto acc = Accelerator<double>::getLHC();
	int freq = int(acc.f_rev);
	int seconds = 300;
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

}; // namespace twodsynch

int main(int argc, char* argv[])
{
    std::vector<std::string> args(argv, argv + argc);

    tbb::task_scheduler_init init;

    typedef twodsynch::ToyModel<double> ToyModel;

    // CHANGE FOR DIFFERENT SIMULATIONS
    twodsynch::RAMP_TYPE type = twodsynch::LHC_RAMP;

    if (args.size() < 2) {
        std::cout << "Not enough arguments specified" << std::endl;
    } else if (args[1] == "animate") {
        std::cout << "Short animation" << std::endl;
        ToyModel tm(500, type, ToyModel::LossAnalysis());
        //tm.takeTimesteps(5000, twodsynch::PATH_FILE, 10);
        tm.takeTimesteps(5000, twodsynch::PATH_FILE, 5, false);
        tm.writeCollHits(twodsynch::COLL_FILE);
		writePhasespaceFrame(ToyModel::Accelerator::getLHC(), twodsynch::LINE_FILE);
	} else if (args[1] == "animate-long") {
        std::cout << "Long animation" << std::endl;
        ToyModel tm(500, type, ToyModel::LossAnalysis());
        //tm.takeTimesteps(5000, twodsynch::PATH_FILE, 10);
        tm.takeTimesteps(300*11245, twodsynch::PATH_FILE, 11245, false);
        tm.writeCollHits(twodsynch::COLL_FILE);
		writePhasespaceFrame(ToyModel::Accelerator::getLHC(), twodsynch::LINE_FILE);
	} else if (args[1] == "animate-background") {
		// In most cases unnecessary if the background lines has already been generated once
		std::cout << "Full animation" << std::endl;
		ToyModel tm(500, type, ToyModel::LossAnalysis());
		tm.takeTimesteps(300*11245, twodsynch::PATH_FILE, 11245, true);
		tm.writeCollHits(twodsynch::COLL_FILE);
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
		writePhasespaceFrame(ToyModel::Accelerator::getLHC(), twodsynch::LINE_FILE);
	} else if (args[1] == "phasespace-mov") {
		twodsynch::generatePhasespaceLines();
    } else if (args[1] == "restart") {
        if (args.size() < 3)
            std::cout << "Must provide file path" << std::endl;
        else {
            ToyModel tm(args[2], type);
            tm.takeTimesteps(20*11245, twodsynch::PATH_FILE, 4);
        }
	} else if (args[1] == "lost") {
		ToyModel tm(1000, type, ToyModel::LossAnalysis());
		tm.takeTimesteps(60*11245);
		tm.writeCollHits(twodsynch::COLL_FILE);
		tm.writeLostTurns("calc/lost.dat");
    } else {
        std::cout << "No action with name '" << args[1] << "' found" << std::endl;
    }

    return 0;
}
