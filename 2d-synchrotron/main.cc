#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <tbb/task_scheduler_init.h>
#include "2d-synchrotron.hh"

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
        tm.takeTimesteps(5000, twodsynch::PATH_FILE, 5);
        tm.writeCollHits(twodsynch::COLL_FILE);
        writePhasespaceFrame(ToyModel::Accelerator::getLHC(), twodsynch::LINE_FILE);
    } else if (args[1] == "animate-long") {
        std::cout << "Long animation" << std::endl;
        ToyModel tm(500, type, ToyModel::LossAnalysis());
        //tm.takeTimesteps(5000, twodsynch::PATH_FILE, 10);
        tm.takeTimesteps(300*11245, twodsynch::PATH_FILE, 11245);
        tm.writeCollHits(twodsynch::COLL_FILE);
        writePhasespaceFrame(ToyModel::Accelerator::getLHC(), twodsynch::LINE_FILE);
    } else if (args[1] == "animate-background") {
        // In most cases unnecessary if the background lines has already been generated once
        int seconds = 300;
        std::cout << "Full animation for " << seconds << " s" << std::endl;
        ToyModel tm(500, type, ToyModel::LossAnalysis());
        tm.takeTimesteps(11245, twodsynch::PATH_FILE, 11245);
        tm.writeCollHits(twodsynch::COLL_FILE);
        twodsynch::generatePhasespaceLines(seconds);
    } else if (args[1] == "energy") {
        std::cout << "Simulating 1 particle" << std::endl;
        ToyModel tm(1, type);
        tm.takeTimesteps(40000, twodsynch::PATH_FILE);
        tm.writeCollHits(twodsynch::COLL_FILE);
    } else if (args[1] == "lossmap-analysis" || args[1] == "lossmap") {
        std::cout << "Loss pattern analysis" << std::endl;
        //ToyModel lossmodel(1000, type, ToyModel::LossAnalysis());
        ToyModel lossmodel("test_lost.dat", type, ToyModel::LossAnalysisMultiplied());
        lossmodel.runLossmap(30);
    } else if (args[1] == "sixtrack-comp") {
        std::cout << "Sixtrack comparison" << std::endl;
        ToyModel tm( (ToyModel::SixTrackTest()) );
        tm.takeTimesteps(30000, twodsynch::SIXTRACK_TEST_FILE);
    } else if (args[1] == "startdist") {
        std::cout << "Start distribution" << std::endl;
        ToyModel lossmodel("test_lost.dat", type, ToyModel::LossAnalysisMultiplied());
        //ToyModel tm(2000, type, ToyModel::LossAnalysis());
    } else if (args[1] == "phasespace") {
        writePhasespaceFrame(ToyModel::Accelerator::getLHC(), twodsynch::LINE_FILE);
    } else if (args[1] == "phasespace-mov") {
        twodsynch::generatePhasespaceLines(300);
    } else if (args[1] == "restart") {
        if (args.size() < 3)
            std::cout << "Must provide file path" << std::endl;
        else {
            ToyModel tm(args[2], type);
            tm.takeTimesteps(20*11245, twodsynch::PATH_FILE, 4);
        }
    } else if (args[1] == "lost") {
        ToyModel tm(1000, type, ToyModel::LossAnalysis());
        tm.takeTimesteps(40*11245);
        tm.writeCollHits(twodsynch::COLL_FILE);
        tm.writeLostTurns("calc/lost.dat");
    } else if (args[1] == "test") {
        auto acc = twodsynch::Accelerator<double>::getLHC();
        acc.setE(7000e9, true);
        std::cout << "Default rev frequency: " << acc.f_rev << std::endl;
    } else {
        std::cout << "No action with name '" << args[1] << "' found" << std::endl;
    }

    return 0;
}
