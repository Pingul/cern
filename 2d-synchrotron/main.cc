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
        //ToyModel tm(10, type, ToyModel::LossAnalysis());
        ToyModel tm(type, ToyModel::LossAnalysisAction2());
        tm.simulateTurns(1000, twodsynch::PATH_FILE, 1);
        tm.writeCollHits(twodsynch::COLL_FILE);
        writePhasespaceFrame(ToyModel::Acc::getLHC(), twodsynch::LINE_FILE);
    } else if (args[1] == "animate-long") {
        std::cout << "Long animation" << std::endl;
        ToyModel tm(500, type, ToyModel::LossAnalysis());
        //tm.simulateTurns(5000, twodsynch::PATH_FILE, 10);
        tm.simulateTurns(300*11245, twodsynch::PATH_FILE, 11245);
        tm.writeCollHits(twodsynch::COLL_FILE);
        writePhasespaceFrame(ToyModel::Acc::getLHC(), twodsynch::LINE_FILE);
    } else if (args[1] == "animate-background") {
        // In most cases unnecessary if the background lines has already been generated once
        int seconds = 300;
        std::cout << "Full animation for " << seconds << " s" << std::endl;
        ToyModel tm(500, type, ToyModel::LossAnalysis());
        tm.simulateTurns(11245, twodsynch::PATH_FILE, 11245);
        tm.writeCollHits(twodsynch::COLL_FILE);
        twodsynch::generatePhasespaceLines(seconds);
    } else if (args[1] == "energy") {
        std::cout << "Simulating 1 particle" << std::endl;
        ToyModel tm(1, type);
        tm.simulateTurns(40000, twodsynch::PATH_FILE);
        tm.writeCollHits(twodsynch::COLL_FILE);
    } else if (args[1] == "lossmap-analysis" || args[1] == "lossmap") {
        std::cout << "Loss pattern analysis" << std::endl;
        ToyModel lossmodel(type, ToyModel::LossAnalysisAction());
        //ToyModel lossmodel(500, type, ToyModel::LossAnalysis());
        lossmodel.runLossmap(50);
    } else if (args[1] == "sixtrack-comp") {
        std::cout << "Sixtrack comparison" << std::endl;
        ToyModel tm( (ToyModel::SixTrackTest()) );
        tm.simulateTurns(30000, twodsynch::SIXTRACK_TEST_FILE);
    } else if (args[1] == "startdist") {
        std::cout << "Start distribution" << std::endl;
        ToyModel lossmodel(type, ToyModel::LossAnalysisAction());
        //ToyModel tm(2000, type, ToyModel::LossAnalysis());
    } else if (args[1] == "phasespace") {
        writePhasespaceFrame(ToyModel::Acc::getLHC(), twodsynch::LINE_FILE);
    } else if (args[1] == "phasespace-mov") {
        twodsynch::generatePhasespaceLines(300);
    } else if (args[1] == "restart") {
        if (args.size() < 3)
            std::cout << "Must provide file path" << std::endl;
        else {
            ToyModel tm(args[2], type);
            tm.simulateTurns(600, twodsynch::PATH_FILE, 1);
            //tm.runLossmap(20);
            //tm.simulateTurns(20*11245, twodsynch::PATH_FILE, 4);
        }
    } else if (args[1] == "lost") {
        ToyModel tm(1000, type, ToyModel::LossAnalysis());
        tm.simulateTurns(40*11245);
        tm.writeCollHits(twodsynch::COLL_FILE);
        tm.writeLostTurns("calc/lost.dat");
    } else if (args[1] == "test") {
        using twodsynch::cnst::pi;
        using twodsynch::levelCurve;
        using twodsynch::hamiltonian;

        auto acc = twodsynch::Accelerator<double>::getLHC();
        std::cout << acc.lag_phase() << std::endl;
        //const double Hsep = -twodsynch::separatrix(acc);
        const double Hsep = 1e8;
        std::cout << "separatrix H = " << Hsep << std::endl;
        for (double ph = 0; ph < 2.0*pi; ph += 0.2*pi) {
            const double de = twodsynch::levelCurve(acc, ph, Hsep);
            const double H = twodsynch::hamiltonian(acc, de, ph);
            std::cout << "H(" << ph << ", " << de << ") = " << H << std::endl;
        }
    } else if (args[1] == "test2") {
        using twodsynch::cnst::pi;
        using twodsynch::levelCurve;
        using twodsynch::hamiltonian;


        auto acc = twodsynch::Accelerator<double>::getLHC();
        const double de = -1.90986e+08;
        std::cout << "de = " << de << std::endl;
        for (double ph = 0; ph < 2.0*pi; ph += 0.2*pi) {
            const double H = twodsynch::hamiltonian(acc, de, ph);
            const double dde = twodsynch::levelCurve(acc, ph, H);
            std::cout << "H-1(" << ph << ", " << H << ") = " << dde << std::endl;
        }
    } else {
        std::cout << "No action with name '" << args[1] << "' found" << std::endl;
    }

    return 0;
}
