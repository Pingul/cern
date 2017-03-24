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
    
    } else if (args[1] == "lossmap" || args[1] == "startdist") {
        // We often work with these two together, so we make sure we have the same
        // ToyModel for both of these
        ToyModel tm(2500, type, ToyModel::LinearDecay());
        if (args[1] == "lossmap")
            tm.runLossmap(50);

    } else if (args[1].find("animate") == 0) {
        ToyModel tm(500, type, ToyModel::AroundSeparatrix());
        if (args[1] == "animate") {
            tm.simulateTurns(1000, twodsynch::PATH_FILE, 1);
        } else if (args[1] == "animate-long") {
            tm.simulateTurns(20*11245, twodsynch::PATH_FILE, 1000);
        } else if (args[1] == "animate-background") {
            tm.simulateTurns(300*11245, twodsynch::PATH_FILE, 11245);
            twodsynch::generatePhasespaceLines(300);
        }
        tm.writeCollHits(twodsynch::COLL_FILE);
        writePhasespaceFrame(ToyModel::Acc::getLHC(), twodsynch::LINE_FILE);

    } else if (args[1] == "sixtrack-comp") {
        std::cout << "Sixtrack comparison" << std::endl;
        double energy = 0.0;
        if (args.size() == 3)
            energy = std::stod(args[2])*1e6;
        std::cout << "∆E = " << std::setprecision(16) << energy << std::endl;
        ToyModel tm(ToyModel::SixTrackTest(), energy );
        tm.simulateTurns(224900, twodsynch::SIXTRACK_TEST_FILE);

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
        }

    } else if (args[1] == "lost") {
        ToyModel tm(1000, type, ToyModel::AroundSeparatrix());
        tm.simulateTurns(40*11245);
        tm.writeCollHits(twodsynch::COLL_FILE);
        tm.writeLostTurns("calc/lost.dat");
    } else if (args[1] == "test") {
        using namespace twodsynch;

        auto acc = Accelerator<double>::getLHC();
        std::cout << acc.lag_phase() << std::endl;
        const double Hsep = 1e8;
        std::cout << "separatrix H = " << Hsep << std::endl;
        for (double ph = 0; ph < 2.0*cnst::pi; ph += 0.2*cnst::pi) {
            const double de = levelCurve(acc, ph, Hsep);
            const double H = hamiltonian(acc, de, ph);
            std::cout << "H(" << ph << ", " << de << ") = " << H << std::endl;
        }
    } else if (args[1] == "test2") {
        using namespace twodsynch;

        auto acc = Accelerator<double>::getLHC();
        const double de = -1.90986e+08;
        std::cout << "de = " << de << std::endl;
        for (double ph = 0; ph < 2.0*cnst::pi; ph += 0.2*cnst::pi) {
            const double H = hamiltonian(acc, de, ph);
            const double dde = levelCurve(acc, ph, H);
            std::cout << "H-1(" << ph << ", " << H << ") = " << dde << std::endl;
        }
    } else if (args[1] == "f_rf") {
        std::cout << std::setprecision(16) << twodsynch::Accelerator<double>::getLHC().f_rf << std::endl;
    } else {
        std::cout << "No action with name '" << args[1] << "' found" << std::endl;
    }

    return 0;
}
