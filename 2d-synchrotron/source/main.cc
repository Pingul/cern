#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <tbb/task_scheduler_init.h>

#include "2d-synchrotron.hh"
#include "accelerator.hh"
#include "settings.hh"
#include "hamiltonian.hh"

void generatePhasespaceLines(int seconds)
{
    using namespace stron;

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


int main(int argc, char* argv[])
{
    std::vector<std::string> args(argv, argv + argc);

    tbb::task_scheduler_init init;

    typedef stron::SimpleSynchrotron<double> SimpleSynchrotron;

    // CHANGE FOR DIFFERENT SIMULATIONS
    stron::RAMP_TYPE type = stron::LHC_RAMP;

    if (args.size() < 2) {
        std::cout << "Not enough arguments specified" << std::endl;
    
    } else if (args[1] == "lossmap" || args[1] == "startdist") {
        // We often work with these two together, so we make sure we have the same
        // SimpleSynchrotron for both of these
        //SimpleSynchrotron tm(2500, type, SimpleSynchrotron::LinearDecay());
        SimpleSynchrotron tm(250, type, SimpleSynchrotron::ActionValues());
        if (args[1] == "lossmap")
            tm.runLossmap(50);

    } else if (args[1].find("animate") == 0) {
        //SimpleSynchrotron tm(500, type, SimpleSynchrotron::AroundSeparatrix());
        SimpleSynchrotron tm(500, type);
        if (args[1] == "animate") {
            tm.simulateTurns(1000, stron::PATH_FILE, 2);
        } else if (args[1] == "animate-long") {
            tm.simulateTurns(20*11245, stron::PATH_FILE, 1000);
        } else if (args[1] == "animate-background") {
            tm.simulateTurns(300*11245, stron::PATH_FILE, 11245);
            generatePhasespaceLines(300);
        }
        tm.writeCollHits(stron::COLL_FILE);
        writePhasespaceFrame(SimpleSynchrotron::Acc::getLHC(), stron::LINE_FILE);

    } else if (args[1] == "sixtrack-comp") {
        std::cout << "Sixtrack comparison" << std::endl;
        double energy = 0.0;
        if (args.size() == 3)
            energy = std::stod(args[2])*1e6;
        std::cout << "∆E = " << std::setprecision(16) << energy << std::endl;
        SimpleSynchrotron tm(SimpleSynchrotron::SixTrackTest(), energy );
        tm.simulateTurns(224900, stron::SIXTRACK_TEST_FILE);

    } else if (args[1] == "phasespace") {
        writePhasespaceFrame(SimpleSynchrotron::Acc::getLHC(), stron::LINE_FILE);
    } else if (args[1] == "phasespace-mov") {
        generatePhasespaceLines(300);

    } else if (args[1] == "test") {
        using namespace stron;

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
        using namespace stron;

        auto acc = Accelerator<double>::getLHC();
        const double de = -1.90986e+08;
        std::cout << "de = " << de << std::endl;
        for (double ph = 0; ph < 2.0*cnst::pi; ph += 0.2*cnst::pi) {
            const double H = hamiltonian(acc, de, ph);
            const double dde = levelCurve(acc, ph, H);
            std::cout << "H-1(" << ph << ", " << H << ") = " << dde << std::endl;
        }
    } else if (args[1] == "f_rf") {
        std::cout << std::setprecision(16) << stron::Accelerator<double>::getLHC().f_rf << std::endl;
    } else {
        std::cout << "No action with name '" << args[1] << "' found" << std::endl;
    }

    return 0;
}
