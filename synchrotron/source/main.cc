#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <tbb/task_scheduler_init.h>

#include "synchrotron.hh"
#include "accelerator.hh"
#include "settings.hh"
#include "hamiltonian.hh"
#include "particles.hh"
#include "ramp_program.hh"

void generatePhasespaceLines(int seconds)
{
    using namespace stron;

    // will generate 1 per second
    std::cout << "Generate phasespace lines" << std::endl;
    auto acc = Accelerator<double>::getLHC();
    int freq = int(acc.f_rev);
    int turns = seconds*freq;

    //std::vector<double> E;
    //readRamp(turns, E, LHC_RAMP);
    auto program = stron::ramp::create(acc, turns, stron::ramp::ProgramType::LHCRamp);

    for (int i = 0; i < seconds; ++i) {
        int turn = i*freq;
        for (int j = 0; j < freq; ++j) program->step();
        
        // Caluclated from LHC_ramp.dat
        const double k = 2.9491187074838457087e-07;
        acc.V_rf = (6 + k*turn)*1e6;
        
        if (turn % 11245 == 0)
            std::cout << i << " " << acc.V_rf << std::endl;
        //std::stringstream ss;
        //ss << "phasespace/" << i << "lines.dat";
        //writePhasespaceFrame(acc, ss.str());
    }
}


int main(int argc, char* argv[])
{
    std::vector<std::string> args(argv, argv + argc);

    tbb::task_scheduler_init init;

    using SimpleSynchrotron = stron::SimpleSynchrotron<double>;
    using Acc = typename SimpleSynchrotron::Acc;
    using ProgramType = stron::ramp::ProgramType;

    //auto LHC = SimpleSynchrotron::Acc::getLHC();

    // CHANGE FOR DIFFERENT SIMULATIONS
    ProgramType progType = ProgramType::AggressiveRamp;

    if (args.size() < 2) {
        std::cout << "Not enough arguments specified" << std::endl;
    
    } else if (args[1] == "lossmap" || args[1] == "startdist") {
        // We often work with these two together, so we make sure we have the same
        // SimpleSynchrotron for both of these
        //SimpleSynchrotron ss(2500, type, SimpleSynchrotron::LinearDecay());

        SimpleSynchrotron ss(Acc::getLHC());
        ss.addParticles(stron::pdist::ActionValues<double>(2500, ss.getAcc()));
        //SimpleSynchrotron ss(250, type, SimpleSynchrotron::ActionValues());
        if (args[1] == "lossmap")
            ss.runLossmap(stron::ramp::create(ss.getAcc(), 50*11245, progType));

    } else if (args[1].find("animate") == 0) {
        SimpleSynchrotron ss(Acc::getLHC());
        ss.addParticles(stron::pdist::AroundSeparatrix<double>(500, ss.getAcc()));
        if (args[1] == "animate") {
            ss.simulateTurns(stron::ramp::create(ss.getAcc(), 1000, progType), stron::PATH_FILE, 2);
        } else if (args[1] == "animate-long") {
            ss.simulateTurns(stron::ramp::create(ss.getAcc(), 20*11245, progType), stron::PATH_FILE, 1000);
        } else if (args[1] == "animate-background") {
            ss.simulateTurns(stron::ramp::create(ss.getAcc(), 300*11245, progType), stron::PATH_FILE, 11245);
            generatePhasespaceLines(300);
        }
        ss.writeCollHits(stron::COLL_FILE);
        writePhasespaceFrame(Acc::getLHC(), stron::LINE_FILE);

    } else if (args[1] == "sixtrack-comp") {
        std::cout << "Sixtrack comparison" << std::endl;
        double momentum = 0.0;
        if (args.size() == 3)
            momentum = std::stod(args[2])*1e6;
        std::cout << "∆E = " << std::setprecision(16) << momentum << std::endl;

        SimpleSynchrotron ss(Acc::getLHC());
        ss.addParticles(stron::pdist::SixTrackTest<double>(momentum));
        ss.simulateTurns(stron::ramp::create(ss.getAcc(), 224900, progType), stron::SIXTRACK_TEST_FILE);

    } else if (args[1] == "phasespace") {
        writePhasespaceFrame(SimpleSynchrotron::Acc::getLHC(), stron::LINE_FILE);
    } else if (args[1] == "phasespace-mov") {
        generatePhasespaceLines(100);

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
