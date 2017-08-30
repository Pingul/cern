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

    using Acc = Accelerator<double>;
    // will generate 1 per second
    std::cout << "Generate phasespace lines" << std::endl;
    auto acc = Acc::getLHC();
    int freq = int(acc.f_rev);
    int turns = seconds*freq;

    auto progGen = stron::ProgramGenerator<Acc>(acc);
    auto program = progGen.create(turns, ProgramType::LHCRamp);

    for (int i = 0; i < seconds; ++i) {
        for (int j = 0; j < freq; ++j) program->step();
        std::stringstream ss;
        ss << "phasespace/" << i << "lines.dat";
        writePhasespaceFrame(acc, ss.str());
    }
}

int main(int argc, char* argv[])
{
    std::vector<std::string> args(argv, argv + argc);

    tbb::task_scheduler_init init;

    using SimpleSynchrotron = stron::SimpleSynchrotron<double>;
    using Acc = typename SimpleSynchrotron::Acc;
    using ParticleGenerator = particles::ParticleGenerator<Acc>;
    using ProgramGenerator = stron::ProgramGenerator<Acc>;

    // CHANGE FOR DIFFERENT SIMULATIONS
    stron::ProgramType progType = stron::ProgramType::LHCRamp;

    SimpleSynchrotron ss(Acc::getLHC());
    ParticleGenerator partGen(ss.getAcc());
    ProgramGenerator progGen(ss.getAcc());

    if (args.size() < 2) {
        std::cout << "Not enough arguments specified" << std::endl;
    
    } else if (args[1] == "lossmap" || args[1] == "startdist" || args[1] == "export") {
        // We often work with these together, so we make sure we have the same
        // particle distribution for both
        auto p = partGen.create(5000, particles::CCombined, particles::DoubleGaussian);
        ss.addParticles(p);
        if (args[1] == "lossmap")
            ss.simulateTurns(progGen.create(20*11245, progType), path::PATH_FILE, 11245);
        else if (args[1] == "export") {
            //ss.simulateTurns(progGen.create(100, progType));
            //for (auto& phi : p->momentum) {
                //auto r = static_cast<int>(std::floor(phi/(2*cnst::pi)));
                //phi -= 2.0*cnst::pi*r;
            //}
            particles::sixtrackExport(ss.getAcc(), *p, "collimat_export.txt");
            particles::sixtrackExport_nonCollimat(ss.getAcc(), *p, "non_collimat_export.txt");
        }

    } else if (args[1] == "inport") {
        std::cout << "Inporting particles from '" << args[2] << "'" << std::endl;
        ss.addParticles(particles::sixtrackInport(ss.getAcc(), args[2]));
        //ss.simulateTurns(progGen.create(20*11245, progType), path::PATH_FILE, 11245);

    } else if (args[1].find("animate") == 0) {
        ss.addParticles(partGen.create(1000, particles::AroundSeparatrix, particles::DoubleGaussian));
        if (args[1] == "animate") {
            ss.simulateTurns(progGen.create(1000, progType), path::PATH_FILE, 2);
        } else if (args[1] == "animate-long") {
            ss.simulateTurns(progGen.create(300*11245, progType), path::PATH_FILE, 11245);
        } else {
            throw std::runtime_error("Invalid function call");
        }
        writePhasespaceFrame(ss.getAcc(), path::LINE_FILE);

    } else if (args[1] == "sixtrack-comp") {
        std::cout << "Sixtrack comparison" << std::endl;
        double momentum = 0.0;
        if (args.size() == 3)
            momentum = std::stod(args[2])*1e6;
        std::cout << "∆E = " << std::setprecision(16) << momentum << std::endl;

        auto p = particles::ParticleCollection<double>::create(1);
        p->momentum[0] = momentum;
        p->phase[0] = cnst::pi;
        ss.addParticles(p);
        ss.simulateTurns(progGen.create(224900, progType), path::SIXTRACK_TEST_FILE);

    } else if (args[1] == "phasespace") {
        writePhasespaceFrame(SimpleSynchrotron::Acc::getLHC(), path::LINE_FILE);
    } else if (args[1] == "phasespace-mov") {
        generatePhasespaceLines(300);

    } else if (args[1] == "onep") {
        int n = 1;
        auto p = particles::ParticleCollection<double>::create(n);
        for (int i = 0; i < n; ++i) {
            p->x[i] = 0.001;
            p->px[i] = 0;
            p->phase[i] = (1 + i)/2.0*cnst::pi;
            p->momentum[i] = -1.0e8 + stron::levelCurve(ss.getAcc(), p->phase[i], stron::separatrix(ss.getAcc()));
        }
        ss.addParticles(p);
        ss.simulateTurns(progGen.create(1000, stron::ProgramType::NoRamp), path::PATH_FILE, 1);


    } else if (args[1] == "x-test") {
        int n = 15;
        auto p = particles::ParticleCollection<double>::create(n);
        for (int i = 0; i < n; ++i) {
            p->x[i] = 0;//i;
            p->px[i] = 0;
            p->phase[i] = cnst::pi;
            p->momentum[i] = 1e8;
        }
        ss.addParticles(p);
        //ss.simulateTurns(progGen.create(1000, stron::ProgramType::NoRamp), path::PATH_FILE, 1);
        ss.simulateTurns(progGen.create(350, stron::ProgramType::LHCRamp), path::PATH_FILE, 1);

    } else if (args[1] == "p-test") {
        int n = 10;
        auto p = particles::ParticleCollection<double>::create(n);
        for (int i = 0; i < n; ++i) {
            p->x[i] = p->px[i] = 0;
            p->phase[i] = cnst::pi;
            p->momentum[i] = 1e9 + 1e8*(i+1);
        }
        ss.addParticles(p);
        ss.simulateTurns(progGen.create(1000, stron::ProgramType::NoRamp), path::PATH_FILE, 1);
        
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
