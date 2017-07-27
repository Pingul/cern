#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <tbb/task_scheduler_init.h>

#include "settings.hh"
#include "accelerator.hh"
#include "hamiltonian.hh"
#include "ramp_program.hh"

int main(int argc, char* argv[])
{
    using Acc = stron::Accelerator<double>;
    std::vector<std::string> args(argv, argv + argc);
    Acc LHC = Acc::getLHC();

    if (args.size() < 2) {
        std::cout << "no arguments" << std::endl;
    } else if (args[1] == "ham-values") {
        int n = 50;
        double maxdE = 1.9e9;
        double de = maxdE/double(n);
        for (int i = 0; i < n; ++i) {
            double d = stron::hamiltonian(LHC, de*double(i), cnst::pi);
            std::cout << d << std::endl;
        }
    } else if (args[1] == "motor") {
    } else if (args[1] == "levelCurve") {
        std::stringstream ss(args[2]);
        double h;
        ss >> h;
        std::cout << "∆H = " << h << std::endl;
        h += stron::separatrix(LHC);
        std::cout << "∆E = " << stron::levelCurve(LHC, cnst::pi, h) << std::endl;
    } else if (args[1] == "hamiltonian") {
        std::stringstream ss(args[2]);
        double e;
        ss >> e;
        std::cout << "∆E = " << e << std::endl;
        std::cout << "∆H(π, ∆E) = " << (stron::hamiltonian(LHC, e, cnst::pi) - stron::separatrix(LHC))  << std::endl;
    } else if (args[1] == "tbb") {
        // This is most likely unreliable
        tbb::task_scheduler_init init;
        std::cout << "Threads: " << init.default_num_threads() << std::endl;
    } else if (args[1] == "acc") {
        std::cout << std::setprecision(16) << stron::Accelerator<double>::getLHC().f_rf << std::endl;
    } else {
        std::cout << "unrecognised action" << std::endl;
    }

}

