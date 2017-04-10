#include <iostream>
#include <vector>
#include <string>

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
            double d = stron::hamiltonian(LHC, de*double(i), stron::cnst::pi);
            std::cout << d << std::endl;
        }
    } else if (args[1] == "motor") {
        unsigned steps = 11246;
        //auto prog = typename Program::Ptr(new stron::CollimatorRamp<Acc>(LHC, steps, "resources/motor_tcp_ir3_f5433b1.txt", 0));
        auto prog = stron::ramp::create(LHC, steps, stron::ramp::ProgramType::LHCRamp);
        for (int i = 0; i < steps - 1; ++i) {
            prog->step();
            std::cout << LHC.E() << " " << LHC.collimators[0].left << " " << LHC.collimators[0].right << std::endl;
        }
        
    }
}

