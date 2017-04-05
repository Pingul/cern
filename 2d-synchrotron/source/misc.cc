#include <iostream>
#include <vector>
#include <string>

#include "settings.hh"
#include "accelerator.hh"
#include "hamiltonian.hh"

int main(int argc, char* argv[])
{
    //using namespace twodsynch;
    std::vector<std::string> args(argv, argv + argc);
    auto acc = twodsynch::Accelerator<double>::getLHC();

    if (args.size() < 2) {
        std::cout << "no arguments" << std::endl;
    } else if (args[1] == "ham-values") {
        int n = 50;
        double maxdE = 1.9e9;
        double de = maxdE/double(n);
        for (int i = 0; i < n; ++i) {
            double d = twodsynch::hamiltonian(acc, de*double(i), twodsynch::cnst::pi);
            std::cout << d << std::endl;
        }
    }
}

