#include "sampled_distribution.hh"
#include <iostream>
#include <fstream>

int main()
{
    auto logFunc = [](double x) { 
        const double k = -1.0;
        const double m = 10;
        return x*(k*std::log(x) + m - k); // PDF(x) = k*log(x) + m
    };
    auto sinFunc = [](double x) { return x + std::cos(x); }; // PDF(x) = 1 - sin(x)
    auto f = [](double x) {
        const double k = -0.905787102751, m = 14.913170454;
        if (x >= -8000 && x < 1) return m*x;
        else if (x >= 1 && x <= 1.4e7) return x*(k*std::log(x) + m - k);
        else throw std::runtime_error("lambda: out of range");
    };

    auto pdf = [](double x){ 
        const double m = 1;
        if (x < 0 || x > m) throw std::runtime_error("CInside_E pdf - invalid value");
        const double a = -6.07;
        const double b = 0.75;
        return std::exp(a*std::pow(x/m, b));
    };

    //Integral<> I([](double x){ return x*x*x; });
    //Integral<> I([](double x){ return std::exp(-x); });
    //double a = 0;
    //double b = 1;
    //std::cout << "Integral (" << a << ", " << b << ") = " << I(a, b) << std::endl;

    std::mt19937 gen;
    //Sampled_distribution<> dist(logFunc, 1.0, 1e4);
    //Sampled_distribution<> dist(sinFunc, 0.0, 6.28);
    Sampled_distribution<> dist(pdf, 0.0, 1.0, Sampled_distribution<>::PDF);

    std::ofstream file("d.txt");
    for (int i = 0; i < 100000; i++) file << dist(gen) << std::endl;
    std::cout << "runs ok" << std::endl;
}
