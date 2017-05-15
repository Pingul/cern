#ifndef HAMILTONIAN
#define HAMILTONIAN


#include "settings.hh"
#include "accelerator.hh"

namespace stron {

template <typename T>
inline T hamiltonian(const Accelerator<T>& acc, T de, T ph)
{
    using std::cos;
    using std::sin;
    using cnst::pi;

    auto p = acc.calcParticleProp(de, ph);
    const T ph_s = acc.lag_phase();
    const T A = -T(0.5)*acc.h_rf*p.eta/(p.b2*acc.E());
    const T H = A*de*de + acc.V_rf/(T(2)*pi)*(cos(ph) - cos(ph_s) + (ph - ph_s)*sin(ph_s));
    return H;
}

template <typename T>
inline T levelCurve(const Accelerator<T>& acc, T ph, T H, T sign=1.0) 
{
    using std::cos;
    using std::sin;
    using std::sqrt;
    using cnst::pi;

    sign = sign/std::abs(sign);

    const T ph_s = acc.lag_phase();
    const T threshold = 0.005;
    T de = 0.0;
    do {
        const auto p = acc.calcParticleProp(de, ph);
        const T A = -T(0.5)*acc.h_rf*p.eta/(p.b2*acc.E());
        de = sign*sqrt(T(1.0)/A*(H - acc.V_rf/(T(2)*pi)*(cos(ph) - cos(ph_s) + (ph - ph_s)*sin(ph_s))));
    } while (std::abs(hamiltonian(acc, de, ph) - H) > threshold);
    return de;
}

template <typename T>
inline T separatrix(const Accelerator<T>& acc)
{
    return hamiltonian<T>(acc, 0.0, cnst::pi - acc.lag_phase());
}

template <typename T>
inline void writePhasespaceFrame(const Accelerator<T>& acc, std::string filePath)
{
    std::cout << "Writing frame to '" << filePath << "'" << std::endl;
    std::ofstream file(filePath.c_str());
    if (!file.is_open())
        throw std::runtime_error("could not open file");

    const T ph_step = 0.005;
    const int ph_steps = std::ceil((cnst::FRAME_X_HIGH - cnst::FRAME_X_LOW)/ph_step);

    const T Hstart = hamiltonian<T>(acc, 0.0, cnst::pi) + 5e4;
    const T Hstep = 2.0e5;
    const T Hdstep = 1.0e5;
    const int Hsteps = 20;
    //const int Hsteps = 0;

    int lines = 2 + 2*Hsteps;

    const T Hseparatrix = separatrix<T>(acc);

    file << lines << "," << ph_steps << std::endl;
    for (T ph = cnst::FRAME_X_LOW; ph <= cnst::FRAME_X_HIGH; ph += ph_step) {
        T de = levelCurve(acc, ph, Hseparatrix);
        file << std::setprecision(16) << de  << "," << ph << "," << Hseparatrix << std::endl
                                      << -de  << "," << ph << "," << Hseparatrix << std::endl;
        
        int n = 0;
        T H = Hstart;
        while (n++ < Hsteps) {
            de = levelCurve(acc, ph, H);
            file << std::setprecision(16) << de  << "," << ph << "," << H << std::endl
                                          << -de  << "," << ph << "," << H << std::endl;
            H += Hstep + T(n)*Hdstep;
        }
    }
}

} // namespace stron


#endif
