#ifndef PARTICLES_HH
#define PARTICLES_HH

#include <string>
#include <vector>
#include <random>

#include "hamiltonian.hh"
#include "accelerator.hh"

//namespace stron {
namespace stron {

template <typename T>
struct ParticleCollection 
{
    using Ptr = std::shared_ptr<ParticleCollection>; // Should maybe be unique?

    ParticleCollection(int n)
        : momentum(n), phase(n), x(n), px(n) {}

    std::vector<T> momentum;
    std::vector<T> phase;
    std::vector<T> x;
    std::vector<T> px;

    //void write(std::string filePath) const;
    size_t size() const { return momentum.size(); }

    // To avoid mistakes
    ParticleCollection(const ParticleCollection&) = delete;
    ParticleCollection(ParticleCollection&&) = delete;
};

namespace pdist {

enum DIST_TYPE
{
    RANDOM,
    AROUND_SEPARATRIX,
    LINEAR,
    EXPONENTIAL,
};

template <typename T, typename Acc = Accelerator<T>>
typename ParticleCollection<T>::Ptr AroundSeparatrix(int n, const Acc& acc)
{
    using Particles = ParticleCollection<T>;
    auto particles = typename Particles::Ptr(new Particles(n));

    std::random_device rdev;
    std::mt19937 generator(rdev());

    std::uniform_real_distribution<> e_dist(-0.5e9, 0.5e9);
    std::uniform_real_distribution<> ph_dist(0, 2*cnst::pi);
    int count = 0;
    const T sep = separatrix<T>(acc);
    while (count < n) {
        const T deltaE = e_dist(generator);
        const T phase = ph_dist(generator);
        const T H = hamiltonian(acc, deltaE, phase);
        if ((sep - 1e6) < H && H < (sep + 1e6)) {
            particles->momentum[count] = deltaE;
            particles->phase[count] = phase;
            count++;
        }
    }
    return particles;
}

template <typename T, typename Acc = Accelerator<T>>
typename ParticleCollection<T>::Ptr ActionValues(int n, const Acc& acc)
{
    // n here is a request for number of particles -- the resulting might be lower

    using Particles = ParticleCollection<T>;

    const T sep = separatrix(acc);
    std::vector<T> d_actions;
    
    int outside = 80;
    double maxdE = 1.6e9;
    double de = maxdE/double(outside);
    for (int i = 1; i <= outside; ++i) {
        d_actions.push_back(hamiltonian(acc, de*double(i), cnst::pi));
    }
    
    int inside = 20;
    for (int i = 1; i <= inside; ++i) {
        d_actions.push_back(T(-8000 + 500*i));
    }

    const int n_per_level = n/(2*d_actions.size());
    const int N = 2*n_per_level*d_actions.size();
    auto particles = typename Particles::Ptr(new Particles(N));

    int count = 0;
    std::random_device rdev;
    std::mt19937 generator(rdev());
    for (T action : d_actions) {
        action += sep;
        std::uniform_real_distribution<> dist(0.0, 2*cnst::pi);
        for (T sign : std::vector<T>({-1.0, 1.0})) {
            for (int i = 0; i < n_per_level; ++i) {
                const T phase = dist(generator);
                const T energy = levelCurve(acc, phase, action, sign);
                if (std::isnan(energy)) { --i; continue; }
                particles->momentum[count] = energy;
                particles->phase[count++] = phase;
            }
        }
    }
    return particles;
}

template <typename T, typename Acc = Accelerator<T>>
typename ParticleCollection<T>::Ptr LinearDecay(int n, const Acc& acc)
{
    using Particles = ParticleCollection<T>;
    auto particles = typename Particles::Ptr(new Particles(n));

    const T sep = separatrix(acc);
    std::random_device rdev;
    std::mt19937 generator(rdev());

    // Distribution as (relative to separatrix)
    //  -15k        -10k         +5k
    //    | constant  | linear dec.|
    std::vector<T> Hs{sep - 15e3, sep - 10e3, sep + 5e3, sep + 2.7e7};
    std::vector<T> prob{100.0, 100.0, 0.1, 0.1};
    std::piecewise_linear_distribution<> H_dist(Hs.begin(), Hs.end(), prob.begin());

    std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);

    for (int i = 0; i < n; ++i) {
        const T phase = uni_dist(generator);
        const T sign = uni_dist(generator) < cnst::pi ? 1.0 : -1.0;
        const T action = H_dist(generator);
        const T energy = levelCurve(acc, phase, action, sign);
        if (std::isnan(energy)) { --i; continue; }
        particles->momentum[i] = energy;
        particles->phase[i] = phase;
    }
    return particles;
}

template <typename T, typename Acc = Accelerator<T>>
typename ParticleCollection<T>::Ptr ExponentialDecay(int n, const Acc& acc)
{
    using Particles = ParticleCollection<T>;
    auto particles = typename Particles::Ptr(new Particles(n));    

    std::random_device rdev;                                       
    std::mt19937 generator(rdev());                                
                                                                   
    std::exponential_distribution<> H_dist(0.0003);                
    std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);    
                                                                   
    const T H_low = separatrix(acc) - 15000;                      
    for (int i = 0; i < n; ++i) {                                  
        const T phase = uni_dist(generator);                       
        const T sign = uni_dist(generator) < cnst::pi ? 1.0 : -    1.0;
        const T action = H_low + H_dist(generator);                
        const T energy = levelCurve(acc, phase, action, sign);    
        if (std::isnan(energy)) { --i; continue; }                 
        particles->momentum[i] = energy;                          
        particles->phase[i] = phase;                              
    }
    return particles;
}

template <typename T>
typename ParticleCollection<T>::Ptr SixTrackTest(T momentum)
{
    using Particles = ParticleCollection<T>;
    auto p = typename Particles::Ptr(new Particles(1));
    p->momentum[0] = momentum;
    p->phase[0] = cnst::pi;
    return p;
}

//template <typename T>
//typename ParticleCollection<T>::Ptr generate(int n, DIST_TYPE) 
//{ return typename ParticleCollection<T>::Ptr(new ParticleCollection<T>(n)); }

} // namespace pdist


} // namespace stron

#endif
