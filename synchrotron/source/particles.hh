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

    static Ptr create(int n ) { return Ptr(new ParticleCollection(n)); }

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

template <typename Acc>
struct ParticleGenerator
{
    using T = typename Acc::ValType;
    using PColl = ParticleCollection<T>;
    using PCollPtr = typename PColl::Ptr;

    ParticleGenerator(const Acc& acc)
        : mAcc(acc), mRDev(), mGenerator(mRDev()) {}
    
    PCollPtr AroundSeparatrix(int n)
    {
        auto particles = PColl::create(n);

        std::uniform_real_distribution<> e_dist(-0.5e9, 0.5e9);
        std::uniform_real_distribution<> ph_dist(0, 2*cnst::pi);
        int count = 0;
        const T sep = separatrix<T>(mAcc);
        while (count < n) {
            const T deltaE = e_dist(mGenerator);
            const T phase = ph_dist(mGenerator);
            const T H = hamiltonian(mAcc, deltaE, phase);
            if ((sep - 1e6) < H && H < (sep + 1e6)) {
                particles->momentum[count] = deltaE;
                particles->phase[count] = phase;
                count++;
            }
        }
        return particles;
    }

    PCollPtr ActionValues(int n)
    {
        // n here is a request for number of particles -- the resulting might be lower
        const T sep = separatrix(mAcc);
        std::vector<T> d_actions;
        
        int outside = 80;
        double maxdE = 1.6e9;
        double de = maxdE/double(outside);
        for (int i = 1; i <= outside; ++i) {
            d_actions.push_back(hamiltonian(mAcc, de*double(i), cnst::pi));
        }
        
        int inside = 20;
        for (int i = 1; i <= inside; ++i) {
            d_actions.push_back(T(-8000 + 500*i));
        }

        const int n_per_level = n/(2*d_actions.size());
        const int N = 2*n_per_level*d_actions.size();
        auto particles = PColl::create(N);

        int count = 0;
        for (T action : d_actions) {
            action += sep;
            std::uniform_real_distribution<> dist(0.0, 2*cnst::pi);
            for (T sign : std::vector<T>({-1.0, 1.0})) {
                for (int i = 0; i < n_per_level; ++i) {
                    const T phase = dist(mGenerator);
                    const T energy = levelCurve(mAcc, phase, action, sign);
                    if (std::isnan(energy)) { --i; continue; }
                    particles->momentum[count] = energy;
                    particles->phase[count++] = phase;
                }
            }
        }
        return particles;
    }

    PCollPtr LinearDecay(int n)
    {
        auto particles = PColl::create(n);    
        const T sep = separatrix(mAcc);
    
        // Distribution as (relative to separatrix)
        //  -15k        -10k         +5k
        //    | constant  | linear dec.|
        std::vector<T> Hs{sep - 15e3, sep - 10e3, sep + 5e3, sep + 2.7e7};
        std::vector<T> prob{100.0, 100.0, 0.1, 0.1};
        std::piecewise_linear_distribution<> H_dist(Hs.begin(), Hs.end(), prob.begin());
    
        std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);
    
        for (int i = 0; i < n; ++i) {
            const T phase = uni_dist(mGenerator);
            const T sign = uni_dist(mGenerator) < cnst::pi ? 1.0 : -1.0;
            const T action = H_dist(mGenerator);
            const T energy = levelCurve(mAcc, phase, action, sign);
            if (std::isnan(energy)) { --i; continue; }
            particles->momentum[i] = energy;
            particles->phase[i] = phase;
        }
        return particles;
    }

    PCollPtr ExponentialDecay(int n)
    {
        auto particles = PColl::create(n);
    
        std::exponential_distribution<> H_dist(0.0003);                
        std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);    
                                                                       
        const T H_low = separatrix(mAcc) - 15000;                      
        for (int i = 0; i < n; ++i) {                                  
            const T phase = uni_dist(mGenerator);                       
            const T sign = uni_dist(mGenerator) < cnst::pi ? 1.0 : -    1.0;
            const T action = H_low + H_dist(mGenerator);                
            const T energy = levelCurve(mAcc, phase, action, sign);    
            if (std::isnan(energy)) { --i; continue; }                 
            particles->momentum[i] = energy;                          
            particles->phase[i] = phase;                              
        }
        return particles;
    }

private:
    const Acc& mAcc;
    std::random_device mRDev;
    std::mt19937 mGenerator;
};

//enum DIST_TYPE
//{
    //RANDOM,
    //AROUND_SEPARATRIX,
    //LINEAR,
    //EXPONENTIAL,
//};



} // namespace stron

#endif
