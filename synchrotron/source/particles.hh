#ifndef PARTICLES_HH
#define PARTICLES_HH

#include <string>
#include <vector>
#include <random>
#include <stdexcept>

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
        : momentum(n), phase(n), x(n), px(n), mActive(n, true) {}

    std::vector<T> momentum;
    std::vector<T> phase;
    std::vector<T> x;
    std::vector<T> px;

    //void write(std::string filePath) const;
    size_t size() const { return momentum.size(); }
    void resize(int n) { momentum.resize(n); phase.resize(n); x.resize(n); px.resize(n); }

    // To avoid mistakes
    ParticleCollection(const ParticleCollection&) = delete;
    ParticleCollection(ParticleCollection&&) = delete;

    bool isActive(int i) const { return mActive[i]; }
    void setActive(int i, bool v) { mActive[i] = v; }
private:
    std::vector<int> mActive; // using int to allow for concurrent writes
};

enum LongitudinalDist
{
    AroundSeparatrix,
    ActionValues,
    LinearDecay,
    ExponentialDecay,
};

enum TransverseDist
{
    Zero,
    //Random,
};

struct DistributionNotFound : public std::runtime_error
{
    DistributionNotFound(const std::string& type)
        : std::runtime_error("particles.hh: " + type) {}
};

template <typename Acc>
struct ParticleGenerator
{
    using T = typename Acc::ValType;
    using PColl = ParticleCollection<T>;
    using PCollPtr = typename PColl::Ptr;

    ParticleGenerator(const Acc& acc)
        : mAcc(acc), mRDev(), mGenerator(mRDev()) {}

    PCollPtr create(int n, LongitudinalDist lDist, TransverseDist tDist = Zero)
    {
        auto p = PColl::create(n);
        switch (lDist) {
            case AroundSeparatrix:
                GenAroundSeparatrix(*p);
                break;
            case ActionValues:
                GenActionValues(*p);
                break;
            case LinearDecay:
                GenLinearDecay(*p);
                break;
            case ExponentialDecay:
                GenExponentialDecay(*p);
                break;
            default:
                throw DistributionNotFound("longitudinal");
        }
        switch (tDist) {
            case Zero:
            {
                for (int i = 0; i < p->size(); ++i) {
                    p->x[i] = 0;
                    p->px[i] = 0;
                }
                break;
            }
            default:
                throw DistributionNotFound("transverse");
        }
        return p;
    }
    
private:
    void GenAroundSeparatrix(PColl& particles)
    {
        std::uniform_real_distribution<> e_dist(-0.5e9, 0.5e9);
        std::uniform_real_distribution<> ph_dist(0, 2*cnst::pi);
        int count = 0;
        const T sep = separatrix<T>(mAcc);
        while (count < particles.size()) {
            const T deltaE = e_dist(mGenerator);
            const T phase = ph_dist(mGenerator);
            const T H = hamiltonian(mAcc, deltaE, phase);
            if ((sep - 1e6) < H && H < (sep + 1e6)) {
                particles.momentum[count] = deltaE;
                particles.phase[count] = phase;
                count++;
            }
        }
    }

    void GenActionValues(PColl& particles)
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

        const int n_per_level = particles.size()/(2*d_actions.size());
        const int N = 2*n_per_level*d_actions.size();
        particles.resize(N);

        int count = 0;
        for (T action : d_actions) {
            action += sep;
            std::uniform_real_distribution<> dist(0.0, 2*cnst::pi);
            for (T sign : std::vector<T>({-1.0, 1.0})) {
                for (int i = 0; i < n_per_level; ++i) {
                    const T phase = dist(mGenerator);
                    const T energy = levelCurve(mAcc, phase, action, sign);
                    if (std::isnan(energy)) { --i; continue; }
                    particles.momentum[count] = energy;
                    particles.phase[count++] = phase;
                }
            }
        }
    }

    void GenLinearDecay(PColl& particles)
    {
        const T sep = separatrix(mAcc);
        // Distribution as (relative to separatrix)
        //  -15k        -10k         +5k
        //    | constant  | linear dec.|
        std::vector<T> Hs{sep - 15e3, sep - 10e3, sep + 5e3, sep + 2.7e7};
        std::vector<T> prob{100.0, 100.0, 0.1, 0.1};
        std::piecewise_linear_distribution<> H_dist(Hs.begin(), Hs.end(), prob.begin());
    
        std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);
    
        for (int i = 0, n = particles.size(); i < n; ++i) {
            const T phase = uni_dist(mGenerator);
            const T sign = uni_dist(mGenerator) < cnst::pi ? 1.0 : -1.0;
            const T action = H_dist(mGenerator);
            const T energy = levelCurve(mAcc, phase, action, sign);
            if (std::isnan(energy)) { --i; continue; }
            particles.momentum[i] = energy;
            particles.phase[i] = phase;
        }
    }

    void GenExponentialDecay(PColl& particles)
    {
        std::exponential_distribution<> H_dist(0.0003);                
        std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);    
                                                                       
        const T H_low = separatrix(mAcc) - 15000;                      
        for (int i = 0, n = particles.size(); i < n; ++i) {                                  
            const T phase = uni_dist(mGenerator);                       
            const T sign = uni_dist(mGenerator) < cnst::pi ? 1.0 : -    1.0;
            const T action = H_low + H_dist(mGenerator);                
            const T energy = levelCurve(mAcc, phase, action, sign);    
            if (std::isnan(energy)) { --i; continue; }                 
            particles.momentum[i] = energy;                          
            particles.phase[i] = phase;                              
        }
    }

private:
    const Acc& mAcc;
    std::random_device mRDev;
    std::mt19937 mGenerator;
};


} // namespace stron

#endif
