#ifndef PARTICLES_HH
#define PARTICLES_HH

#include <string>
#include <vector>
#include <random>
#include <stdexcept>

#include "hamiltonian.hh"
#include "accelerator.hh"
#include "sampled_distribution.hh"

//namespace stron {
namespace stron {

template <typename T>
struct ParticleCollection 
{
    using Ptr = std::shared_ptr<ParticleCollection>; // Should maybe be unique?

    static Ptr create(int n ) { return Ptr(new ParticleCollection(n)); }

    ParticleCollection(int n)
        : momentum(n), phase(n), x(n), px(n), mActive(n, true) {}

    std::vector<T> momentum; // [eV]
    std::vector<T> phase; // [rad]
    std::vector<T> x; // [m]
    std::vector<T> px; // 

    struct HCoord { T x, px; }; // Horizontal coordinate
    HCoord xBeta(int i, T alpha, T beta, T mu_b, T mu_g) {
        // p. 165 in Wiedemann             | relativistic beta and gamma
        T sb = std::sqrt(beta*cnst::emittance/(mu_b*mu_g));
        HCoord xc{
            x[i]*sb,
            px[i]/sb - x[i]*alpha/sb
        };
        return xc;
    }

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

static const char* LONGITUDINAL_DIST_NAMES[] = {"AroundSeparatrix", "AVFull", "AVInside", "LinearDecay", "ExponentialDecay", "LogDecay", "LogLinear"};
enum LongitudinalDist
{
    AroundSeparatrix,
    AVFull,
    AVInside,
    LinearDecay,
    ExponentialDecay,
    LogDecay,
    LogLinear,
};
std::ostream& operator<<(std::ostream& os, LongitudinalDist e)
{
    os << LONGITUDINAL_DIST_NAMES[e];
    return os;
}

static const char* TRANSVERSE_DIST_NAMES[] = {"Zero", "DoubleGaussian"};
enum TransverseDist
{
    Zero,
    DoubleGaussian,
};
std::ostream& operator<<(std::ostream& os, TransverseDist e)
{
    os << TRANSVERSE_DIST_NAMES[e];
    return os;
}

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
                aroundSep(*p);
                break;
            case AVFull:
                AVFullRange(*p);
                break;
            case AVInside:
                AVInsideBucket(*p);
                break;
            case LinearDecay:
                linDecay(*p);
                break;
            case ExponentialDecay:
                expDecay(*p);
                break;
            case LogDecay:
                logDist(*p);
                break;
            case LogLinear:
                logLinDist(*p);
                break;
            default:
                throw DistributionNotFound("longitudinal");
        }
        switch (tDist) {
            case Zero:
            {
                for (size_t i = 0; i < p->size(); ++i) {
                    p->x[i] = 0;
                    p->px[i] = 0;
                }
                break;
            }
            case DoubleGaussian:
            {
                std::normal_distribution<> d(0, 1);
                for (size_t i = 0; i < p->size(); ++i) {
                    p->x[i] = d(mGenerator);
                    p->px[i] = d(mGenerator);
                }
                break;
            }
            default:
                throw DistributionNotFound("transverse");
        }
        std::cout << "Particle distribution: " << p->size() << " particles" << std::endl
                  << "\tLongitudinal: " << lDist << std::endl
                  << "\tTransverse  : " << tDist << std::endl;
        return p;
    }
    
private:
    using Generator = std::mt19937;

    template <typename AVGen>
    void generateAVDist(PColl& particles, AVGen nextAV)
    {
        std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);
    
        for (int i = 0, n = particles.size(); i < n; ++i) {
            const T phase = uni_dist(mGenerator);
            const T sign = uni_dist(mGenerator) < cnst::pi ? 1.0 : -1.0;
            const T action = nextAV(mGenerator);
            const T energy = levelCurve(mAcc, phase, action, sign);
            if (std::isnan(energy)) { --i; continue; }
            particles.momentum[i] = energy;
            particles.phase[i] = phase;
        }
    }   

    void generateAVDist(PColl& particles, std::vector<T>& actions, int n_per_level)
    {
        std::uniform_real_distribution<> dist(0.0, 2*cnst::pi);
        int count = 0;
        for (T av : actions) {
            for (T sign : std::vector<T>({-1.0, 1.0})) {
                for (int i = 0; i < n_per_level; ++i) {
                    const T phase = dist(mGenerator);
                    const T energy = levelCurve(mAcc, phase, av, sign);
                    if (std::isnan(energy)) { --i; continue; }
                    particles.momentum.at(count) = energy;
                    particles.phase.at(count++) = phase;
                }
            }
        }
    }

    void AVFullRange(PColl& particles)
    {
        // n here is a request for number of particles -- the resulting might be lower
        std::vector<T> d_actions;
        
        const T sep = separatrix(mAcc);
        int outside = 80;
        double maxdE = 1.2e9;
        double de = maxdE/double(outside);
        for (int i = 1; i <= outside; ++i) {
            d_actions.push_back(sep + hamiltonian(mAcc, de*double(i), cnst::pi));
        }
        
        int inside = 20;
        for (int i = 1; i <= inside; ++i) {
            d_actions.push_back(sep + T(-8000 + 500*i));
        }

        const int n_per_level = particles.size()/(2*d_actions.size());
        const int N = 2*n_per_level*d_actions.size();
        particles.resize(N);
        generateAVDist(particles, d_actions, n_per_level);
    }

    void AVInsideBucket(PColl& particles)
    {
        std::vector<T> actions;
        const T sep = separatrix(mAcc);
        int n = 20;
        for (int i = 0; i < n; ++i) {
            actions.emplace_back(sep + 500*i - 8000);
        }
        const int n_per_level = particles.size()/(2*actions.size());
        const int N = 2*n_per_level*actions.size();
        particles.resize(N);
        generateAVDist(particles, actions, n_per_level);
    }

    void aroundSep(PColl& particles)
    {
        const T sep = separatrix(mAcc);
        std::uniform_real_distribution<> hdist(sep - 1e6, sep + 1e6);
        auto avGen = [&](Generator& g) { return hdist(g); };
        generateAVDist(particles, avGen);
    }


    void linDecay(PColl& particles)
    {
        const T sep = separatrix(mAcc);
        // Distribution as (relative to separatrix)
        std::vector<T> Hs{sep - 8e3, sep, sep + 1.8e7};
        std::vector<T> prob{1000.0, 1.0, 0.01};
        std::piecewise_linear_distribution<> H_dist(Hs.begin(), Hs.end(), prob.begin());

        auto avGen = [&](Generator& g) { return H_dist(g); };
        generateAVDist(particles, avGen);
    }

    void expDecay(PColl& particles)
    {
        std::exponential_distribution<> H_dist(0.0003);                
        std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);    
                                                                       
        const T H_low = separatrix(mAcc) - 15000;                      
        auto avGen = [&](Generator& g) { return H_low + H_dist(mGenerator); };
        generateAVDist(particles, avGen);
    }

    
    void logDist(PColl& particles)
    {
        auto cdf = [](T x) {
            const T k = -0.905787102751;
            const T m = 14.913170454;
            return x*(k*std::log(x) + (m - k));
        };
        Sampled_distribution<T> logDist(cdf, 1.0, 1.4e7);
        const T sep = separatrix(mAcc);
        auto avGen = [&](Generator& g) { return sep + logDist(g); };
        generateAVDist(particles, avGen);
    }

    void logLinDist(PColl& particles)
    {
        auto cdf = [](double x) {
            const double k = -0.905787102751, m = 14.913170454;
            if (x >= -8000 && x < 1) return m*x;
            else if (x >= 1 && x <= 1.4e7) return x*(k*std::log(x) + m - k);
            else throw std::runtime_error("logLinDist: out of range");
        };
        Sampled_distribution<T> dist(cdf, -8000, 1.4e7, /*resolution*/20000);
        const T sep = separatrix(mAcc);
        auto avGen = [&](Generator& g) { return sep + dist(g); };
        generateAVDist(particles, avGen);
    }
    
    
    const Acc& mAcc;
    std::random_device mRDev;
    Generator mGenerator;
};


} // namespace stron

#endif
