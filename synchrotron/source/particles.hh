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

static const char* LONGITUDINAL_DIST_NAMES[] = {
    "AroundSeparatrix", 
    "AVFull", 
    "AV inside, uniform H", 
    "AV inside, uniform E",
    "Cont. outside, exponential",
    "Cont. outside, non-uniform exponential dist. above/below bucket",
    "Cont. inside, linear",
    "AV outside, uniform H",
    "AV outside, unfform E",
};
enum LongitudinalDist
{
    AroundSeparatrix,
    AVFull,
    AVInside_H,
    AVInside_E,
    COutside_exp,
    COutside_ab,
    CInside_lin,
    AVOutside_H,
    AVOutside_E,
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
            case AroundSeparatrix: aroundSep(*p); break;
            case AVFull: AVFullRange(*p); break;
            case AVInside_H: AVInRange(*p, -7000, 1000, 100, /*uniform_in_H*/true); break;
            case AVInside_E: AVInRange(*p, -7000, 1000, 100, /*uniform_in_H*/false); break;
            case COutside_exp: {
                auto pdf = [](T x){ 
                    const T m = 1;
                    if (x < 0 || x > m) throw std::runtime_error("CInside_E - range");
                    const T a = -6.07;
                    const T b = 0.75;
                    return std::exp(a*std::pow(x/m, b));
                };
                Sampled_distribution<T> dist(pdf, 0, 1, Sampled_distribution<T>::PDF);

                const T sep = separatrix(mAcc);
                auto pnext = [&](Generator& g) { return 1.75e7*dist(g) + sep - 7000; };
                generateAVDist(*p, pnext);
            } break;
            case COutside_ab: {
                const T max = 3.20e7;
                auto pdf = [](T x, T a, T b, T c) {
                    if (x < 0 || x > 1) throw std::runtime_error("COutside_ab - range");
                    return std::exp(a*std::pow(x, b) + c);
                };

                //auto a_pdf = [&](T x) { return pdf(x, -8.06, 0.75, 0.28); };
                //auto b_pdf = [&](T x) { return pdf(x, -5.47, 0.36, 1.18); };
                //const T a_ratio = 0.53;

                auto a_pdf = [&](T x) { return pdf(x, -10.997, 0.785, 0.401); };
                auto b_pdf = [&](T x) { return pdf(x, -10.032, 0.703, 0.488); };
                const T a_ratio = 0.49;
                
                Sampled_distribution<T> a_dist(a_pdf, 0, 1, Sampled_distribution<T>::PDF);
                Sampled_distribution<T> b_dist(b_pdf, 0, 1, Sampled_distribution<T>::PDF);

                const T sep = separatrix(mAcc);
                std::uniform_real_distribution<> dist(0.0, 1.0);
                for (int i = 0; i < p->size(); ++i) {
                    const T phase = dist(mGenerator)*2.0*cnst::pi;
                    const T sign = dist(mGenerator) > a_ratio ? 1 : -1;
                    const T action = max*(sign > 0 ? a_dist(mGenerator) : b_dist(mGenerator)) + sep;
                    const T energy = levelCurve(mAcc, phase, action, sign);
                    if (std::isnan(energy)) { --i; continue; }
                    p->momentum[i] = energy;
                    p->phase[i] = phase;
                }
            } break;
            case CInside_lin: {
                const T sep = separatrix(mAcc);
                std::vector<T> q{0.8152908985876791, 7.14053036791e-05};
                std::vector<T> d{-7000 + sep, 0 + sep};
                std::piecewise_linear_distribution<> dist(d.begin(), d.end(), q.begin());
                auto pnext = [&](Generator& g) { return dist(g); };
                generateAVDist(*p, pnext);
            } break;
            case AVOutside_H: AVInRange(*p, 0, 1.75e7, 100, /*uniform_in_H*/true); break;
            //case AVOutside_E: AVInRange(*p, 0, 1.75e7, 15, [>uniform_in_H<]false); break;
            case AVOutside_E: AVInRange(*p, 0, 3.20e7, 35, /*uniform_in_H*/false); break;
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

    template <typename AVGen>
    void generateAVDist(PColl& particles, AVGen nextAV, int n, int sign) 
    {
        std::uniform_real_distribution<> uni_dist(0.0, 2*cnst::pi);

        int s = particles.size();
        particles.resize(s + n);
        for (int i = s; i < s + n; ++i) {
            const T phase = uni_dist(mGenerator);
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

    void AVInRange(PColl& particles, T minH, T maxH, T n, bool uniform_in_H)
    {
        if (minH > maxH)
            throw std::runtime_error("AVInRange: minH > maxH");

        const T sep = separatrix(mAcc);
        std::vector<T> actions;
        T max = uniform_in_H ? maxH : levelCurve(mAcc, cnst::pi, sep + maxH);
        T min = uniform_in_H ? minH : levelCurve(mAcc, cnst::pi, sep + minH);
        for (int i = 0; i < n; ++i) {
            const T v = min + (max - min)/(n - 1)*i;
            if (uniform_in_H)
                actions.emplace_back(sep + v);
            else 
                actions.emplace_back(hamiltonian(mAcc, v, cnst::pi));
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
    
    const Acc& mAcc;
    std::random_device mRDev;
    Generator mGenerator;
};


} // namespace stron

#endif
