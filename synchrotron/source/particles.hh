#ifndef PARTICLES_HH
#define PARTICLES_HH

#include <string>
#include <vector>
#include <random>
#include <stdexcept>
#include <fstream>
#include <iomanip>

#include "hamiltonian.hh"
#include "accelerator.hh"
#include "sampled_distribution.hh"

namespace particles {

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
    HCoord xBeta(int i, T alpha, T beta, T mu_b, T mu_g) const {
        // p. 165 in Wiedemann             | relativistic beta and gamma
        //const T eg = cnst::emittance/(mu_b*mu_g);
        //const T gt = (1 + alpha*alpha)/beta;
        //const T sigmap = std::sqrt(gt*eg/beta);
        const T sb = std::sqrt(beta*cnst::emittance/(mu_b*mu_g));
        HCoord xc{
            x[i]*sb,
            px[i]/sb - x[i]*alpha/sb, // WARNING: This is probably not correct
            //-alpha*sigmap*x[i] + sigmap*px[i],
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

template <typename T>
void sixtrackExport(const stron::Accelerator<T>& acc, const ParticleCollection<T>& p, const std::string& file)
{
    std::ofstream f(file.c_str());
    for (size_t i = 0; i < p.size(); ++i) {
        auto prop = acc.calcParticleProp(p.momentum[i]);
        const T z = (- p.phase[i])*acc.C/(2*cnst::pi*acc.h_rf*prop.b)*1e3;
        const T e = (p.momentum[i] + acc.E())*1e-6;
        
        // hard coded for ip1
        const T x = p.x[i]*std::sqrt(cnst::emittance/prop.g*11);
        const T px = p.px[i]*std::sqrt(cnst::emittance/prop.g/11);
        const T y = 0;
        const T py = 0;
        f << std::fixed << std::setprecision(16) << x << " " << std::fixed << std::setprecision(16) << px << " " // x  px
          << std::fixed << std::setprecision(16) << y << " " << std::fixed << std::setprecision(16) << py << " " // y  py
          << std::fixed << std::setprecision(16) << z << " " << std::fixed << std::setprecision(16) << e << std::endl;
    }
    std::cout << "Exported data to '" << file << "'" << std::endl;
}

template <typename T>
void sixtrackExport_nonCollimat(const stron::Accelerator<T>& acc, const ParticleCollection<T>& p, const std::string& file)
{
    std::ofstream f(file.c_str());
    for (size_t i = 0; i < p.size(); i += 2) {
        for (size_t j = i; j < i + 2; ++j) {
            // closed orbit correction from SixTrack twiss function
            const T xcorr = -1.999999996185254147462700000000000;
            const T pxcorr = 0.000000000150460499225583966566850;
            const T ycorr = -0.000000006411949440904753153327300;
            const T pycorr = -0.169999999800309692377100000000000;

            auto prop = acc.calcParticleProp(p.momentum[j]);
            const T z = (- p.phase[j])*acc.C/(2*cnst::pi*acc.h_rf*prop.b)*1e3;
            const T dp = (p.momentum[j] - acc.E())/acc.E();

            // currently hard coded for ip1 
            const T x = p.x[j]*std::sqrt(cnst::emittance/prop.g*11)*1e3 + xcorr;
            const T px = p.px[j]*std::sqrt(cnst::emittance/prop.g/11)*1e3 + pxcorr;
            f << std::fixed << std::setprecision(16) << x << std::endl
              << std::fixed << std::setprecision(16) << px << std::endl
              << std::fixed << std::setprecision(16) << ycorr << std::endl // y
              << std::fixed << std::setprecision(16) << pycorr << std::endl // py
              << std::fixed << std::setprecision(16) << z << std::endl // path length
              << std::fixed << std::setprecision(16) << dp << std::endl;
        }
        f << std::fixed << std::setprecision(16) << acc.E()*1e-6 << std::endl;
        for (size_t j = i; j < i + 2; ++j) {
            f << std::fixed << std::setprecision(16) << (p.momentum[j] + acc.E())*1e-6 << std::endl;
        }
    }
    std::cout << "Exported data to '" << file << "'" << std::endl;
}

template <typename T>
typename ParticleCollection<T>::Ptr sixtrackInport(const stron::Accelerator<T>& acc, const std::string& file) 
{
    // count lines first
    std::string line;
    std::ifstream f(file.c_str());
    int count = 0;
    while (std::getline(f, line)) ++count;
    std::cout << "Lines: " << count << std::endl;

    f.clear();
    f.seekg(0, std::ios::beg);

    auto p = ParticleCollection<T>::create(count);
    for (int i = 0; i < count; ++i) {
        T x, px, y, py, z, p_tot;
        f >> x >> px >> y >> py >> z >> p_tot;
        p->momentum[i] = p_tot*1e6 - acc.E();
        auto prop = acc.calcParticleProp(p->momentum[i]);
        p->phase[i] = cnst::pi - z*(2*cnst::pi*acc.h_rf*prop.b)/acc.C;
        p->x[i] = 0; 
        p->px[i] = 0;
        //p->x[i] = x/std::sqrt(cnst::emittance/479.0*11);
        //p->px[i] = x/std::sqrt(cnst::emittance/479.0/11);
    }
    return p;
}

template <typename T>
typename ParticleCollection<T>::Ptr sixtrackInport_nonCollimat(const stron::Accelerator<T>& acc, const std::string& file, int howMany = std::numeric_limits<int>::max()) 
{
    std::string line;
    std::ifstream f(file.c_str());
    int count = 0;
    while (std::getline(f, line)) ++count;
    int n = 2*count/15;
    std::cout << "Lines: " << count << " --> " << n << " particles" << std::endl;
    n = howMany < n ? howMany : n;
    std::cout << "Imported " << n << " particles" << std::endl;
    
    f.clear();
    f.seekg(0, std::ios::beg);

    auto p = ParticleCollection<T>::create(n);
    for (int i = 0; i < n; i += 2) {
        T x[2], px[2], y, py, z[2], e[2], dp, e_ref;
        f >> x[0] >> px[0] >> y >> py >> z[0] >> dp 
          >> x[1] >> px[1] >> y >> py >> z[1] >> dp 
          >> e_ref >> e[0] >> e[1];
        
        // closed orbit correction from SixTrack twiss function
        const T xcorr = -1.999999996185254147462700000000000;
        const T pxcorr = 0.000000000150460499225583966566850;
        for (int j = 0; j < 2; ++j) {
            p->momentum[i + j] = (e[j] - e_ref)*1e6;
            auto prop = acc.calcParticleProp(p->momentum[i + j]);
            p->x[i + j] = (x[j] - xcorr)/std::sqrt(cnst::emittance/prop.g*11)*1e-3;
            p->px[i + j] = (px[j] - pxcorr)/std::sqrt(cnst::emittance/prop.g/11)*1e-3;
            p->phase[i + j] = -z[j]*1e-3*(2*cnst::pi*acc.h_rf*prop.b)/acc.C;
        }
    }
    return p;
}


static const char* LONGITUDINAL_DIST_NAMES[] = {
    "AroundSeparatrix", 
    "AVFull", 
    "AV inside, uniform H", 
    "AV inside, uniform E",
    "Cont. outside, exponential",
    "Cont. outside, non-uniform exponential dist. above/below bucket",
    "Cont. inside, linear",
    "Cont. inside, linearly increasing (towards separatrix)",
    "Cont. inside, linearly decreasing (towards separatrix)",
    "Cont. inside, constant",
    "Cont. combined: non-uniform exponential outside, linear inside",
    "AV outside, uniform H",
    "AV outside, uniform E",
    "OutsideColl",
    "Zero",
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
    CInside_LI,
    CInside_LD,
    CInside_C,
    CCombined,
    AVOutside_H,
    AVOutside_E,
    OutsideColl,
    LZero,
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

                //auto a_pdf = [&](T x) { return pdf(x, -10.997, 0.785, 0.401); };
                auto a_pdf = [&](T x) { return pdf(x, -20.984, 1.000, 0.320); };
                //auto b_pdf = [&](T x) { return pdf(x, -10.032, 0.703, 0.488); };
                auto b_pdf = [&](T x) { return pdf(x, -10.032, 0.703, 0.488); };
                const T a_ratio = 0.40;
                
                Sampled_distribution<T> a_dist(a_pdf, 0, 1, Sampled_distribution<T>::PDF);
                Sampled_distribution<T> b_dist(b_pdf, 0, 1, Sampled_distribution<T>::PDF);

                const T sep = separatrix(mAcc);
                std::uniform_real_distribution<> dist(0.0, 1.0);
                for (int i = 0; i < p->size(); ++i) {
                    const T phase = dist(mGenerator)*2.0*cnst::pi;
                    const T sign = dist(mGenerator) < a_ratio ? 1 : -1;
                    const T action = max*(sign > 0 ? a_dist(mGenerator) : b_dist(mGenerator)) + sep;
                    const T energy = levelCurve(mAcc, phase, action, sign);
                    if (std::isnan(energy)) { --i; continue; }
                    p->momentum[i] = energy;
                    p->phase[i] = phase;
                }
            } break;
            case CInside_LD:
            case CInside_lin: {
                const T sep = separatrix(mAcc);
                std::vector<T> q{0.8152908985876791, 7.14053036791e-05};
                std::vector<T> d{-7000 + sep, sep};
                std::piecewise_linear_distribution<> dist(d.begin(), d.end(), q.begin());
                auto pnext = [&](Generator& g) { return dist(g); };
                generateAVDist(*p, pnext);
            } break;
            case CInside_LI: {
                const T sep = separatrix(mAcc);
                std::vector<T> q{7.14053036791e-05, 0.8152908985876791};
                std::vector<T> d{-7000 + sep, sep};
                std::piecewise_linear_distribution<> dist(d.begin(), d.end(), q.begin());
                auto pnext = [&](Generator& g) { return dist(g); };
                generateAVDist(*p, pnext);
            } break;
            case CInside_C: {
                const T sep = separatrix(mAcc);
                std::vector<T> q{1.0, 1.0};
                std::vector<T> d{-7000 + sep, sep};
                std::piecewise_linear_distribution<> dist(d.begin(), d.end(), q.begin());
                auto pnext = [&](Generator& g) { return dist(g); };
                generateAVDist(*p, pnext);
            } break;
            case CCombined: {
                const T max = 3.20e7;
                auto pdf = [](T x, T a, T b, T c) {
                    if (x < 0 || x > 1) throw std::runtime_error("COutside_ab - range");
                    return std::exp(a*std::pow(x, b) + c);
                };

                //auto a_pdf = [&](T x) { return pdf(x, -10.997, 0.785, 0.401); };
                auto a_pdf = [&](T x) { return pdf(x, -20.984, 1.000, 0.320); };
                //auto b_pdf = [&](T x) { return pdf(x, -10.032, 0.703, 0.488); };
                auto b_pdf = [&](T x) { return pdf(x, -13.032, 0.703, 0.488); };
                const T a_ratio = 0.40;
                
                Sampled_distribution<T> oa_dist(a_pdf, 0, 1, Sampled_distribution<T>::PDF);
                Sampled_distribution<T> ob_dist(b_pdf, 0, 1, Sampled_distribution<T>::PDF);

                const T i_ratio = 0.20;
                std::vector<T> q{0.8152908985876791, 7.14053036791e-05};
                std::vector<T> d{-7000, 0};
                std::piecewise_linear_distribution<> i_dist(d.begin(), d.end(), q.begin());

                int pin = 0;
                int pout = 0;
                const T sep = separatrix(mAcc);
                std::uniform_real_distribution<> dist(0.0, 1.0);
                for (int i = 0; i < p->size(); ++i) {
                    const T phase = dist(mGenerator)*2.0*cnst::pi;
                    T action;
                    T sign;
                    if (dist(mGenerator) < i_ratio) {
                        // generate particle inside separatrix
                        sign = dist(mGenerator) > 0.5 ? 1 : -1;
                        action = i_dist(mGenerator);
                        ++pin;
                    } else {
                        // outside
                        sign = dist(mGenerator) < a_ratio ? 1 : -1;
                        action = max*(sign > 0 ? oa_dist(mGenerator) : ob_dist(mGenerator));
                        ++pout;
                    }   
                    action += sep;
                    const T energy = levelCurve(mAcc, phase, action, sign);
                    if (std::isnan(energy)) { --i; continue; }
                    p->momentum[i] = energy;
                    p->phase[i] = phase;
                }
                std::cout << "Particles" << std::endl;
                std::cout << "\tinside: " << pin << ", " << T(pin)/(pin+pout) << std::endl;
                std::cout << "\toutside: " << pout << ", " << T(pout)/(pin+pout) << std::endl;
            } break;
            case AVOutside_H: AVInRange(*p, 0, 1.75e7, 100, /*uniform_in_H*/true); break;
            case AVOutside_E: AVInRange(*p, 0, 3.20e7, 35, /*uniform_in_H*/false); break;
            case OutsideColl: AVInRange(*p, 3.3e7, 4.0e7, 20, false); break;
            case LZero: {
                for (size_t i = 0; i < p->size(); ++i) {
                    p->momentum[i] = 0;
                    p->phase[i] = cnst::pi;
                }
                break;
            }
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
                // Orthogonal
                //std::normal_distribution<> d(0, 1);
                //for (size_t i = 0; i < p->size(); ++i) {
                    //p->x[i] = d(mGenerator);
                    //p->px[i] = d(mGenerator);
                //}

                // Bi-gaussian
                std::normal_distribution<> rd(0, 1);
                std::uniform_real_distribution<> phid(0, 2*cnst::pi);
                for (size_t i = 0; i < p->size(); ++i) {
                    const T r = rd(mGenerator);
                    const T phi = phid(mGenerator);
                    p->x[i] = r*std::cos(phi);
                    p->px[i] = r*std::sin(phi);
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


} // namespace particles

#endif
