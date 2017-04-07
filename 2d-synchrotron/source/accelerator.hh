#ifndef ACCELERATOR
#define ACCELERATOR

#include <fstream>
#include <iomanip>
#include <cmath>
#include "settings.hh"

namespace stron {

template <typename T>
struct Collimator
{
    enum Type {
        TCP_IR3,
        TCPc_IR7,
    };

    Collimator(enum Type t, T l, T r)
        : type(t), left(l), right(r) {}

    enum Type type;
    // in mm
    T left;
    T right;
};

template <typename T>
struct Accelerator
{
    using ValType = T;
    using Acc = Accelerator<T>;
    using Collimator = Collimator<T>;

private:
    // Use accessor methods for these instead
    T mE_ref;
    T mE_pref;
public:
    T C;
    T V_rf;
    T f_rf;
    T h_rf;
    T k_rf;
    T m_compaction;
    T coll_top;
    T coll_bot;
    T f_rev;
    T w_rev;
    std::vector<Collimator> collimators;

    static Acc getLHC() 
    {
        // Parameters for LHC can be found at http://atlas.physics.arizona.edu/~kjohns/downloads/linac/UTclass2lhc.pdf
        Acc acc;
        acc.C = 26658.8832;
        acc.mE_ref = T(450e9); // eV
        acc.mE_pref = acc.mE_ref;
        acc.V_rf = T(6e6); // V
        acc.h_rf = T(35640);
        acc.k_rf = acc.h_rf*T(2)*cnst::pi/acc.C;
        acc.m_compaction = T(0.0003225);
        acc.coll_top = T(0.5e9); // ∆eV
        acc.coll_bot = T(-0.5e9);
        
        auto p = acc.calcParticleProp(0.0, 0.0); // This is ok: we only want b
        acc.f_rev = p.b*cnst::c/acc.C;
        acc.f_rf = acc.f_rev*acc.h_rf;

        acc.w_rev = 2*cnst::pi*acc.f_rev; // Hz
        
        // Raw data from Timber measured in mm
        acc.collimators.emplace_back(Collimator::Type::TCP_IR3, -7.385, 8.285);
        acc.collimators.emplace_back(Collimator::Type::TCPc_IR7, -6.485, 5.425);

        return acc;
    }

    static Acc getLHC_NOCOLL()
    {
        Acc a = getLHC();
        a.coll_top = std::numeric_limits<T>::max();
        a.coll_bot = std::numeric_limits<T>::min();;
        return a;
    }

    struct ParticleProp 
    { 
        T pc;
        T g;    // gamma
        T g_2;    // 1/gamma^2
        T eta;
        T b2;    // beta^2
        T b;    // beta
        T W2;    // Omega2
    };
    ParticleProp calcParticleProp(T de, T phase) const 
    {
        ParticleProp p;
        //p.g = (E() + de)/cnst::m_proton;

        // We treat ∆E really as ∆pc
        const T pc_E0 = (E() + de)/cnst::m_proton;
        p.g = std::sqrt(T(1) + pc_E0*pc_E0);
        p.g_2 = T(1)/(p.g*p.g);
        p.eta = p.g_2 - m_compaction;
        p.b2 = T(1) - p.g_2;
        p.b = std::sqrt(p.b2);
        p.W2 = w_rev*w_rev*h_rf*p.eta*V_rf/(T(2)*cnst::pi*p.b2*E())*std::cos(lag_phase()); 
        return p;
    }


    void setE(T v, bool reset = false) { mE_pref = mE_ref; if (reset) mE_pref = v; mE_ref = v; }
    T E() const { return mE_ref; }
    T E_prev() const { return mE_pref; }
    T lag_phase() const { return cnst::pi - (std::asin((E() - E_prev())/V_rf)); }
    //T lag_phase() const { return cnst::pi - 0.3; }
};

} // namespace stron

#endif
