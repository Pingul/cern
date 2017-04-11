#ifndef RAMP_PROGRAM
#define RAMP_PROGRAM

#include "common.hh"
#include "accelerator.hh"
#include "settings.hh"

namespace stron {
namespace ramp {

template <typename Acc>
class Program
{
public:
    using Ptr = std::shared_ptr<Program>;

    Program(Acc& acc, unsigned steps)
        : mAcc(acc), mSteps(steps) {}
    virtual ~Program() {}
    virtual void setup() {}
    virtual void step() {}
    virtual unsigned steps() { return mSteps; }
protected:
    Acc& mAcc;
    unsigned mSteps;
};


/*
 * Energy programs
 *
 */
template <typename Acc, typename Prog = Program<Acc>>
class EnergyProgram : virtual public Prog
{
public:
    EnergyProgram(Acc& acc, unsigned steps)
        : Prog(acc, steps), mIndex(0) {}

    EnergyProgram(Acc& acc, unsigned steps, std::string energyProgram)
        : Prog(acc, steps), mIndex(0)
    {
        using common::skip;
        
        mEnergy.reserve(steps);
        std::cout << "Reading '" << energyProgram << "'..." << std::endl;
        std::ifstream file(energyProgram);
        if (!file.is_open())
            throw std::runtime_error("Could not open file");
        for (unsigned i = 0; i <= steps; ++i) {
            T data;
            file >> skip >> data;
            mEnergy.emplace_back(data*1e6);
        }
    }

    // Parameter passing for constructors
    virtual void setup() override { mIndex = 0; Prog::mAcc.setE(mEnergy.at(mIndex), true); }
    virtual void step() override { Prog::mAcc.setE(mEnergy.at(++mIndex)); }
protected:
    using T = typename Acc::ValType;

    unsigned mIndex;
    std::vector<T> mEnergy;
};

template <typename Acc, typename EnergyProgram = EnergyProgram<Acc>>
class DefaultEnergyProgram : public EnergyProgram
{
    using Prog = Program<Acc>;
public:
    DefaultEnergyProgram(Acc& acc, unsigned steps)
        : Prog(acc, steps), EnergyProgram(acc, steps, LHC_RAMP_FILE) {}
};

template <typename Acc, typename EnergyProgram = EnergyProgram<Acc>>
class AggressiveEnergyProgram : public EnergyProgram
{
    using Prog = Program<Acc>;
public:
    AggressiveEnergyProgram(Acc& acc, unsigned steps)
        : Prog(acc, steps), EnergyProgram(acc, steps)
    {
        this->mEnergy.reserve(steps);
        for (unsigned i = 0; i <= steps; ++i) this->mEnergy.emplace_back(450e9 + i*1e7);
    }
};


/*
 * Collimator programs
 *
 */
template <typename Acc, typename Prog = Program<Acc>>
class CollimatorProgram : virtual public Prog
{
public:
    CollimatorProgram(Acc& acc, unsigned steps, std::string motorProgram, unsigned collIndex)
        : Prog(acc, steps), mCollIndex(collIndex) 
    {
        using common::skip;

        std::cout << "Reading '" << motorProgram << "'" << std::endl;
        std::ifstream coll_file(motorProgram);
        if (!coll_file.is_open())
            throw std::runtime_error("Could not open file");

        mData.reserve(steps);
        for (unsigned i = 0; i <= steps; ++i) {
            T left, right;
            coll_file >> skip >> left >> right;
            mData.emplace_back(std::make_pair(left, right));
        }
    }
    virtual void setup() override { mIndex = 0; set(); }
    virtual void step() override { ++mIndex; set(); }
protected:
    using T = typename Acc::ValType;

    void set()
    {
        // The data is in seconds but mIndex is in turns -- we need to convert and interpolate
        T s = mIndex/cnst::s_to_turn;
        int i_s = s;
        T frac = s - std::floor(s);

        T left, right;
        if (mIndex < Prog::mSteps) {
            // linear interpolation
            left = mData.at(i_s).first*(T(1) - frac) + mData.at(i_s + 1).first*frac;
            right = mData.at(i_s).second*(T(1) - frac) + mData.at(i_s + 1).second*frac;
        } else {
            left = mData.at(i_s).first;
            right = mData.at(i_s).second;
        }

        typename Acc::Collimat& coll = Prog::mAcc.collimators.at(mCollIndex);
        coll.left = left; 
        coll.right = right; 
    }

    unsigned mIndex;
    unsigned mCollIndex;
    std::vector<std::pair<T, T>> mData;
};

template <typename Acc, typename CollimatorProgram = CollimatorProgram<Acc>>
class TCP_IR3Program : public CollimatorProgram
{
    using Prog = Program<Acc>;
public:
    TCP_IR3Program(Acc& acc, unsigned steps)
        : Prog(acc, steps), CollimatorProgram(acc, steps, COLL_MOTOR_FILE, 0) {}
};


/*
 * Voltage reference programs
 *
 */
template <typename Acc, typename Prog = Program<Acc>>
class VoltageProgram : virtual public Prog
{
public:
    VoltageProgram(Acc& acc, unsigned steps)
        : Prog(acc, steps), mIndex(0)
    {}
    virtual void setup() override { mIndex = 0; }
    virtual void step() override { ++mIndex; Prog::mAcc.V_rf = (6 + V_k*mIndex)*1e6; }
protected:
    using T = typename Acc::ValType;

    // Calculated from the LHC ramp
    const T V_k = 2.9491187074838457087e-07;
    unsigned mIndex;
};


/*
 * Combined programs
 *
 */
template <typename Acc, typename ...Programs>
class ProgramGenerator : public Programs...
{
    using Base = Program<Acc>;
public:
    ProgramGenerator(Acc& acc, unsigned steps) : Base(acc, steps), Programs(acc, steps)... {} 
    
    virtual void setup() override { int dummy[] = {0, (Programs::setup(), void(), 0)...}; static_cast<void>(dummy); }
    virtual void step() override { int dummy[] = {0, (Programs::step(), void(), 0)...}; static_cast<void>(dummy); }
};

enum ProgramType
{
    NoRamp,
    LHCRamp,
    AggressiveRamp,
};

template <typename Acc>
typename Program<Acc>::Ptr create(Acc& acc, unsigned steps, ProgramType type)
{
    using Ptr = typename Program<Acc>::Ptr;
    switch (type)
    {
        case NoRamp:
            return Ptr(new Program<Acc>(acc, steps));
        case LHCRamp:
            return Ptr(new ProgramGenerator< 
                    Acc, 
                    DefaultEnergyProgram<Acc>, 
                    TCP_IR3Program<Acc>,
                    VoltageProgram<Acc>
                >(acc, steps));
        case AggressiveRamp:
            return Ptr(new ProgramGenerator<Acc, AggressiveEnergyProgram<Acc>>(acc, steps));
        default:
            throw std::runtime_error("ProgramType not supported");
    }
}


} // namespace ramp
} // namespace stron

#endif