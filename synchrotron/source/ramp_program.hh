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
template <typename Acc, typename Program = Program<Acc>>
class EnergyRamp : virtual public Program
{
public:
    EnergyRamp(Acc& acc, unsigned steps)
        : Program(acc, steps), mIndex(0) {}

    EnergyRamp(Acc& acc, unsigned steps, std::string energyProgram)
        : Program(acc, steps), mIndex(0)
    {
        using common::skip;
        
        mEnergy.reserve(steps);
        std::cout << "Reading '" << energyProgram << "'..." << std::endl;
        std::ifstream file(energyProgram);
        if (!file.is_open())
            throw std::runtime_error("Could not open file");
        for (int i = 0; i <= steps; ++i) {
            T data;
            file >> skip >> data;
            mEnergy.emplace_back(data*1e6);
        }
    }

    // Parameter passing for constructors
    virtual void setup() override { mIndex = 0; Program::mAcc.setE(mEnergy.at(mIndex), true); }
    virtual void step() override { Program::mAcc.setE(mEnergy.at(++mIndex)); }
protected:
    using T = typename Acc::ValType;

    unsigned mIndex;
    std::vector<T> mEnergy;
};

template <typename Acc, typename EnergyRamp = EnergyRamp<Acc>>
class DefaultEnergyRamp : public EnergyRamp
{
    using Program = Program<Acc>;
public:
    DefaultEnergyRamp(Acc& acc, unsigned steps)
        : Program(acc, steps), EnergyRamp(acc, steps, LHC_RAMP_FILE) {}
};

template <typename Acc, typename EnergyRamp = EnergyRamp<Acc>>
class AggressiveEnergyRamp : public EnergyRamp
{
    using Program = Program<Acc>;
public:
    AggressiveEnergyRamp(Acc& acc, unsigned steps)
        : Program(acc, steps), EnergyRamp(acc, steps)
    {
        this->mEnergy.reserve(steps);
        for (int i = 0; i <= steps; ++i) this->mEnergy.emplace_back(450e9 + i*1e7);
    }
};


/*
 * Collimator programs
 *
 */
template <typename Acc, typename Program = Program<Acc>>
class CollimatorRamp : virtual public Program
{
public:
    CollimatorRamp(Acc& acc, unsigned steps, std::string motorProgram, unsigned collIndex)
        : Program(acc, steps), mCollIndex(collIndex) 
    {
        using common::skip;

        std::cout << "Reading '" << motorProgram << "'" << std::endl;
        std::ifstream coll_file(motorProgram);
        if (!coll_file.is_open())
            throw std::runtime_error("Could not open file");

        mData.reserve(steps);
        for (int i = 0; i <= steps; ++i) {
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
        if (mIndex < Program::mSteps) {
            // linear interpolation
            left = mData.at(i_s).first*(T(1) - frac) + mData.at(i_s + 1).first*frac;
            right = mData.at(i_s).second*(T(1) - frac) + mData.at(i_s + 1).second*frac;
        } else {
            left = mData.at(i_s).first;
            right = mData.at(i_s).second;
        }

        typename Acc::Collimator& coll = Program::mAcc.collimators.at(mCollIndex);
        coll.left = left; 
        coll.right = right; 
    }

    unsigned mIndex;
    unsigned mCollIndex;
    std::vector<std::pair<T, T>> mData;
};

template <typename Acc, typename CollimatorRamp = CollimatorRamp<Acc>>
class TCP_IR3Ramp : public CollimatorRamp
{
    using Program = Program<Acc>;
public:
    TCP_IR3Ramp(Acc& acc, unsigned steps)
        : Program(acc, steps), CollimatorRamp(acc, steps, COLL_MOTOR_FILE, 0) {}
};



/*
 * Combined programs
 *
 */
template <typename Acc, typename ...Programs>
class RampGenerator : public Programs...
{
    using Base = Program<Acc>;
public:
    RampGenerator(Acc& acc, unsigned steps) : Base(acc, steps), Programs(acc, steps)... {} 
    
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
            return Ptr(new RampGenerator<Acc, DefaultEnergyRamp<Acc>, TCP_IR3Ramp<Acc>>(acc, steps));
        case AggressiveRamp:
            return Ptr(new RampGenerator<Acc, AggressiveEnergyRamp<Acc>>(acc, steps));
        default:
            throw std::runtime_error("ProgramType not supported");
    }
}


} // namespace ramp
} // namespace stron

#endif
