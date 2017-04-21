#ifndef RAMP_PROGRAM
#define RAMP_PROGRAM

#include <stdexcept>
#include "common.hh"
#include "accelerator.hh"
#include "settings.hh"

namespace stron {

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
    virtual unsigned steps() const { return mSteps; }

    Program() = delete;
    Program(const Program&) = delete;
protected:
    Acc& mAcc;
    unsigned mSteps;
};

namespace program {

// Custom exceptions for this module
namespace except {
struct ProgramOutOfBounds : public std::runtime_error 
{ ProgramOutOfBounds(const std::string& program) : std::runtime_error("ramp_program.hh: " + program ) {} };
struct FileNotFound : public std::runtime_error 
{ FileNotFound(const std::string& filename) : std::runtime_error("ramp_program.hh: " + filename) {} };
struct ProgramTypeNotFound : public std::runtime_error 
{ ProgramTypeNotFound() : std::runtime_error("ramp_program.hh: ProgramType not found") {} };

} // namespace except

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
        
        std::cout << "Reading '" << energyProgram << "'..." << std::endl;
        std::ifstream file(energyProgram);
        if (!file.is_open())
            throw except::FileNotFound(energyProgram);
        
        mEnergy.reserve(steps);
        T data;
        while (file >> skip >> data && !valid()) {
            mEnergy.emplace_back(data*1e6);
        }

        if (!valid()) 
            throw except::ProgramOutOfBounds("energy");
    }

    virtual void setup() override { mIndex = 0; Prog::mAcc.setE(mEnergy.at(mIndex), true); }
    virtual void step() override { Prog::mAcc.setE(mEnergy.at(++mIndex)); }
    virtual bool valid() const { return mEnergy.size() >= Prog::mSteps; } // 1 extra step for the setup
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
            throw except::FileNotFound(motorProgram);

        mData.reserve(steps);
        T left, right;
        while (coll_file >> skip >> left >> right && !valid()) {
            mData.emplace_back(std::make_pair(left, right));
        }

        if (!valid()) 
            throw except::ProgramOutOfBounds("collimator");
    }
    virtual void setup() override { mIndex = 0; set(); }
    virtual void step() override { ++mIndex; set(); }
    virtual bool valid() const { return mData.size() > std::ceil(Prog::mSteps/cnst::s_to_turn); } 
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
class ProgramCombiner : public Programs...
{
    using Base = Program<Acc>;
public:
    ProgramCombiner(Acc& acc, unsigned steps) : Base(acc, steps), Programs(acc, steps)... {} 
    
    // Explanation at http://stackoverflow.com/questions/43322854/multiple-inheritance-with-variadic-templates-how-to-call-function-for-each-base/43322961?noredirect=1#comment73711445_43322961
    virtual void setup() override { int dummy[] = {0, (Programs::setup(), void(), 0)...}; static_cast<void>(dummy); }
    virtual void step() override { int dummy[] = {0, (Programs::step(), void(), 0)...}; static_cast<void>(dummy); }
};


} // namespace program

enum ProgramType
{
    NoRamp,
    LHCRamp,
    LHCWithoutEnergyRamp,
    AggressiveRamp,
    Voltage,
    Energy,
};

template <typename Acc>
struct ProgramGenerator
{

    using Program = Program<Acc>;
    using ProgramPtr = typename Program::Ptr;

    ProgramGenerator(Acc& acc)
        : mAcc(acc) {}

    ProgramPtr create(unsigned steps, ProgramType type) 
    {
        using namespace program;
        switch (type)
        {
            case NoRamp:
                return ProgramPtr(new Program(mAcc, steps));
            case LHCRamp:
                return ProgramPtr(new ProgramCombiner< 
                        Acc, 
                        DefaultEnergyProgram<Acc>, 
                        TCP_IR3Program<Acc>,
                        VoltageProgram<Acc>
                    >(mAcc, steps));
            case LHCWithoutEnergyRamp:
                return ProgramPtr(new ProgramCombiner<
                        Acc,
                        TCP_IR3Program<Acc>,
                        VoltageProgram<Acc>
                    >(mAcc, steps));
            case AggressiveRamp:
                return ProgramPtr(new ProgramCombiner<Acc, AggressiveEnergyProgram<Acc>>(mAcc, steps));
            case Voltage:
                return ProgramPtr(new VoltageProgram<Acc>(mAcc, steps));
            case Energy:
                return ProgramPtr(new DefaultEnergyProgram<Acc>(mAcc, steps));
            default:
                throw except::ProgramTypeNotFound();
        }
    }

private:
    Acc& mAcc;
};


} // namespace stron

#endif
