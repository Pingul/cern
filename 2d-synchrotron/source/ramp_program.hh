#ifndef RAMP_PROGRAM
#define RAMP_PROGRAM

#include "common.hh"
#include "accelerator.hh"
#include "settings.hh"

namespace stron {

template <typename Accelerator>
class Program
{
public:
    using Ptr = std::unique_ptr<Program>;

    Program(Accelerator& acc)
        : mAcc(acc) {}
    virtual ~Program() {}
    virtual void restart() = 0;
    virtual void step() = 0;
protected:
    Accelerator& mAcc;
};

template <typename Accelerator, typename Program = Program<Accelerator>>
class EnergyRamp : public Program
{
    using T = typename Accelerator::ValType;

    size_t mIndex;
    std::vector<T> mEnergy;
public:
    EnergyRamp(Accelerator& acc, std::string energyProgram, int steps)
        : Program(acc), mIndex(0)
    {
        using common::skip;
        
        mEnergy.reserve(steps);
        std::cout << "Reading '" << energyProgram << "'..." << std::endl;
        std::ifstream file(energyProgram);
        if (!file.is_open())
            throw std::runtime_error("Could not open file");
        for (int i = 0; i < steps; ++i) {
            T data;
            file >> skip >> data;
            mEnergy.emplace_back(data*1e6);
        }
    }
    virtual void restart() override { mIndex = 0; Program::mAcc.setE(mEnergy[mIndex], true); }
    virtual void step() override { Program::mAcc.setE(mEnergy[++mIndex]); }
    virtual ~EnergyRamp() {}
};

template <typename Accelerator, typename Program = Program<Accelerator>>
class CollimatorRamp : public Program
{
    using T = typename Accelerator::ValType;
public:
    CollimatorRamp(Accelerator& acc, std::string motorProgram)
        : Program(acc)
    {
    }
    virtual ~CollimatorRamp() {}
};

template <typename Accelerator, typename EnergyRamp = EnergyRamp<Accelerator>>
class LHCRamp : public EnergyRamp
{
    using T = typename Accelerator::ValType;
public:
    LHCRamp(Accelerator& acc, int steps)
        : EnergyRamp(acc, LHC_RAMP_FILE, steps) {}
};

} // namespace stron

#endif
