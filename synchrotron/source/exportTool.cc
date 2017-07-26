#include <iostream>
#include <vector>
#include <string>
#include <tbb/task_scheduler_init.h>
#include <sys/stat.h>

#include "accelerator.hh"
#include "particles.hh"
#include "synchrotron.hh"
#include "ramp_program.hh"

// Arguments
std::string outputFormat;
std::string outputDirectory;
int nbrOutputFiles;
int pPerFile;
std::string longDist;
std::string transDist;
int nbrTurnsSimulate;

void printInput()
{
    std::cout 
        << "Received paramters" << std::endl
        << "\toutputFormat : " << outputFormat << std::endl
        << "\toutputDirectory : " << outputDirectory << std::endl
        << "\tnbrOutputFiles : " << nbrOutputFiles << std::endl
        << "\tpPerFile : " << pPerFile << std::endl
        << "\tlongDist : " << longDist << std::endl
        << "\ttransDist : " << transDist << std::endl
        << "\tnbrTurnsSimulate : " << nbrTurnsSimulate << std::endl;
}

bool validateArguments() 
{
    std::vector<std::string> errors;
    if (outputFormat != "coll" && outputFormat != "nocoll")
        errors.emplace_back("outputFormat: has to be 'coll' or 'nocoll'");
    
    struct stat sb;
    if (!(stat(outputDirectory.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
        errors.emplace_back("outputDirectory: is not a valid directory.");

    if (!(nbrOutputFiles > 0))
        errors.emplace_back("nbrOutputFiles: must be > 0");
    
    if (!(pPerFile > 0 && pPerFile % 2 == 0))
        errors.emplace_back("pPerFile: must be > 0 and even");

    if (nbrTurnsSimulate < 0)
        errors.emplace_back("nbrTurnsSimulate: must be > 0");

    if (longDist != "lin+exp" && longDist != "constant")
        errors.emplace_back("longDist: must be 'lin+exp' or 'constant', is '" + longDist + "'");

    if (errors.size() > 0)  {
        std::cout << "Errors:" << std::endl;
        for (auto& e : errors) std::cout << "\t" << e << std::endl;
    }
    
    std::cout << "transDist not used" << std::endl;

    return errors.size() == 0;
}

// The return 'Acc' here is a little bit of a hack -- SimpleSynchrotron is initialized by value
template <typename Acc, typename ParticlesPtr>
Acc simulateTurns(Acc& acc, ParticlesPtr& p, int turns)
{
    if (turns < 1) return acc;
    stron::SimpleSynchrotron<double> ss(acc);
    stron::ProgramGenerator<Acc> progGen(ss.getAcc());

    ss.addParticles(p);
    ss.simulateTurns(progGen.create(turns, stron::ProgramType::LHCRamp));

    // make sure the particles are between 0 and 2π so we can import into SixTrack
    for (auto& phi : p->phase) {
        auto r = static_cast<int>(std::floor(phi/(2*cnst::pi)));
        phi -= 2.0*cnst::pi*r;
    }
    return ss.getAcc();
}

void generateFiles()
{
    using Acc = typename stron::Accelerator<double>;
    using PColl = particles::ParticleCollection<double>;
    using ParticleGenerator = particles::ParticleGenerator<Acc>;

    Acc lhc = Acc::getLHC();
    ParticleGenerator pgen(lhc);

    void (*exportFunc)(const Acc&, const PColl&, const std::string&) = particles::sixtrackExport_nonCollimat;
    if (outputFormat == "coll") exportFunc = particles::sixtrackExport;


    typename ParticleGenerator::PCollPtr all_particles;
    auto n = nbrOutputFiles*pPerFile;
    if (longDist == "lin+exp")
        all_particles = pgen.create(n, particles::CCombined, particles::DoubleGaussian);
    else if (longDist == "constant")
        all_particles = pgen.create(n, particles::CConstant, particles::DoubleGaussian);
    else
        throw std::runtime_error("Unexpected 'longDist' argument");
    lhc = simulateTurns(lhc, all_particles, nbrTurnsSimulate);

    exportFunc(lhc, *all_particles, outputDirectory + "/all.txt");
    for (int file_i = 0; file_i < nbrOutputFiles; ++file_i) {
        std::vector<int> indices(pPerFile);
        int c = 0;
        for (int& i : indices) i = file_i*pPerFile + c++;

        PColl part(*all_particles, indices);
        exportFunc(lhc, part, outputDirectory + "/" + std::to_string(file_i + 1) + ".txt");
    }
}

/*
    Input format:
        1. 'coll' or 'nocoll'
            Exports for SixTrack collimation version or non-collimation version
        2. Directory file path for output (relative or absolute)
            E.g.: 'relative/path', '/absolute/path'
        3. Number of output files
        4. Number of particles per file (should be an even number)
        5. Longitudinal distribution
        6. Horizontal (transversal?) distribution
        7. Numbr of turns to simulate
*/
int main(int argc, char* argv[])
{
    tbb::task_scheduler_init init;

    std::vector<std::string> args(argv, argv + argc);

    if (args.size() != 8)
        throw std::runtime_error("Expected 7 arguments, received " + std::to_string(args.size()));
    
    outputFormat = args[1];
    outputDirectory = args[2];
    nbrOutputFiles = std::stoi(args[3]);
    pPerFile = std::stoi(args[4]);
    longDist = args[5];
    transDist = args[6];
    nbrTurnsSimulate = std::stoi(args[7]);

    printInput();
    if (!validateArguments()) {
        std::cout << "Errors in input, try again" << std::endl;
        return 1;
    }
    generateFiles();
}
