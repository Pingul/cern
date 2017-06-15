#include <iostream>
#include <vector>
#include <string>
#include <tbb/task_scheduler_init.h>
#include <sys/stat.h>

#include "accelerator.hh"
#include "particles.hh"

// Arguments
std::string outputFormat;
std::string outputDirectory;
int nbrOutputFiles;
int pPerFile;
std::string longDist;
std::string transDist;

void printInput()
{
    std::cout 
        << "Received paramters" << std::endl
        << "\toutputFormat : " << outputFormat << std::endl
        << "\toutputDirectory : " << outputDirectory << std::endl
        << "\tnbrOutputFiles : " << nbrOutputFiles << std::endl
        << "\tpPerFile : " << pPerFile << std::endl
        << "\tlongDist : " << longDist << std::endl
        << "\ttransDist : " << transDist << std::endl;
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

    if (errors.size() > 0)  {
        std::cout << "Errors:" << std::endl;
        for (auto& e : errors) std::cout << "\t" << e << std::endl;
    }
    
    std::cout << "longDist and transDist not used" << std::endl;

    return errors.size() == 0;
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

    auto all_particles = pgen.create(nbrOutputFiles*pPerFile, particles::CCombined, particles::DoubleGaussian);
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
        3. # output files
        4. # particles per file (should be an even number)
        5. Longitudinal distribution
        6. Horizontal (transversal?) distribution
*/
int main(int argc, char* argv[])
{
    tbb::task_scheduler_init init;

    std::vector<std::string> args(argv, argv + argc);

    if (args.size() != 7)
        throw std::runtime_error("Expected 6 arguments, received " + std::to_string(args.size()));
    
    outputFormat = args[1];
    outputDirectory = args[2];
    nbrOutputFiles = std::stoi(args[3]);
    pPerFile = std::stoi(args[4]);
    longDist = args[5];
    transDist = args[6];

    printInput();
    if (!validateArguments()) {
        std::cout << "Errors in input, try again" << std::endl;
        return 1;
    }
    generateFiles();
}
