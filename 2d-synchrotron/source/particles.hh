#ifndef PARTICLES_HH
#define PARTICLES_HH

#include <string>
#include <vector>

//namespace stron {
namespace twodsynch {

template <typename T>
struct ParticleCollection 
{
    ParticleCollection(int n)
        : momentum(n), phase(n), x(n), px(n) {}

    std::vector<T> momentum;
    std::vector<T> phase;
    std::vector<T> x;
    std::vector<T> px;

    //void write(std::string filePath) const;
    size_t size() const { return momentum.size(); }

    // To avoid mistakes
    ParticleCollection(const ParticleCollection&) = delete;
    ParticleCollection(ParticleCollection&&) = delete;
};


} // namespace stron

#endif
