#ifndef SAMPLED_DISTRIBUTION
#define SAMPLED_DISTRIBUTION

#include <algorithm>
#include <vector>
#include <random>
#include <stdexcept>

/*
 * Given a CDF (un-normalised/normalised), generate values according to that distribution
 *
 * The class samples the given CDF to calculate the inverse.
 *
 * Use as:
 *  auto f = [](double x){ return x^/2; } // linear distribution; integral of f(x) = x.
 *  std::mt19937 gen;
 *  Sampled_distribution<> dist(f, 0.0, 1.0);
 *  double v = dist(gen);
 */
template <typename T = double>
class Sampled_distribution
{
public:
    using CDFFunc = T (*)(T);

    Sampled_distribution(CDFFunc cdfFunc, T low, T high, unsigned resolution = 200) 
        : mLow(low), mHigh(high), mRes(resolution), mDist(0.0, 1.0)
    {
        if (mLow >= mHigh) throw InvalidBounds();

        mSampledCDF.resize(mRes);
        const T cdfLow = cdfFunc(low);
        const T cdfHigh = cdfFunc(high);
        for (int i = 0; i < mSampledCDF.size(); ++i) {
            const T x = (i + 1)/double(mRes)*(mHigh - mLow) + mLow;
            const T prob = (cdfFunc(x) - cdfLow)/(cdfHigh - cdfLow); // normalising
            mSampledCDF[i] = Sample{prob, x};
            std::cout << i << " " << mSampledCDF[i].prob << " " << mSampledCDF[i].value << std::endl;
        }
    }

    template <typename Generator>
    T operator()(Generator& g) 
    {
        T cdf = mDist(g);
        auto search = std::upper_bound(mSampledCDF.begin(), mSampledCDF.end(), cdf);
        return search->value;
    }

private:
    struct InvalidBounds : public std::runtime_error { InvalidBounds() : std::runtime_error("") {} };

    const T mLow, mHigh;
    const unsigned mRes;
    
    struct Sample { 
        T prob, value; 
        friend bool operator<(T p, const Sample& s) { return p < s.prob; }
    };

    std::vector<Sample> mSampledCDF;
    std::uniform_real_distribution<> mDist;
};

#endif
