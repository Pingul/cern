#ifndef SAMPLED_DISTRIBUTION
#define SAMPLED_DISTRIBUTION

#include <algorithm>
#include <vector>
#include <random>
#include <stdexcept>

#include <iostream>

template <typename T = double>
class Integral
{
public:
    struct MaxIterationsReached : public std::runtime_error { MaxIterationsReached() : std::runtime_error("") {} };

    using Func = T (*)(T);

    Integral(Func f, T threshold = 1e-5, int maxIterations = 15)
        : mF(f), mThresh(threshold), mMaxIter(maxIterations)
    {}

    T operator()(T a, T b)
    {
        if (b < a)
            return -operator()(b, a);
        else if (b == a)
            return 0;

        const T r = b - a;

        int n = 128; // power of 2
        std::vector<T> x(n);
        std::vector<T> y(n);
        for (int i = 0; i < n; ++i) {
            x[i] = a + r*T(i)/n;
            y[i] = mF(x[i]);
        }

        int iter = 0;
        T v, lv = 0;
        do {
            lv = v;
            v = 0;
            for (int i = 1; i < n; ++i)
                v += (y[i] + y[i - 1])/2*(x[i] - x[i - 1]);
            if (std::abs(v - lv) < mThresh)
                break;

            // double size and fill in gaps
            int old_n = n;
            n *= 2;
            x.resize(n);
            y.resize(n);

            for (int i = old_n - 1; i > 0; --i) {
                x[i*2] = x[i]; 
                y[i*2] = y[i];
            }
            for (int i = 1; i < n; i += 2) {
                x[i] = a + r*T(i)/n;
                y[i] = mF(x[i]);
            }
        } while (iter++ == 0 || iter < mMaxIter);

        if (iter >= mMaxIter) 
            throw MaxIterationsReached();
        return v;
    }
private:
    Func mF;
    T mThresh;
    int mMaxIter;
};

template <typename T = double, bool Interpolate = true>
class Sampled_distribution
{
public:
    struct InvalidBounds : public std::runtime_error { InvalidBounds() : std::runtime_error("") {} };
    struct CDFNotMonotonic : public std::runtime_error { CDFNotMonotonic() : std::runtime_error("") {} };

    using CDFFunc = T (*)(T);
    using PDFFunc = T (*)(T);

    Sampled_distribution(CDFFunc cdfFunc, T low, T high, unsigned resolution = 200) 
        : mLow(low), mHigh(high), mRes(resolution), mDist(0.0, 1.0)
    {
        if (mLow >= mHigh) throw InvalidBounds();

        mSampledCDF.resize(mRes + 1);
        const T cdfLow = cdfFunc(low);
        const T cdfHigh = cdfFunc(high);
        T last_p = 0;
        for (unsigned i = 0; i < mSampledCDF.size(); ++i) {
            const T x = i/mRes*(mHigh - mLow) + mLow;
            const T p = (cdfFunc(x) - cdfLow)/(cdfHigh - cdfLow); // normalising 
            if (! (p >= last_p)) throw CDFNotMonotonic();
            mSampledCDF[i] = Sample{p, x};
            last_p = p;
        }
    }

    template <typename Generator>
    T operator()(Generator& g) 
    {
        T cdf = mDist(g);
        auto s = std::upper_bound(mSampledCDF.begin(), mSampledCDF.end(), cdf);
        auto bs = s - 1;
        if (Interpolate && bs >= mSampledCDF.begin()) { 
            const T r = (cdf - bs->prob)/(s->prob - bs->prob);
            return r*bs->value + (1 - r)*s->value;
        }
        return s->value;
    }

private:
    const T mLow, mHigh;
    const double mRes;

    struct Sample { 
        T prob, value; 
        friend bool operator<(T p, const Sample& s) { return p < s.prob; }
    };

    std::vector<Sample> mSampledCDF;
    std::uniform_real_distribution<> mDist;
};

#endif
