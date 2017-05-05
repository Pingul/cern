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

    Integral(Func f, T threshold = 1e-5, int maxIterations = 20)
        : mF(f), mThresh(threshold), mMaxIter(maxIterations)
    {}

    T operator()(T a, T b)
    {
        if (b < a)
            return -operator()(b, a);
        else if (b == a)
            return 0;

        const T r = b - a;

        int n = 1;
        int iter = 0;

        T sum = mF((a + b)/2);
        T fl = 0.5*(mF(a) + mF(b));
        T lv, v = r/n*(fl + sum);
        do {
            lv = v;
            n *= 2;

            for (int i = 1; i < n; i += 2)
                sum += mF(a + i*r/n);
            v = r/n*(fl + sum);

            if (std::abs(v - lv) < mThresh)
                break;
        } while (iter++ == 0 || iter < mMaxIter);
        //std::cout << "Iterations: " << iter << std::endl;
        //std::cout << "Threshold: " << std::abs(v - lv) << std::endl;

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
    using Func = T (*)(T);

    // PDF requires numerical integration -- better performance and result if possible to provide analytic CDF instead
    enum DistType { PDF, CDF };

    Sampled_distribution(Func f, T low, T high, DistType dtype = CDF, unsigned resolution = 200) 
        : mL(low), mH(high), mR(resolution), mDist(0.0, 1.0)
    {
        if (mL >= mH) throw InvalidBounds();

        mSampledCDF.resize(mR + 1);
        const T fl = dtype == CDF ? f(mL) : 0;
        const T fh = dtype == CDF ? f(mL) : Integral<T>(f)(mL, mH);
        T last_p = 0;
        for (unsigned i = 0; i < mSampledCDF.size(); ++i) {
            const T x = i/mR*(mH - mL) + mL;

            T cdf;
            if (dtype == CDF)
                cdf = f(x);
            else
                cdf = Integral<T>(f)(mL, x); 
            const T p = (cdf - fl)/(fh - fl); // normalising 

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
    const T mL, mH;
    const double mR;

    struct Sample { 
        T prob, value; 
        friend bool operator<(T p, const Sample& s) { return p < s.prob; }
    };

    std::vector<Sample> mSampledCDF;
    std::uniform_real_distribution<> mDist;
};

#endif
