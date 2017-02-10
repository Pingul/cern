#include "timer.hh"
#include <iostream>

namespace common
{

template <bool WithMessage>
void MessageTimer<WithMessage>::start(const std::string message)
{
	if (WithMessage) std::cout << message;
	mStart = std::chrono::high_resolution_clock::now();
}

template <bool WithMessage>
double MessageTimer<WithMessage>::stop() const
{
	// This should really stop the clock as well, but is unecessary. 'stop'ing a clock without a 'start' is undefined behaviour.
	double t = check();
	if (WithMessage) std::cout << "..." << t << " ms." << std::endl;
	return t;
}

template <bool WithMessage>
double MessageTimer<WithMessage>::restart(const std::string message)
{
	double t = stop();
	start(message);
	return t;
}

template <bool WithMessage>
double MessageTimer<WithMessage>::check() const
{
	TimePoint end = std::chrono::high_resolution_clock::now();
	double t = std::chrono::duration_cast<std::chrono::microseconds>(end - mStart).count()/1000.0;
	return t;
}

}; // namespace common

template class common::MessageTimer<true>;
template class common::MessageTimer<false>;

// #include <thread>

// void f(int n)
// {
// 	std::this_thread::sleep_for(std::chrono::seconds(n));
// }


// int main()
// {
// 	common::SilentTimer t;
// 	t.start("Test");
// 	f(1);
// 	t.restart("Test2");
// 	f(2);
// 	t.stop();
// }
