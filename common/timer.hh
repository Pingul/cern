#ifndef TIMER_HH
#define TIMER_HH

#include <chrono>
#include <string>
#include <iostream>

namespace jwc {

template <bool WithMessage>
class MessageTimer
{
public:
	MessageTimer() = default;
	~MessageTimer() = default;

	void start(const std::string message = "");
	double restart(const std::string message = "");
	double check() const;
	double stop() const;

	MessageTimer(const MessageTimer& other) = delete;
private:
	using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
	TimePoint mStart;
};

typedef MessageTimer<true> Timer;
typedef MessageTimer<false> SilentTimer;


// Helper functionality for the timer
struct TimeDuration { unsigned d, h, m, s, ms; };
// std::ostream& operator<<(std::ostream& os, )

inline TimeDuration MillisecondsToDuration(unsigned ms)
{
	TimeDuration td;
	const unsigned DAY = 24*3600*1000;
	const unsigned HOUR = 3600*1000;
	const unsigned MIN = 60*1000;
	const unsigned SEC = 1000;
	td.d  = ms/DAY;
	td.h  = (ms - td.d*DAY)/HOUR;
	td.m  = (ms - td.d*DAY - td.h*HOUR)/MIN;
	td.s  = (ms - td.d*DAY - td.h*HOUR - td.m*MIN)/SEC;
	td.ms = (ms - td.d*DAY - td.h*HOUR - td.m*MIN - td.s*SEC);
	return td;
}


}; // namespace jwc

#endif