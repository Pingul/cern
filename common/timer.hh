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

}; // namespace jwc

#endif