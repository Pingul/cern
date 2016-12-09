#include <istream>
#include <string>

template <typename T = std::string, typename CharT = char>
std::basic_istream<CharT>& skip(std::basic_istream<CharT>& in) 
{
	T ignoredValue;
	return in >> ignoredValue;
}

