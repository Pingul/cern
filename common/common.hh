#include <istream>
#include <string>

/*
 * Makes it possible to write
 *
 *   T a, b;
 *   std::cin >> skip >> a >> b;
*/
template <typename T = std::string, typename CharT = char>
std::basic_istream<CharT>& skip(std::basic_istream<CharT>& in) 
{
	T ignoredValue;
	return in >> ignoredValue;
}

