// A function to read in data

//#pragma once
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <map>
#include <utility>
#include <vector>
#include <type_traits>
#include <cstdlib>

#include "boost/variant.hpp"

using namespace std;

// A class that can:
// preprocess the predictors by:
// getting its Name
// identifing the type of each field, 
// the range of its values (min, max) or elements if categorical
// and storing these in a container. 
// 
// Multiple lines can be stored, each as a map, indexed
// by a vector.

class Dataline {
	private:
		enum type { STRING_T, LDOUBLE_T, NONE_T };
		typedef boost::variant<string, long double> variant_t;
	public:
		typedef std::pair<variant_t, type> value;
		typedef std::map<std::string, value> line;
		typedef std::vector<line> dframe;
		typedef std::map<std::string, value>::const_iterator map_it;
		typedef std::vector<line>::const_iterator vec_it;
