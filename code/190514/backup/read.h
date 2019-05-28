// A CLASS TO READ IN DATA

#pragma once
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
#include <Eigen/Dense>
#include "dist.h"

#include "boost/variant.hpp"

using namespace std;

// A class that can read a line of input,
// identifies the type of each field, and stores it, 
// with its name.
// Multiple lines can be stored, each as a map, indexed
// by a vector.

class Dataline {
	public:
		enum type {STRING_T, DOUBLE_T, INT_T};
		typedef boost::variant<string, double, int> variant_t;
		typedef std::pair<variant_t, type> value;
		typedef std::map<std::string, variant_t> line;
		typedef std::vector<line> dframe;
		typedef std::map<std::string, variant_t>::const_iterator line_it;
		typedef std::vector<line>::const_iterator frame_it;
		
		Dist Di;
		int p;
		std::vector<std::string> names;
		std::vector<int> catvec;
		std::map<string, string> allcats;
		std::map<string, vector<double>> dbounds;
		std::map<string, vector<int>> ibounds;
		std::map<string, type> alltypes;
		int ncats;
		dframe data;
		line dline;
		variant_t dvar;
		value dvalue;
		type dtype; 
		
		// Constructors
		Dataline() {};
		Dataline(int dl) : p(dl){};

		// A function to read the types of the variables
		void read_alltypes(std::string& input);

		// A function to read the bounds of the variables
		void read_bounds(std::string& input);
		
		// A function to read a line of data from a file
		line read_line(std::string& input);
		
		// A function to create variable names
		void make_names();
		
		// A container to store lines - i.e. dataframe
		void read_frame(std::ifstream& file, int nlines);

		// Functions to get particular types
		std::string get_string(const std::string& key,
				Dataline::line& fline);
		double get_double(const std::string& key,
				Dataline::line& fline);
		int get_int(const std::string& key,
				Dataline::line& fline);
		auto get_value(const std::string& key, 
				Dataline::line& fline);

		// A function to partition a set of categorical varaibles
//		void part_cats(string ssplit, string cats,
//			string lcat, string rcat);

		// Iterators of line type
		line_it get_key(const std::string& key, 
				Dataline::line& fline) const 
		{return fline.find(key);}
		line_it begin() const {return dline.begin();}
		line_it end() const {return dline.end();}
		VectorXd read_vec(string line, const int& DIM); 
		MatrixXd read_matbyrow(string line, const int& NDIM, const int& MDIM);
};

void Dataline::read_alltypes(std::string& input) {
	string word;
	int index = 0;

	std::istringstream instr(input);

	while(instr >> word) {
		if (word == "c" || word == "s") {
			alltypes[names[index]] = Dataline::STRING_T;
			catvec.push_back(0);
			index++;
		}
		else if (word == "d" || word == "cont") {
			alltypes[names[index]] = Dataline::DOUBLE_T;
			catvec.push_back(1);
			index++;
		}
		else {
			alltypes[names[index]] = Dataline::INT_T;
			catvec.push_back(2);
			index++;
		}

	}
}

void Dataline::read_bounds(std::string& input) {
	char * end;
	int index = 0;

	istringstream iss(input);
	vector<string> words{istream_iterator<string>{iss},
			istream_iterator<string>{}};
	//for(int i = 0; i < words.size(); i++) {
	//	cout << " " << words[i];
	//}
	cout << endl;
	for(pair<string, type> item : alltypes) {
		if(item.second == Dataline::STRING_T) {
			allcats[item.first] = words[index]; 
			index++;
		}
		else if(item.second == Dataline::DOUBLE_T) {
			size_t pos = 0;
			std::string token;
			dbounds[item.first];
			auto it = dbounds.find(item.first);
			istringstream ws(words[index]);
			while(getline(ws, token, ',')) {
				double num = strtod(token.c_str(), 
						&end);
				it->second.push_back(num);
			}	
			index++;
		}
		else {
			size_t pos = 0;
			std::string token;
			ibounds[item.first];
			auto it = ibounds.find(item.first);
			istringstream ws(words[index]);
			while(getline(ws, token, ',')) {
				int num = stoi(token);
				it->second.push_back(num);
			}	
			index++;

		}
	}
}

Dataline::line Dataline::read_line(std::string& input) {
	Dataline::line fline;

	char * end;
	int index = 0;

	istringstream iss(input);
	vector<string> words{istream_iterator<string>{iss},
			istream_iterator<string>{}};

	for(pair<string, type> item : alltypes) {
		if(item.second == Dataline::STRING_T) {
			fline[item.first] = words[index];
			index++;
		}
		else if(item.second == Dataline::DOUBLE_T) {
			double num = strtod(words[index].c_str(), 
						&end);
			fline[item.first] = num;
			index++;
		}
		else {
			int num = stoi(words[index]);
			fline[item.first] = num;
			index++;
		}
	}
	return fline;
}

// A function that makes a vector of lines (matrix of maps) 
// containing data.

void Dataline::read_frame(std::ifstream& file, int nlines = 10) {
	string newline; 
	int i = 0;
	while(i < nlines) {
		while(getline(file, newline)) {
			data.push_back(read_line(newline));
			i++;
		}
		if (!getline(file, newline)) {
			break;
		}
	}
}

VectorXd Dataline::read_vec(string line, const int& DIM) {
	istringstream bss(line);
	VectorXd outvec(DIM);
	double val;
	for (int i = 0; i < DIM; ++i) {
		bss >> val;
		outvec(i) = val;
	}
	return outvec;
}

MatrixXd Dataline::read_matbyrow(string line, const int& NDIM, const int& MDIM) {
	istringstream bss(line);
	MatrixXd outmat(NDIM, MDIM);
	double val;
	for(int i = 0; i < NDIM; ++i) {
		for(int j = 0; j < MDIM; ++j) {
			bss >> val;
			outmat(i, j) = val;
		}
	}
	return outmat;
}

std::string Dataline::get_string(const std::string& key, 
		Dataline::line& fline) {
	auto it = fline.find(key);
	return boost::get<std::string>(it->second);
}
double Dataline::get_double(const std::string& key, 
		Dataline::line& fline) {
	auto it = fline.find(key);
	return boost::get<double>(it->second);
}
int Dataline::get_int(const std::string& key, 
		Dataline::line& fline) {
	auto it = fline.find(key);
	return boost::get<int>(it->second);
}
auto Dataline::get_value(const std::string& key, Dataline::line& fline) {
	line_it it = get_key(key, fline);
	return it->second;
}

// A function that makes varaible names Xi based on number of predictors
void Dataline::make_names() {
	for (int i = 1; i <= p; i++) {
	std::string pred = "X"+to_string(i);
	Dataline::names.push_back(pred);
	}	
}

// Get the index of a string in a vector of strings
int get_pred_index(string pred, vector<string> names) {
	int index = 0;
    for (int i = 0; i < names.size(); i++) {
		if (pred == names[i]) {
			index = i;
		}
	}
	return index;
}
/*
int main() {
	ifstream tin("testdata.txt");
	//string tst;
	//getline(tin, tst);
	Dataline rdata(9);
	rdata.make_names();
	string firstline;
	getline(tin, firstline);
	rdata.read_alltypes(firstline);
	string secondline;
	getline(tin, secondline);
	rdata.read_bounds(secondline);
	for(char t : secondline) {
		cout << t << " ";
	}
	cout << endl;

	for(auto t : rdata.allcats) {
		cout << "second: " << t.second << " ";
	}
	cout << endl;
	
	for(auto t : rdata.dbounds) {
		cout << "first: " << t.second[0] << " second: " << t.second[1] << " ";
	}
	cout << endl;

	for(auto t : rdata.ibounds) {
		cout << "first: " << t.second[0] << " second: " << t.second[1] << " ";
	}
	cout << endl;

	for (int i = 0; i < 10; ++i) {
		string newline;
		getline(tin, newline);
		Dataline::line nline;
		nline = rdata.read_line(newline);
		for(auto l : nline) {
		       cout << l.first << " " 
			       << rdata.get_value(l.first, nline) << " ";
		}
		cout << endl;
	}
}
	//rdata.dline = rdata.read_line(tst, 9);
	//rdata.read_frame(tin);
//	rdata.dline = rdata.data[2];
//	auto it = rdata.dline.begin();
//	it = it++;
//	cout << it->second.first << endl;
	//rdata.dline = rdata.data[0];
	//rdata.get_cats();

	cout << rdata.get_string("X6", rdata.data[0]) << endl;
	cout << rdata.get_string("X6", rdata.data[1]) << endl;
	cout << rdata.get_string("X6", rdata.data[2]) << endl;
	cout << rdata.ncats << endl;
	cout << rdata.data.size() << endl;

	vector<string> tcats = {"a", "b", "c", "d"};
	string tssplit = "a";
	vector<string> tlcats;
	vector<string> trcats;
	rdata.part_cats(tssplit, tcats, tlcats, trcats);
	cout << "lcatsize " << tlcats.size() << endl;
	for(int i = 0; i < tlcats.size(); i++) {
		if(tlcats.empty()) {cout << "empty";}
		else {cout << "lcat " << tlcats[i] << " ";}
	}
	cout << endl;
	cout << "rcatsize " << trcats.size() << endl;
	for(int i = 0; i < trcats.size(); i++) {
		if(trcats.empty()) {cout << "empty";}
		else {cout << "rcat " << trcats[i] << " ";}
	}
	cout << endl;
}
*/
