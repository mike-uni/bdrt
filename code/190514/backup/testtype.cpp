// Test data types
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <iterator>
#include <type_traits>
#include <cstdlib>
#include <set>
#include "boost/variant.hpp"
#include "dist.h"
#include "read.h"
#include <Eigen/Dense>
#include <functional>

using namespace std;

template <class T>
inline bool isPrimitiveType(const T& data) {
	return std::is_fundamental<T>::value;
}

vector<pair<string, int>>::const_iterator prev_it(const string& pred, 
		const vector<pair<string, int>>& prevpreds) {
	vector<pair<string, int>>::const_reverse_iterator rit;
	for(rit = prevpreds.rbegin(); 
			rit != prevpreds.rend(); ++rit) {
		cout << rit->first << endl;
		if (rit->first == pred) {
			return (rit+1).base();}
	}
		return prevpreds.end();
}

vector<int> parent_pairs(vector<int> leaves) {
	set<int>lparents;
	set<int>rparents;
	for(int i = 0; i<leaves.size(); ++i) {
		if(leaves[i]%2 == 0) {
			lparents.insert(leaves[i]/int(2));
		}
		else {
			rparents.insert(leaves[i]/int(2));
		}
	}
	vector<int> parent_pairs;
	set_intersection(lparents.begin(), lparents.end(),
			rparents.begin(), rparents.end(),
			back_inserter(parent_pairs));
	return parent_pairs;
}

void mackey(Ref<MatrixXd>sigpup, const MatrixXd& sigpmat, const double& alpha, const double& beta, const double& gamma) {
	MatrixXd FMAT(2,2); 
	FMAT << (1-gamma), 0, 0, beta;
	VectorXd temp(2);
	double ztemp;
	for (int i = 0; i < 5; ++i) {
		ztemp = sigpmat(1,i);
		ztemp = ztemp/(double)(1+pow(ztemp, alpha));
		temp << sigpmat(0,i), ztemp; 
		sigpup.col(i) << FMAT*temp;
	}
}	

double mackeysim(const double& xt, const double& xtm, const double& alpha, const double& beta,
		const double& gamma, const double& V) {
	double xt1 = xt*(1-gamma) + beta*(xtm/(double)pow(xtm, alpha));
	return xt1;
}

int main() {
	double alpha = 10;
	double beta = 0.2;
	double gamma = 0.1;
	double gamma2 = 1.4;
	VectorXd zhat(2);
	zhat << 1,1;
	MatrixXd SMAT(2,2); 
	SMAT << 1,0.5,0.5,1;

	cout << SMAT.transpose()*zhat << endl;
	cout << SMAT*zhat.transpose() << endl;

	LLT<MatrixXd> chSMAT(SMAT);
	chSMAT.rankUpdate((SMAT.col(0) - SMAT.col(1)), 2);
	MatrixXd LSMAT = chSMAT.matrixL();
	cout << "LSMAT " << endl;
	cout << LSMAT << endl;
	LSMAT = gamma2*LSMAT;
	cout << "LSMAT " << endl;
	cout << LSMAT << endl;
	MatrixXd sigpu(2, 5);
	MatrixXd sigmt(2, 5);
	sigmt.col(0) << zhat;
	sigmt.block(0,1,2,2) << LSMAT.colwise() + zhat;
	cout << "sigmt " << endl;
	cout << sigmt << endl;
	sigmt.block(0,3,2,2) << (-LSMAT).colwise() + zhat;
	cout << "sigmt " << endl;
	cout << sigmt << endl;
	mackey(sigpu, sigmt, alpha, beta, gamma);
	cout << "sigmt " << endl;
	cout << sigmt << endl;
	cout << "sigpu " << endl;
	cout << sigpu << endl;
}
	/*
	Dataline D;
	string tests2 = "1 2 3 4";
	istringstream ins(tests2);
	istringstream ins1(tests2);
	Eigen::VectorXd vec;
	Eigen::MatrixXd mat;

	vec = D.read_vec(ins, 2);
	mat = D.read_matbyrow(ins1, 2);

	cout << vec << endl;
	cout << mat << endl;
	Dist Di;
	for(int i = 0; i < 20; ++i) {
		int rand_int = Di.runif(0,2);
		cout << rand_int << " ";
	}
	cout << endl;
	for(int i = 0; i < 20; ++i) {
		double rand_dub = Di.cunif(0,1);
		cout << rand_dub << " ";
	}
	cout << endl;
	vector<int>leaves = {4,6,10,11,14,30,31};
	vector<int>ppairs1 = parent_pairs(leaves);
	//set<int> lparent {2,3,5,7,15};
	//set<int> rparent {5,15};

	vector<int>leaves2 = {7,8,9,11,12,20,21,26,27};
	vector<int>ppairs2 = parent_pairs(leaves2);
	// lparent = {4,6,10,13}
	// rparent = {3,4,5,10,13}
	
	for(auto s : ppairs1) {
		cout << s << " ";
	}
	cout << endl;

	for(auto s : ppairs2) {
		cout << s << " ";
	}
	cout << endl;
	vector<pair<string, int>> test;
	for(int i = 0; i <= 5; ++i) {
		pair<string, int> mypair;
		mypair = make_pair("X"+to_string(i%4+1), i+1);
		test.emplace_back(mypair);
		cout << test[i].first << " " << test[i].second << " ";
	}
	cout << endl;

	string tpred = "X"+to_string(6);

	vector<pair<string, int>>::const_iterator tit;

	tit = prev_it(tpred, test);
	if (tit != test.end()) {
		cout << tit->first << " " << tit->second << endl;
	}
	else {cout << "This is the end." << endl;}

	ifstream infile("testdata.txt");
	string line;
	string word;
	boost::variant<string, int, long double> var; 
	char* end = NULL;
	
	// getline reads a line of data as stores it as a string
	while (getline(infile, line)) {	
		istringstream inl(line);
			while(inl >> word) {
			cout << setprecision(17);
			cout << "original word " << word << endl;
			long double f = strtold(word.c_str(), &end);
			var = f;
			cout << "original f " << f << endl;
			if(*end) {
				var = word;
			}
			cout << var << endl;
			}
		cout << " This is the line." << endl;
	}
	*/


/*	boost::variant<int, double> var;
	int num = 1;
	double dble;
	var = 1;

	cout << boolalpha << isPrimitiveType(boost::get<int>(var)) << endl;
}

#include <map>
#include <algorithm>
#include <string>
#include <memory>
#include <functional>
#include <iostream>

void setRecursion( int dint ){} 
void setDepth( float ){} 
void statusReport( std::string job_name ){}

template<typename R> R readA(std::istream& stream){
	R   value;
	stream >> value;
	return value;
}

int main() {
	std::map<std::string, std::function<void(std::istream&)> > const  action {
		{"enable_recursion", [](std::istream& value){setRecursion(readA<bool>(value));}},
		{"max_search_depth", [](std::istream& value){setDepth(readA<int>(value));}},
		{"job_name",         [](std::istream& value){statusReport(readA<std::string>(value));}}
};

	std::string     key;
	while(std::cin >> key) {
		auto find = action.find(key);
		if (find != action.end()) { 
			find->second(std::cin);
		}
		else {
		}
		std::string ignoreLine;
		std::getline(std::cin, ignoreLine);
	}
}
*/
