// DISTRIBUTIONS USED

# pragma once
# include <iostream>
# include <cmath>
# include <Eigen/Dense>
# include <random>
# include <algorithm>
# include <vector>
# include <map>
# include <iterator>

using namespace std;
using namespace Eigen;


class Dist {
	int seed;
	public:
	Dist(){};
	Dist(int seednum) : seed(seednum){};
	int runif(const int min, const int max);
	double cunif(const double min, const double max);
	double rnorm(const double mu, const double var);
	double rgamma(const double alphg, const double betg);
	MatrixXd mvn(int n, const VectorXd& mu, const MatrixXd& sig);
	char rchar(const string& cats);
	pair<string, string> part_cats(const char& ssplit, string cats);

	// Uniform Sample from STL container. 
	// Usage: random_selector<> selector{"container"}
	// Returns a reference to container value type.
	template <typename RandomGenerator = std::default_random_engine>
	struct random_selector { 
		random_selector(RandomGenerator g = 
				RandomGenerator(std::random_device()()))
			: gen(g) {}

		template <typename Iter>
		Iter select(Iter start, Iter end) {
			std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
			std::advance(start, dis(gen));
			return start;
		}

		//convenience function
		template <typename Iter>
		Iter operator()(Iter start, Iter end) {	
			return select(start, end);
		}

		//convenience function that works on anything with a sensible begin() and end()
		//Returns with a ref to the value type
		template <typename Container>
		auto operator()(const Container& c) -> decltype(*begin(c))& {
			return *select(begin(c), end(c));
		}
		
		private:
			RandomGenerator gen;
	};
};	

// Sample from a discrete UNIFORM distribution
int Dist:: runif(const int min, const int max) {
	random_device rand_dev;
	mt19937 generator(rand_dev());
	uniform_int_distribution<int> distr(min, max);
	int runif = distr(generator);
	return runif;
}

// Sample from a UNIFORM distribution
double Dist:: cunif(const double min, const double max) {
	random_device rand_dev;
	mt19937 generator(rand_dev());
	uniform_real_distribution<double> distr(min, max);
	double cunif = distr(generator);
	return cunif;
}

// Sample from a NORMAL distribution
double Dist:: rnorm(const double mu, const double var) {
	random_device rand_dev;
	mt19937 generator(rand_dev());
	normal_distribution<double> distr(mu, var);
	double rnorm = distr(generator);
	return rnorm;
}

// Sample from a GAMMA distribution
double Dist:: rgamma(const double alphg, const double betg) {
	random_device rand_dev;
	mt19937 generator(rand_dev());
	gamma_distribution<double> distr(alphg, betg);
	double rgamma= distr(generator);
	return rgamma;
}

// Sample from the multivariate normal distribution
// Assume convariance matrix is positive semi-definite
MatrixXd Dist:: mvn(int n, const VectorXd& mu, const MatrixXd& sig) {
	int DIMS = mu.size();
	MatrixXd outmat(DIMS, n);
	LLT<MatrixXd> llt(sig);
	MatrixXd l = llt.matrixL();
	VectorXd z(DIMS);
	VectorXd x(DIMS);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < DIMS; i++) {
			double u1 = runif(0, RAND_MAX)/(double)RAND_MAX;
			double u2 = runif(0, RAND_MAX)/(double)RAND_MAX;
			double z0 = pow((-2*log(u1)), 0.5)*cos(2*M_PI*u2);
			z(i) = z0;
		}
		x = mu + (l*z);
		outmat.col(j) = x;
	}
	return outmat;
}

// Sample from a string
char Dist:: rchar(const string& cats) {
	int min = 0;
	int max = cats.length()-1;
	random_device rand_dev;
	mt19937 generator(rand_dev());
	uniform_int_distribution<int> distr(min, max);
	int rint = distr(generator);
	char split = cats[rint];
	return split;
}

// Random partition of string given a letter
pair<string, string> Dist::part_cats(const char& ssplit, string cats) {
	string sub = cats;
	string token;
	vector<string> tokens;
	random_shuffle(sub.begin(), sub.end());
	istringstream ins(sub);
	for(int i = 0; i<2; ++i) {
		getline(ins, token, ssplit);
		if(token.length() == 0) {
			token = "";
		}
		tokens.push_back(token);
	}
	pair<string, string> outpair = std::make_pair(tokens[0], tokens[1]);
	return outpair;
}
/*
int main() {
	Dist Di;
	string test = "abcde";
	for(char t : test)  {
		cout << t << " ";
	}
	cout << endl;
	char tsplit = Di.rchar(test);
	cout << tsplit << endl;

	pair<string, string> lrcats;
	lrcats = Di.part_cats(tsplit, test);

	cout << "first " << lrcats.first << " second " << lrcats.second << endl;

}

	vector<int> test = {10,20,30,40,50,60};
	Dist::random_selector<> selector{};
	auto rind = selector(test);
	cout << "value: " << rind << endl;

	map<int, int> test1;
	for(int i = 0; i < test.size(); i++) {
		test1[i] = test[i];
	}

	auto rind1 = selector(test1);
	cout << "value2: " << rind1.second << endl;

	Dist rn(2233);
	int rnum = rn.runif(0, 10);
	cout << rnum << endl;
	Dist d;
	int i = 0;
	while(i<10){
		double rnum1 = d.rnorm(0, 1);
		cout << rnum1 << endl;
		i++;
	}
	VectorXd test(1);
	MatrixXd sigt(2,2);
	sigt << 1,0,0,1;
	VectorXd mut(2);
	mut << 0,0;
	for (int j = 0; j<10; j++) {
		test = d.mvn(1, mut, sigt);
		cout << test.transpose() << endl;
	}
}*/
