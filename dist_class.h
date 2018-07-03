// DISTRIBUTIONS USED

# pragma once
# include <iostream>
# include <cmath>
# include <Eigen/Dense>
# include <random>

using namespace std;
using namespace Eigen;


class Dist {
	int seed;
	public:
	Dist(){};
	Dist(int seednum) : seed(seednum){};
	int runif(const int min, const int max);
	MatrixXd mvn(int n, const VectorXd& mu, const MatrixXd& sig);
};	

// Sample from a uniform distribution
int Dist:: runif(const int min, const int max) {
	random_device rand_dev;
	mt19937 generator(rand_dev());
	uniform_int_distribution<int> distr(min, max);
	int runif = distr(generator);
	return runif;
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
/*
int main() {
	Dist rn(2233);
	int rnum = rn.runif(0, 10);
	cout << rnum << endl;
	Dist d;
	int i = 0;
	while(i<10){
		int rnum1 = d.runif(0, 10);
		cout << rnum1 << endl;
		i++;
	}
}*/
