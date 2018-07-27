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
	double rnorm(const double mu, const double var);
	double rgamma(const double alphg, const double betg);
	MatrixXd mvn(int n, const VectorXd& mu, const MatrixXd& sig);
};	

// Sample from a UNIFORM distribution
int Dist:: runif(const int min, const int max) {
	random_device rand_dev;
	mt19937 generator(rand_dev());
	uniform_int_distribution<int> distr(min, max);
	int runif = distr(generator);
	return runif;
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
/*
int main() {
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
