// functions to calculate the priors

# pragma once
# include <iostream>
# include <cmath>
# include <Eigen/Dense>

using namespace std;
using namespace Eigen;

double psplitf(double alpha, double beta, int depth) {
	double psplit = alpha*pow(1 + depth, -beta);
	return psplit;
}

double nsplit(int iter, int spread) {
	if (iter > spread) 
		return spread;
	else
	 return iter;
}

double prule(double npred, double nsplits) {
	double prule = (1/(double)npred)*(1/(double)nsplits);
	return prule;
}

/*
int main() {
	double ps = psplit(0.95, 1, 0);
	double spr = nsplit(20, 10);
	double pr = prule(10, spr);
	cout << "The probability of splitting is: " << ps << endl;
	cout << "The prule splitting value is: " << pr << endl;
}
*/
