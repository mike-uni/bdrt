// functions to calculate the priors

# pragma once
# include <iostream>
# include <cmath>

using namespace std;

double psplitf(double alpha, double beta, int depth) {
	double psplit = alpha*(1/double(pow(1 + depth, beta)));
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
