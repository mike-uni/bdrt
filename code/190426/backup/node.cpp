// NODE CLASS

#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <limits>
#include <map>
#include "read.h"

#include "boost/variant.hpp"

using namespace std;
using namespace Eigen;

class Node {
	public:
		// For Doubles
		double dsplit = 0;
		double min = 0;
		double max = 0;
		double bwidth = 0;
	
		// For Strings
		string ssplit ="";
		string cats ="";
		string lcats ="";
		string rcats ="";
		int ncats = 0;
		int nlcats = 0;
		int nrcats = 0;

		// General
		Dataline::type ntype;
		string pred = "";
		vector<pair<string, int>>pathpreds;
		int depth = 0;
		int nnode;
		int dim;
		double psplit = 0;
		double prule = 0;
		VectorXd state;
		MatrixXd svar;
		MatrixXd atmat;
		MatrixXd aimat;
		VectorXd cvec;
		VectorXd dvec;
		double ppost = 0;
		double lpost = 0;
		bool isleaf = true;
		bool isleft = false;
		bool isright = false;

		// Constructors
		Node() {};
		Node(int item, int MDIM) : 
			nnode(item), dim(MDIM), state(MDIM), 
			svar(MDIM, MDIM), atmat(MDIM, MDIM), 
			aimat(MDIM, MDIM), cvec(MDIM), dvec(MDIM){};
};
/*
int main() {
	map<int, Node*> tree;
	vector<string> cats = {"a", "b", "c"};
	for (int i = 1; i < 4; ++i) {
		tree[i] = new Node(i, 1, 0, 0, 10);
		tree[i]->dsplit = rand()%10;
		tree[i]->pred = "X1";
	}

	for (int i = 4; i < 7; ++i) {
		tree[i] = new Node(i, 1, 3, 0, 0);
		tree[i]->cats = cats;
		tree[i]->ssplit = 
			tree[i]->cats[rand()%tree[i]->ncats];
		tree[i]->pred = "X2";
	}

	map<int, Node*>::iterator it;
	cout << "My tree contains:\n";
	for (it = tree.begin(); it != tree.end(); ++it) {
		cout << it->second->nnode << 
			" " << it->second->pred <<
			" " << it->second->dsplit << 
			" " << it->second->ssplit << endl;
	}
}
*/
