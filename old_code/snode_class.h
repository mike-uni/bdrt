// A node class

# pragma once
# include <iostream>
# include <vector>
# include <Eigen/Dense>
# include <algorithm>

using namespace std;
using namespace Eigen;

class sNode {
	public:
		string pred = "NONE";
		double split = 0.0;
		int depth = 0;
		int nnode;
		int dim;
		double psplit = 0.0;
		double prule = 0.0;
		VectorXd zcur;
		VectorXd zprev;
		VectorXd ycur;
		MatrixXd wvar;
		MatrixXd vvar;
		bool isleaf = true;
		
		sNode * left = nullptr;
		sNode * right = nullptr;
		sNode * parent = this;
		
		// Constructors
		sNode() {};
		sNode(int item, int DIM) : nnode(item), dim(DIM),
							zcur(DIM), zprev(DIM),
							ycur(DIM), wvar(DIM, DIM),
							vvar(DIM, DIM) {};
		~sNode();
		sNode(const sNode &other); 
		sNode operator = (const sNode &rhs);
};

sNode::sNode(const sNode &other) {
	pred = other.pred;
	split = other.split;
	depth = other.depth;
	nnode = other.nnode;
	psplit = other.psplit;
	prule = other.prule;
	zcur = other.zcur;
	zprev = other.zprev;
	ycur = other.ycur;
	wvar = other.wvar;
	vvar = other.vvar;
	isleaf = other.isleaf;

  	if (other.left != nullptr ) {
        left = new sNode(other.left->nnode, other.left->dim);
    }
    if (other.right != nullptr) {
        right = new sNode(other.right->nnode, other.right->dim);
    }
    if (other.parent != nullptr) {
        parent = new sNode(other.parent->nnode, other.parent->dim);
    }
 
}

sNode sNode:: operator = (const sNode &rhs) {
	pred = rhs.pred;
	split = rhs.split;
	depth = rhs.depth;
	nnode = rhs.nnode;
	psplit = rhs.psplit;
	prule = rhs.prule;
	zcur = rhs.zcur;
	zprev = rhs.zprev;
	ycur = rhs.ycur;
	wvar = rhs.wvar;
	vvar = rhs.vvar;
	isleaf = rhs.isleaf;
	
	sNode * newleft = nullptr;
	sNode * newright = nullptr;
	sNode * newparent = nullptr;
	swap(left, newleft);
	swap(right, newright);
    return *this;
}

sNode::~sNode() {
	delete left;
	delete right;
}
/*
int main() {
	sNode a(1);
	cout << a.left << endl;
	cout << a.right << endl;
	cout << a.parent << endl;
	cout << endl;

	a.prule = 2;
	a.split = 4;

	sNode b = a;
	
	cout << b.parent << endl;
	cout << b.left << endl;
	cout << b.right << endl;
	cout << a.parent << endl;
	
	cout << a.prule << " " << b.prule << " " << a.prule << endl;
}*/
