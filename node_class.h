// A node class

# pragma once
# include <iostream>
# include <vector>
# include <Eigen/Dense>
# include <algorithm>

using namespace std;
using namespace Eigen;

class Node {
	public:
		string pred = "NONE";
		double split = 0.0;
		int depth = 0;
		int nnode;
		int dim;
		double psplit = 0.0;
		double prule = 0.0;
		VectorXd state;
		MatrixXd svar;
		MatrixXd atmat;
		MatrixXd aimat;
		VectorXd cvec;
		VectorXd dvec;
		double ppost;
		double lpost;
		bool isleaf = true;
		
		Node * left = nullptr;
		Node * right = nullptr;
		Node * parent = nullptr;
		
		// Constructors
		Node() {};
		Node(int item, int DIM) : nnode(item), dim(DIM),
							state(DIM), svar(DIM, DIM),
							atmat(DIM, DIM), aimat(DIM, DIM),
							cvec(DIM), dvec(DIM){};
		~Node();
		Node(const Node &other); 
		Node operator = (const Node &rhs);
};

Node::Node(const Node &other) {
	pred = other.pred;
	split = other.split;
	depth = other.depth;
	nnode = other.nnode;
	psplit = other.psplit;
	prule = other.prule;
	state = other.state;
	svar = other.svar;
	atmat = other.atmat;
	aimat = other.aimat;
	cvec = other.cvec;
	dvec = other.dvec;
	ppost = other.ppost;
	lpost = other.lpost;
	isleaf = other.isleaf;

  	if (other.left != nullptr ) {
        left = new Node(other.left->nnode, other.left->dim);
    }
    if (other.right != nullptr) {
        right = new Node(other.right->nnode, other.right->dim);
    }
//    if (other.parent != nullptr) {
//        parent = new Node(other.parent->nnode, other.parent->dim);
//    }
}

Node Node:: operator = (const Node &rhs) {
	pred = rhs.pred;
	split = rhs.split;
	depth = rhs.depth;
	nnode = rhs.nnode;
	psplit = rhs.psplit;
	prule = rhs.prule;
	state = rhs.state;
	svar = rhs.svar;
	aimat = rhs.aimat;
	atmat = rhs.atmat;
	cvec = rhs.cvec;
	dvec = rhs.dvec;
	ppost = rhs.ppost;
	lpost = rhs.lpost;
	isleaf = rhs.isleaf;
	
	Node * newleft = nullptr;
	Node * newright = nullptr;
	//Node * newparent = nullptr;
	swap(left, newleft);
	swap(right, newright);
	//swap(parent, newparent);
    return *this;
}

Node::~Node() {
	//cout << "Node destructor called.\n";
	delete left;
	delete right;
//	if (parent != nullptr) {
//	delete parent;
//	}
}
/*
int main() {
	Node a(1, 1);
	cout << a.left << endl;
	cout << a.right << endl;
	cout << a.parent << endl;
	cout << endl;

	a.prule = 2;
	a.split = 4;

	Node b = a;
	
	cout << b.parent << endl;
	b.prule = 3;
	cout << b.left << endl;
	cout << b.right << endl;
	cout << a.parent << endl;
	
	cout << a.prule << " " << b.prule << " " << a.prule << endl;
}*/
