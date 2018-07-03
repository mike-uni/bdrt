// A Dynamic Regression sTree class

# pragma once
# include <iostream>
# include <string>
# include <Eigen/Dense>
# include <vector>
# include <algorithm>
# include "read_data.h"
# include "snode_class.h"
# include "priorf.h"

using namespace std;
using namespace Eigen;

class sTree {
	public:
		sNode * root;
	
		// Constructors
		sTree() {};	
		sTree(int nnode, int DIMS) : 
				root(new sNode(nnode, DIMS)) {};
		sTree(const sTree& other);
		~sTree() {};
		sTree& operator = (const sTree& rhs);

		sNode * find_leaf(sNode * leaf, const RowVectorXd &vec,
						vector<string> names);
		void display(sNode * leaf);
		int size(sNode * leaf) const;
		int nleaves(sNode * leaf); 
		sNode * find_node(sNode * leaf, int node_num);
		void copy_nodes(sNode * oroot, sNode *troot);
		void assign_nodes(sNode * oroot, sNode *troot);
		void find_nnodes(sNode * leaf, vector<int>& outvec);
		void find_prior(sNode * leaf); 
		void find_leaf_nodes(sNode * leaf, vector<int>& leafvec);
		void find_int_nodes(sNode * leaf, vector<int>& intvec);
		void find_post(sNode * leaf);
		void sgrow(sNode * leaf, const int DIM,  string lpred, 
					double lsplit, string rpred, double rsplit);
};

sTree::sTree(const sTree& other) {
	root = new sNode(other.root->nnode, other.root->dim);
	*root = *other.root;
	copy_nodes(other.root, root);
}

sTree& sTree:: operator = (const sTree& rhs) {
	delete root;
	root = new sNode(rhs.root->nnode, rhs.root->dim);
	*root = *rhs.root;
	assign_nodes(rhs.root, root);	
	return *this;
}

void sTree:: copy_nodes(sNode * oroot, sNode *troot) {
	if (oroot->left != nullptr) {
		troot->left = new sNode(oroot->left->nnode, oroot->left->dim);
		*troot->left = *oroot->left;
		troot->left->parent = troot;
		copy_nodes(oroot->left, troot->left);	
	}
	if (oroot->right != nullptr) {
		troot->right = new sNode(oroot->right->nnode, oroot->right->dim);
		*troot->right = *oroot->right;
		troot->right->parent = troot;
		copy_nodes(oroot->right, troot->right);
	}
}

void sTree:: assign_nodes(sNode * oroot, sNode *troot) {
	if (oroot->left != nullptr) {
		delete troot->left;
		troot->left = new sNode(oroot->left->nnode, oroot->left->dim);
		*troot->left = *oroot->left;
		troot->left->parent = troot;
		copy_nodes(oroot->left, troot->left);
	}
	if (oroot->right != nullptr) {
		delete troot->right;
		troot->right = new sNode(oroot->right->nnode, oroot->right->dim);
		*troot->right = *oroot->right;
		troot->right->parent = troot;
		copy_nodes(oroot->right, troot->right);
	}
}

sNode * sTree:: find_leaf(sNode * leaf, const RowVectorXd &vec,
						vector<string> names) {
	string npred = leaf->pred;
	double nsplit = leaf->split;
	double temp_split;
	int index = get_pred_index(npred, names);
	temp_split = vec(index);
	if (leaf!= NULL) {
		if (leaf->isleaf) {
			return leaf;
		}
		else if (temp_split < nsplit) {
			return find_leaf(leaf->left, vec, names);
		}
		else {
			return find_leaf(leaf->right, vec, names);
		}
	}	
	else {
		return NULL;
	}
}

sNode * sTree:: find_node(sNode * leaf, int node_num) {
	static sNode * temp;
 	if (leaf != NULL) {
		if (leaf->nnode == node_num) {
    		temp = leaf;
  		} 
		else {
    		find_node(leaf->left, node_num);
			find_node(leaf->right, node_num);
		}
  	}
	return temp;
}

void sTree::find_leaf_nodes(sNode * leaf, vector<int>& leafvec) {
 	if (leaf != NULL) {
		if (leaf->isleaf) {
    		return leafvec.push_back(leaf->nnode);
  		} 
		else {
    		find_leaf_nodes(leaf->left, leafvec);
			find_leaf_nodes(leaf->right, leafvec);
		}
  	}
}

void sTree::find_int_nodes(sNode * leaf, vector<int>& intvec) {
 	if (leaf != NULL) {
		if (!leaf->isleaf) {
    		return intvec.push_back(leaf->nnode);
  		} 
		else {
    		find_int_nodes(leaf->left, intvec);
			find_int_nodes(leaf->right, intvec);
		}
  	}
}

void sTree:: find_nnodes(sNode * leaf, vector<int>& outvec) {
	if (leaf != NULL) {
		outvec.push_back(leaf->nnode);
		find_nnodes(leaf->left, outvec);
		find_nnodes(leaf->right, outvec);
	}
}

int sTree:: size(sNode * leaf) const {
    if(!leaf) 
        return 0;
    else
        return size(leaf->left) + 1 + size(leaf->right); 
}

int sTree:: nleaves(sNode * leaf){
 	if (leaf == NULL)
    	return 0;
	if (leaf->left == NULL && leaf->right == NULL) {
    	return 1;
  	} 
	else {
    	return nleaves(leaf->left) + nleaves(leaf->right);
  	}
}

void sTree:: display(sNode * leaf) {
	if (leaf != NULL) {
		cout << "The nnode is: " << leaf->nnode << endl; 
		//				<< " at leaf address " << leaf << endl;
		//cout << "The zprev is:\n" << leaf->zprev << endl; 
		//cout << "The wvar is:\n" << leaf->wvar << endl; 
		cout << "The split is: " << leaf->split << endl;
		cout << "The pred is: " << leaf->pred << endl << endl;
	//	cout << "The parent leaf nnode is: " << leaf->parent->nnode<<endl; 
		//cout << "The psplit is: " << leaf->psplit << endl;
		//cout << "Lpost is: " <<  leaf->lpost << endl;
		//cout << "The parent address is: " << leaf->parent << endl;
		//				" at address" << leaf->parent << endl;

		display(leaf->left);
		display(leaf->right);
	}
}

void sTree:: sgrow(sNode * leaf, const int DIM, string lpred, double lsplit,
					string rpred, double rsplit) {

	leaf->isleaf = false;
		
	sNode * templ, * tempr;
	templ = leaf->left;
	tempr = leaf->right;

	templ = new sNode(2*leaf->nnode, DIM);
	templ->pred = lpred;
	templ->split = lsplit;
	templ->dim = leaf->dim;
	templ->depth = leaf->depth+1;
	templ->parent = leaf;
	templ->zcur = VectorXd::Zero(DIM);
	templ->zprev = VectorXd::Zero(DIM);
	templ->ycur = VectorXd::Zero(DIM);
	leaf->left = templ;

	tempr = new sNode(2*leaf->nnode+1, DIM);
	tempr->pred = rpred;
	tempr->split = rsplit;
	tempr->dim = leaf->dim;
	tempr->depth = leaf->depth+1;
	tempr->parent = leaf;
	tempr->zcur = VectorXd::Zero(DIM);
	tempr->zprev = VectorXd::Zero(DIM);
	tempr->ycur = VectorXd::Zero(DIM);
	leaf->right = tempr;
	
}

/*
int main() {
	int const DIM = 2;
	sTree b(1, 1);
	
	ifstream file("xsim_data.txt", ios::in);
	VectorXd dataLine = read_line(0, 10, file);
	vector<string> myNames = make_names(100, 10);
	pair<string, double> nrandp = rand_pair(dataLine, myNames);
	MatrixXd Amat(DIM, 10);
	for (int i = 0; i < Amat.size(); i++) {
		Amat(i) = rand()%2;
	}
	MatrixXd tsig(DIM, DIM);
	tsig << 1,0,0,1;
	MatrixXd HT(DIM, DIM);
	HT << 1,0,0,1; 
	//cout << Amat << endl;
	//cout << tsig << endl;
	//cout << Amat*dataLine << endl;
	cout << endl;
	
	b.root->pred = nrandp.first;
	b.root->split = nrandp.second;
	b.root->psplit = psplit(0.95, 1, b.root->depth);
	b.root->prule = 0.2;
	b.prior = 10;
	b.display(b.root);
	VectorXd yi(DIM);
	yi = sim_y(dataLine, tsig, Amat, HT);
	//cout << yi << endl;
	cout << b.prior <<  endl;

	sTree c = b;
	c.display(c.root);
	cout << endl;
	c.root->pred = "X5";
	c.root->split = 34;
	c.prior = 30;
	cout << c.root << endl;
	cout << b.root << endl;
	b.grow(b.root, 0.95, 1, 0.2);

	cout << endl;
	c.display(c.root);
	cout << c.prior << endl;
	b.display(b.root);
	cout << b.prior << endl;

	b = c;
	b.display(b.root);

	b.grow(b.root, 0.95, 1, 0.2);
	b.display(b.root);
	cout << endl; 

	VectorXd dataLine1;
	dataLine1 = read_line(10, 20, file);
	cout << dataLine1 << endl;

	sNode * found;
	found = b.find_leaf(b.root, dataLine1, myNames);
	cout << "Found address is: " << found << endl;
	cout << endl;

	pair<string, double> nrandp1 = rand_pair(dataLine1, myNames);
	found->pred = nrandp1.first;
	found->split = nrandp1.second;
	b.display(b.root);
	cout << endl;
	cout << "sNode depth is: " << found->depth << endl;
	cout << "The address of leaf is: " << found << endl;
	cout << "The address of found->parent is: " << found->parent 
				<< endl;
	cout << "The address of parent node of root is: " << 
				b.root->parent << endl;
	cout << "The found->parent->nnode is: " << found->parent->nnode
				<< endl;
	cout << endl;

	sNode * found_int;
	found_int = found->parent; 
	b.swap(found_int);
	b.display(b.root);
	cout << endl; 
	
	b.change(found, dataLine1, myNames);
	b.display(b.root);
	cout << endl;
	
	cout << "Number of leaves is: " << b.nleaves(b.root) << endl;

	vector<int> outvec;
 	b.find_nnodes(b.root, outvec);
	RowVectorXi nodevec(b.size(b.root));
	for (int j = 0; j < nodevec.size(); j++) {
		nodevec(j) = outvec[j];
	}
	cout << "sNode numbers are: " << nodevec << endl;

	vector<int> leaf_nodes;
 	b.find_leaf_nodes(b.root, leaf_nodes);
	RowVectorXi leafvec(b.nleaves(b.root));
	for (int j = 0; j < leafvec.size(); j++) {
		leafvec(j) = leaf_nodes[j];
	}
	cout << "Leaf node numbers are: " << leafvec << endl;

	cout << "Leaf 2 address is: " << b.root->left << endl;
	cout << "Leaf 3 address is: " << b.root->right << endl;
	
	cout << b.find_node(b.root, 2) << endl;
	cout << b.find_node(b.root, 3) << endl;

	VectorXd dataLine2 = read_line(20, 30, file);

	b.grow(found, 0.95, 1, 0.2);
	b.display(b.root);
	cout << endl; 

	sNode * found2 = b.find_leaf(b.root, dataLine2, myNames);
	cout << "Found2 nnode is: " << found2->nnode << endl;
	pair<string, double> nrandp2 = rand_pair(dataLine2, myNames);
	found2->pred = nrandp2.first;
	found2->split = nrandp2.second;
	b.display(b.root);
	cout << endl;
	
	b.grow(found2, 0.95, 1, 0.2);
	b.display(b.root);

	int	num_leaves = b.nleaves(b.root);
	cout << "num_leaves are: " << num_leaves << endl;

	vector<int> leaf_nodes1;
 	b.find_leaf_nodes(b.root, leaf_nodes1);
	cout << "The leaf nodes1 are: " << endl;
	for (int j = 0; j < num_leaves; ++j) {
		cout << leaf_nodes1[j] << " ";
	}
	cout << endl;
	
	sNode * found3;
	found3 = b.find_node(b.root, 7);
	cout << found3->nnode << endl;	

	sNode * found4;
	found4 = b.find_node(b.root, 6);
	cout << found4->nnode << endl;

	b.grow(found4, 0.95, 1, 0.2);
	b.display(b.root);
	sNode * found5;
	found5 = b.find_node(b.root, 12);
	cout << found5->nnode << endl;

	VectorXd dataLine6 = read_line(60, 70, file);
	sNode * found6 = b.find_leaf(b.root, dataLine6, myNames);
	cout << found6->nnode << endl;

	sTree d = b;	
	d.display(d.root);
	cout << endl;

	b.find_prior(b.root);
	cout << "The prior is: " << b.prior << endl;
	cout << endl;

	//cout << b.root->left->parent << endl;
	//cout << b.root->right->parent << endl;
	//cout << b.root->right->right->parent << endl;
	//cout << b.root->right->left->parent << endl;

	sTree e = b;
	cout << "sTree b display is: " << endl;
	b.display(b.root);
	cout << endl;
	
	cout << "sTree e display is: " << endl;
	e.display(e.root);
	cout << endl;
	

	sNode * founde6 = e.find_leaf(e.root, dataLine6, myNames);
	cout << "founde6 nnode is: " << founde6->nnode << endl;
	cout << "found6 nnode is: " << found6->nnode << endl;
	sNode * founde2 = e.find_leaf(e.root, dataLine2, myNames);
	cout << "founde2 nnode is: " << founde2->nnode << endl;
	cout << "found2 nnode is: " << found2->nnode << endl;
	cout << endl;

	//e.prune(founde6->parent);
	//e.display(e.root);	
	//cout << endl;

	for (int i = 0; i < 5; i++) {
	e.display(e.root);
	cout << endl; 
	//cout << "Found2 nnode is: " << found2->nnode << endl;
	cout << "e parent is: " << founde6->parent << endl;
	cout << "e parent nnode is: " << founde6->parent->nnode << endl;
	cout << "b parent is: " << found6->parent << endl;
	cout << "b parent nnode is: " << found6->parent->nnode << endl;
	e.prune(founde6->parent);
	e.display(e.root);
	cout << endl; 
	//cout << "Found6 nnode is: " << found6->nnode << endl;
	sNode * founde2 = e.find_leaf(e.root, dataLine2, myNames);
	e.grow(founde2, 0.95, 1, 0.2);
	//cout << "Found6 Parent is: " << found6->parent->nnode << endl;
	//b = e;
	}
	cout << endl;
}*/	
