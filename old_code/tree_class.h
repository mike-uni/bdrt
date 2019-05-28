// A Dynamic Regression Tree class

# pragma once
# include <iostream>
# include <string>
# include <Eigen/Dense>
# include <vector>
# include <algorithm>
# include "read_data.h"
# include "node_class.h"
# include "priorf.h"

using namespace std;
using namespace Eigen;

class Tree {
	public:
		Node * root;
		int tname;
		double prior = 1.0;
		double prop_new = 0.0;
		double prop_current = 0.0;
		double tree_post = 0.0;
		int move = 0;
	
		// Constructors
		Tree() {};	
		Tree(int nnode, int DIMS) : root(new Node(nnode, DIMS)){};
		Tree(const Tree& other);
		~Tree() { delete root;
		//		cout << "Tree destructor called.\n";
		};
		Tree& operator = (const Tree& rhs);

		Node * find_leaf(Node * leaf, const RowVectorXd &vec,
						vector<string> names);
		double get_prior() {return prior;}; 
		void set_prior(double &new_prior) {prior = new_prior;};
		double get_post() {return tree_post;}; 
		void set_post(double &new_post) {tree_post = new_post;};
/*		void grow(Node * leaf, const double ALPHA, const double BETA,
					double in_rule, int dim); 
		void prune(Node * leaf, int dim);
		void change(Node * leaf, const MatrixXd& data, 
					vector<string> namesi, int iter); 
		void nswap(Node * leaf);
*/		void display(Node * leaf);
		int size(Node * leaf) const;
		int nleaves(Node * leaf); 
		Node * find_node(Node * leaf, int node_num);
		void copy_nodes(Node * troot, Node * oroot);
		void assign_nodes(Node * troot, Node * oroot);
//		void destroy_nodes(Node * dnode);
		void find_nnodes(Node * leaf, vector<int>& outvec);
		void find_prior(Node * leaf); 
		void find_leaf_nodes(Node * leaf, vector<int>& leafvec);
		void find_int_nodes(Node * leaf, vector<int>& intvec);
		void find_post(Node * leaf);
};

Tree::Tree(const Tree& other) {
	//delete root;
	root = new Node(other.root->nnode, other.root->dim);
	*root = *other.root;
	prior = other.prior;
	prop_new = other.prop_new;
	prop_current = other.prop_current;
	tree_post = other.tree_post;
	move = other.move;

//	Node * temp(new Node);
//	*temp = *other.root;
//	swap(root, temp);
	copy_nodes(root, other.root);
	//delete temp;
}

Tree& Tree:: operator = (const Tree& rhs) {
	delete root;
	root = new Node(rhs.root->nnode, rhs.root->dim);
	*root = *rhs.root;
	assign_nodes(root, rhs.root);
	prior = rhs.prior;
	prop_new = rhs.prop_new;
	prop_current = rhs.prop_current;
	tree_post = rhs.tree_post;
	move = rhs.move;

	return *this;
}

void Tree:: copy_nodes(Node * troot, Node *oroot) {
//	cout << "copy called\n\n\n";
	if (oroot->left != nullptr) {
		//delete troot->left;
		//troot->left = new Node(oroot->left->nnode, oroot->left->dim);
		Node * temp(new Node);
		*temp = *oroot->left;
		swap(troot->left, temp);
	//	delete temp;
		if (troot->left) {
			troot->left->parent = troot;
		}
		copy_nodes(troot->left, oroot->left);
		//delete troot;	
	}
	if (oroot->right != nullptr) {
		//delete troot->right;
		//troot->right = new Node(oroot->right->nnode, oroot->right->dim);
		Node * temp(new Node);
		*temp = *oroot->right;
		swap(troot->right, temp);
	//	delete temp;
		if (troot->right) {
			troot->right->parent = troot;
		}
		copy_nodes(troot->right, oroot->right);
		//delete troot;
	}
}

void Tree:: assign_nodes(Node * troot, Node * oroot) {
	if (oroot->left != nullptr) {
		Node * temp(new Node);
		*temp = *oroot->left;
		swap(troot->left, temp);
//		delete temp;
		//troot->left->parent = troot->parent;
		assign_nodes(troot->left, oroot->left);
	}
	if (oroot->right != nullptr) {
		Node * temp(new Node);
		*temp = *oroot->right;
		swap(troot->right, temp);
//		delete temp;
		//troot->right->parent = troot->parent;
		assign_nodes(troot->right, oroot->right);
	}
}

/*
void Tree:: assign_nodes(Node * troot, Node * oroot) {
	if (oroot->left != nullptr) {
		delete troot->left;
		troot->left = new Node(oroot->left->nnode, oroot->left->dim);
		*troot->left = *oroot->left;
		troot->left->parent = troot;
		assign_nodes(troot->left, oroot->left);
	}
	if (oroot->right != nullptr) {
		delete troot->right;
		troot->right = new Node(oroot->right->nnode, oroot->right->dim);
		*troot->right = *oroot->right;
		troot->right->parent = troot;
		assign_nodes(troot->right, oroot->right);
	}
}

void Tree:: destroy_nodes(Node * dnode) {
	if (dnode->parent == nullptr) {
		destroy_nodes(dnode->left);
		destroy_nodes(dnode->right);
		delete dnode;
	}
}
*/
Node * Tree:: find_leaf(Node * leaf, const RowVectorXd &vec,
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

Node * Tree:: find_node(Node * leaf, int node_num) {
	static Node * temp;
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

void Tree::find_leaf_nodes(Node * leaf, vector<int>& leafvec) {
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

void Tree::find_int_nodes(Node * leaf, vector<int>& intvec) {
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
/*
void Tree::find_post(Node * leaf) {
	vector<int> lnodes;
	find_leaf_nodes(leaf, lnodes);
	for (int i = 0; i < lnodes.size(); i++) {
		find_node
}
*/
void Tree:: find_post(Node * leaf) {
	double mypost = get_post(); 
	if (leaf != NULL) {
		if (leaf->isleaf == true) {
			mypost = mypost + leaf->lpost;
			set_post(mypost);
			cout << "isleaf is true and tree_post is: " << mypost << endl;
		}
		else {
			mypost = mypost;
			set_post(mypost);
			cout << "isleaf is false and tree_post is: " << mypost << endl;
		}
		find_post(leaf->left);
		find_post(leaf->right);
	}
}

void Tree:: find_nnodes(Node * leaf, vector<int>& outvec) {
	if (leaf != NULL) {
		outvec.push_back(leaf->nnode);
		find_nnodes(leaf->left, outvec);
		find_nnodes(leaf->right, outvec);
	}
}

void Tree:: find_prior(Node * leaf) {
	double myprior = get_prior(); 
	if (leaf != NULL) {
		if (leaf->isleaf == true) {
			myprior = myprior*(1-(leaf->psplit));
			set_prior(myprior);
		}
		else {
			myprior = myprior*(leaf->psplit)*(leaf->prule);
			set_prior(myprior);
		}
		find_prior(leaf->left);
		find_prior(leaf->right);
	}
}

int Tree:: size(Node * leaf) const {
    if(!leaf) 
        return 0;
    else
        return size(leaf->left) + 1 + size(leaf->right); 
}

int Tree:: nleaves(Node * leaf){
 	if (leaf == NULL)
    	return 0;
	if (leaf->left == NULL && leaf->right == NULL) {
    	return 1;
  	} 
	else {
    	return nleaves(leaf->left) + nleaves(leaf->right);
  	}
}

void Tree:: display(Node * leaf) {
	if (leaf != NULL) {
		cout << "The leaf nnode is: " << leaf->nnode << endl; 
		cout << "isleaf is: " << leaf->isleaf << endl; 
		//				<< " at leaf address " << leaf << endl;
		//cout << "The split is: " << leaf->split << endl;
		//cout << "The pred is: " << leaf->pred << endl;
		//cout << "The psplit is: " << leaf->psplit << endl;
		//if (leaf->parent == nullptr) {
		//	cout << "leaf parent is nullptr\n";
		//}
		//else {
		//	cout << "Parent nnode is " << leaf->parent->nnode << endl;
		//}
		//cout << "The parent address is: " << leaf->parent << endl;
		//				" at address" << leaf->parent << endl;
		cout << "lpost is: " << leaf->lpost << endl;
		cout << "ppost is: " << leaf->ppost << endl;

		display(leaf->left);
		display(leaf->right);
	}
}
/*
void Tree:: grow(Node * leaf, const double alph, const double bet,
					double in_rule, int dim) {
	leaf->isleaf = false;

	Node * templ, * tempr;
	templ = leaf->left;
	tempr = leaf->right;

	templ = new Node(2*(leaf->nnode), dim);
	templ->depth = leaf->depth+1;
	templ->psplit = psplitf(alph, bet, templ->depth);
	templ->prule = in_rule;
	templ->parent = leaf;
	templ->state = leaf->state;
	templ->svar = leaf->svar;
	leaf->left = templ;

	tempr = new Node(2*(leaf->nnode)+1, dim);
	tempr->depth = leaf->depth+1;
	tempr->psplit = psplitf(alph, bet, tempr->depth);
	tempr->prule = in_rule;
	tempr->parent = leaf;
	tempr->state = leaf->state;
	tempr->svar = leaf->svar;
	leaf->right = tempr;

	//delete templ;
	//delete tempr;
}

void Tree:: prune(Node * leaf, int dim) {
	VectorXd left_state(dim);
	left_state = leaf->left->state;
	MatrixXd left_svar(dim, dim);
	left_svar = leaf->left->svar;
	VectorXd right_state(dim);
	right_state = leaf->right->state;
	MatrixXd right_svar(dim, dim);
	right_svar = leaf->right->svar;
	delete leaf->left;
	leaf->left = nullptr;
	delete leaf->right;
	leaf->right = nullptr;
//	leaf->parent = nullptr;
	leaf->isleaf = true;
//	leaf->state = 0.5*(left_state + right_state);
//	leaf->svar = 0.5*(left_svar + right_svar);
}

void Tree:: change(Node * leaf, const MatrixXd& data, 
		vector<string> names, int iter) {
	int dlength = data.cols();
	int dsplits = data.rows();
	pair<string, double> newpair;
	newpair = rand_pair(data, names, dlength, dsplits, iter);
	leaf->pred = newpair.first;
	leaf->split = newpair.second; 
}

void Tree:: nswap(Node * leaf) {
	if (leaf->nnode == 1) {
		double temp_split = leaf->left->split;
		string temp_pred = leaf->left->pred;
		leaf->left->split = leaf->right->split;
		leaf->left->pred = leaf->right->pred;
		leaf->right->split = temp_split;
		leaf->right->pred = temp_pred; 
	}	
	else {
		double temp_split = leaf->split;
		string temp_pred = leaf->pred;
		leaf->split = leaf->parent->split;
		leaf->pred = leaf->parent->pred;
		leaf->parent->split = temp_split;
		leaf->parent->pred = temp_pred;
	}
}

int main() {
	static int const DIM = 2;
	Tree b(1, DIM);
	
	ifstream file("xsim_data2.txt", ios::in);
	MatrixXd data(5, 10);
	VectorXd dataline(10);
	for (int i = 0; i < 5; i++) {
		dataline = read_line(i*10, 10*(i+1), file);
		data.row(i) = dataline;
	}

	VectorXd dataLine(10);
	dataLine = data.row(1);
	vector<string> myNames = make_names(100, 10);
	pair<string, double> nrandp = rand_pair(data, myNames, 10, 5, 1);
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
	b.root->psplit = psplitf(0.95, 1, b.root->depth);
	b.root->prule = 0.2;
	b.prior = 10;
	b.display(b.root);
	VectorXd yi(DIM);
	yi = sim_y(data.row(1), tsig, Amat, HT);
	//cout << yi << endl;
	cout << b.prior <<  endl;

	Tree c = b;
	c.display(c.root);
	cout << endl;
	c.root->pred = "X5";
	c.root->split = 34;
	c.prior = 30;
	cout << c.root << endl;
	cout << b.root << endl;
	b.grow(b.root, 0.95, 1, 0.2, DIM);

	cout << endl;
	cout << "c display is:\n";
	c.display(c.root);
	cout << c.prior << endl << endl;
	cout << "b display is:\n";
	b.display(b.root);
	cout << b.prior << endl << endl;

	cout << "Assign c to b:\n";
	b = c;
	b.display(b.root);

	b.grow(b.root, 0.95, 1, 0.2, DIM);
	cout << "b display is:\n";
	b.display(b.root);
	cout << endl; 

	VectorXd dataLine1(10);
	dataLine1 = data.row(2);
	//cout << dataLine1 << endl;

	Node * found;
	found = b.find_leaf(b.root, dataLine1, myNames);
	cout << "Found address is: " << found << endl;
	cout << endl;

	pair<string, double> nrandp1 = rand_pair(data, myNames, 10, 5, 2);
	found->pred = nrandp1.first;
	found->split = nrandp1.second;
	b.display(b.root);
	cout << endl;
	cout << "Node depth is: " << found->depth << endl;
	cout << "The address of leaf is: " << found << endl;
	cout << "The address of found->parent is: " << found->parent 
				<< endl;
	cout << "The address of parent node of root is: " << 
				b.root->parent << endl;
	cout << "The found->parent->nnode is: " << found->parent->nnode
				<< endl;
	cout << endl;

	Node * found_int;
	found_int = found->parent; 
	b.nswap(found_int);
	b.display(b.root);
	cout << endl; 
	
	b.change(found, dataLine1, myNames, 1);
	b.display(b.root);
	cout << endl;
	
	cout << "Number of leaves is: " << b.nleaves(b.root) << endl;

	vector<int> outvec;
 	b.find_nnodes(b.root, outvec);
	RowVectorXi nodevec(b.size(b.root));
	for (int j = 0; j < nodevec.size(); j++) {
		nodevec(j) = outvec[j];
	}
	cout << "Node numbers are: " << nodevec << endl;

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

	VectorXd dataLine2(10);
	dataLine2 = data.row(3);

	b.grow(found, 0.95, 1, 0.2, DIM);
	b.display(b.root);
	cout << endl; 

	Node * found2 = b.find_leaf(b.root, dataLine2, myNames);
	cout << "Found2 nnode is: " << found2->nnode << endl;
	pair<string, double> nrandp2 = rand_pair(data, myNames, 10, 5, 3);
	found2->pred = nrandp2.first;
	found2->split = nrandp2.second;
	b.display(b.root);
	cout << endl;
	
	b.grow(found2, 0.95, 1, 0.2, DIM);
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
	
	Node * found3;
	found3 = b.find_node(b.root, 7);
	cout << found3->nnode << endl;	

	Node * found4;
	found4 = b.find_node(b.root, 6);
	cout << found4->nnode << endl;

	b.grow(found4, 0.95, 1, 0.2, DIM);
	b.display(b.root);
	Node * found5;
	found5 = b.find_node(b.root, 12);
	cout << found5->nnode << endl;

	VectorXd dataLine6(10);
	dataLine6 = data.row(4);
	Node * found6 = b.find_leaf(b.root, dataLine6, myNames);
	cout << found6->nnode << endl;

	Tree d = b;	
	d.display(d.root);
	cout << endl;

	b.find_prior(b.root);
	cout << "The prior is: " << b.prior << endl;
	cout << endl;

	//cout << b.root->left->parent << endl;
	//cout << b.root->right->parent << endl;
	//cout << b.root->right->right->parent << endl;
	//cout << b.root->right->left->parent << endl;

	Tree e = b;
	cout << "Tree b display is: " << endl;
	b.display(b.root);
	cout << endl;
	
	cout << "Tree e display is: " << endl;
	e.display(e.root);
	cout << endl;
	

	Node * founde6 = e.find_leaf(e.root, dataLine6, myNames);
	cout << "founde6 nnode is: " << founde6->nnode << endl;
	cout << "found6 nnode is: " << found6->nnode << endl;
	Node * founde2 = e.find_leaf(e.root, dataLine2, myNames);
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
	e.prune(founde6->parent, DIM);
	e.display(e.root);
	cout << endl; 
	//cout << "Found6 nnode is: " << found6->nnode << endl;
	Node * founde2 = e.find_leaf(e.root, dataLine2, myNames);
	e.grow(founde2, 0.95, 1, 0.2, DIM);
	//cout << "Found6 Parent is: " << found6->parent->nnode << endl;
	//b = e;
	}
	cout << endl;
}*/
