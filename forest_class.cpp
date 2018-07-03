// A FOREST CLASS

# include <iostream>
# include <Eigen/Dense>
# include <fstream>
# include <algorithm>
# include <vector>
# include <time.h>
# include "read_data.h"
# include "kalman_class.h"
# include "priorf.h"
# include "posterior_class.h"
# include "node_class.h"
# include "tree_class.h"
# include "move_class.h"
# include "proposal_class.h"
# include "dist_class.h"

using namespace std;
using namespace Eigen;

class Forest {
	public:	
	Tree * fptr;
	vector<Tree> fvec;
	
	// Constructors
	Forest() {};
	//Forest(const int FOREST_SIZE){} : fvec(FOREST_SIZE) { 
		//*cout << "Forest constructor called.\n";*/};
	~Forest() {/*cout << "Forest destructor called.\n";*/ };	

	// Declarations
	void init_tree(Tree * tptr, Ref<MatrixXd> ATMAT, 
			Ref<MatrixXd> AIMAT, double ALPHA, double BETA,  
			Ref<VectorXd> CVEC, Ref<VectorXd> DVEC, double& ppost, 
			double& lpost, 
			const vector<string> names, const MatrixXd& data, 
			const VectorXd& MU, const MatrixXd& SIG);
	vector<int> lnodes(Tree * tptr);
	Node * ch_node(Tree * tptr, const VectorXd& data, 
					const vector<string> names);
	void update(Tree * tptr, Node * uptr, /*const VectorXd& data,
					const vector<string> names,*/ vector<int> nodes, 
					const int DIM, const VectorXd& yi,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, double Vdet, double Wdet);
	void proposal(Tree * tptr, Node * fptr,
					const vector<int> leaf_nodes,
					const vector<int> int_nodes,
					const double ALPH, const double BET,
					double in_rule, const int DIM,
					const vector<string> names,
					const MatrixXd& data, int iter,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, const VectorXd& yi,
					double Vdet, double Wdet);
	int ch_tree(Tree * pcur, Tree * pnew);
	void treeout(Node * nptr, ofstream& outfile);
};

void Forest::init_tree(Tree * tptr, Ref<MatrixXd> ATMAT, 
			Ref<MatrixXd> AIMAT, double ALPHA, double BETA, 
			Ref<VectorXd> CVEC, Ref<VectorXd> DVEC, double& ppost,
			double& lpost, 
			const vector<string> names, const MatrixXd& data, 
			const VectorXd& MU, const MatrixXd& SIG) {
	pair<string, double> randpair;
	randpair = rand_pair(data, names, data.cols(), data.rows(), 0);
	tptr->root->pred = randpair.first;
	tptr->root->split = randpair.second;
	tptr->root->psplit = psplitf(ALPHA, BETA, tptr->root->depth);
	tptr->root->prule = (1/(float)data.cols());
	tptr->root->state = MU;
	tptr->root->svar = SIG;
	tptr->root->atmat = ATMAT;
	tptr->root->aimat = AIMAT;
	tptr->root->cvec = CVEC;
	tptr->root->dvec = DVEC;
	tptr->root->ppost = ppost;
	tptr->root->lpost = lpost;
}

vector<int> Forest::lnodes(Tree * tptr) {
	vector<int> outvec;
	tptr->find_leaf_nodes(tptr->root, outvec);
	return outvec;
}
	
Node * Forest::ch_node(Tree * tptr, const VectorXd& data, 
						const vector<string> names) {
	Node * temp;
	temp = tptr->find_leaf(tptr->root, data, names);
	return temp;
}

void Forest::update(Tree * tptr, Node * uptr, /*const VectorXd& data,
					const vector<string> names,*/ vector<int> nodes, 
					const int DIM, const VectorXd& yi,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, double Vdet, double Wdet) {
	int const num_leaves = nodes.size();
	//int const data_size = data.size();

	Kalman k(DIM);
	Lpost l;
	if (num_leaves > 1) {
		int rem_leaf = uptr->nnode;
		nodes.erase(remove(nodes.begin(), nodes.end(), rem_leaf),
						nodes.end());
		Node * temp;
		for (int j = 0; j < (num_leaves-1); j++) {
			temp = tptr->find_node(tptr->root, nodes[j]);
			k.nzhat(temp->state, FMAT);
			k.nsighat(temp->svar, WMAT);
	//		cout << "temp->ppost is: " << temp->ppost
	//			<< "and temp->lpost is: " << temp->lpost << endl;
			l.lpostnup(temp->atmat, temp->aimat, 
				temp->cvec, temp->dvec, 
				temp->ppost, temp->lpost, FWF, WI,
				WF, WFT, Wdet);
		}
	}
	// Predict
	k.zvec(k.z_vec, FMAT, uptr->state);
	k.rmat(k.r_mat, FMAT, WMAT, uptr->svar);
	// Update
	k.yvec(k.y_vec, yi, HMAT, k.z_vec);
	k.smat(k.s_mat, HMAT, k.r_mat, VMAT);
	k.kalg(k.kal_g, k.r_mat, FMAT, k.s_mat);
	k.zhat(uptr->state, k.z_vec, k.kal_g, k.y_vec);
	k.sighat(uptr->svar, k.r_mat, HMAT, k.kal_g);

	l.lpostup(uptr->atmat, uptr->aimat, uptr->cvec, uptr->dvec,
		uptr->ppost, uptr->lpost, yi, VI, FWF, WI, WF,
		WFT, HVH, VH, Vdet, Wdet);

	//cout << "uptr->ppost is: " << uptr->ppost
	//	<< " and uptr->lpost is: " << uptr->lpost << endl;

	tptr->find_post(tptr->root);
}

void Forest::proposal(Tree * tptr, Node * fptr,
					const vector<int> leaf_nodes,
					const vector<int> int_nodes,
					const double ALPH, const double BET,
					double in_rule, const int DIM,
					const vector<string> names,
					const MatrixXd& data, int iter,
					const MatrixXd& WMAT, const MatrixXd& VMAT,
					const MatrixXd& WI, const MatrixXd& VI,
					const MatrixXd& FMAT, const MatrixXd& HMAT,
					const MatrixXd& HVH, const MatrixXd& VH,
					const MatrixXd& WFT, const MatrixXd& WF,
					const MatrixXd& FWF, const VectorXd& yi,
					double Vdet, double Wdet) {
	int num_leaves = leaf_nodes.size();

	Prop p;
	Dist d;	
	int move;
	float pmove;
	if (num_leaves == 1) {
		move = d.runif(1, 2);
		pmove = 1/(float)2;
	}
	else {
		move = d.runif(1, 4);
		pmove = 1/(float)4;
	}
	tptr->move = move; 
	switch (move) {
		case 1: {
		//cout << "The move is grow.\n";
		p.gprop(tptr, fptr, pmove, leaf_nodes, ALPH, BET, in_rule,
					DIM, names, data, iter, WMAT, VMAT, WI, VI,
					FMAT, HMAT, HVH, VH, WFT, WF, FWF, yi,
					Vdet, Wdet);
	/*	cout << "log Tree Posterior is: " << tptr->tree_post << endl;
		cout << "log Tree Prior is: " << log(tptr->prior) << endl;
		cout << "log Tree prop_new is: " << log(tptr->prop_new) <<endl;
		cout << "log Tree prop_current is: "<< log(tptr->prop_current) 
					<< endl << endl; 
	*/	break;
		}
		case 2: {
		//cout << "The move is change.\n";
		p.cprop(tptr, fptr, pmove, int_nodes, iter, data, names);
	/*	cout << "log Tree Posterior is: " << tptr->tree_post << endl;
		cout << "log Tree Prior is: " << log(tptr->prior) << endl;
		cout << "log Tree prop_new is: " << log(tptr->prop_new) <<endl;
		cout << "log Tree prop_current is: "<< log(tptr->prop_current) 
					<< endl << endl; 
	*/	break;
		}
		case 3: {
		//cout << "The move is swap.\n";
		p.sprop(tptr, fptr, pmove, int_nodes, iter);
		break;
		}
		case 4: {
		//cout << "The move is prune.\n";
		p.pprop(tptr, fptr, pmove, int_nodes, ALPH, BET, DIM, names,
					data, iter, WMAT, VMAT, WI, VI, FMAT, HMAT,
					HVH, VH, WFT, WF, FWF, yi, Vdet, Wdet);
		break;
		}
	};
}

int Forest::ch_tree(Tree * pcur, Tree * pnew) {
	Dist d;
	double qcur = log(pcur->prior)+pcur->tree_post;
//	cout << "Current prior is: " << pcur->prior << 
//	" and Current tree posterior is: " << pcur->tree_post << endl; 
	double qnew = log(pnew->prior)+pnew->tree_post;
//	cout << "New prior is: " << pnew->prior << 
//	" and New tree posterior is: " << pnew->tree_post << endl; 
	double prop_ratio = qnew+log(pnew->prop_new)
						-(qcur+log(pnew->prop_current));
/*	if (isnan(lprop_ratio)) {
		cout << "The lprop_ratio is nan." << endl;
		return 0;
	}
	cout << endl;
	if (isinf(lprop_ratio)) {
		cout << "The lprop_ratio is inf." << endl;
		return 0;
	}
*/	double u_init = log(d.runif(0, RAND_MAX)/(double)RAND_MAX);
//	 cout << "log Prop_ratio is: " << prop_ratio << endl;
//	 cout << "log u_init is: " << u_init << endl;
	if (u_init > prop_ratio) {
		//cout << "The move was not accepted." << endl << endl;
		return 1;
	}
	else {
		//cout << "The move was accepted." << endl << endl;
		return 2;
	}
}

void Forest:: treeout(Node * nptr, ofstream& outfile)  {
	if (nptr != nullptr) {
		outfile << nptr->nnode << " " << nptr->pred << " " << 
				nptr->split << " ";
		treeout(nptr->left, outfile);
		treeout(nptr->right, outfile);
	}
}

		
/*	
int main() {
	const int LOOPS = 500;
	const int FOREST = 20;
	const int DIMS = 2;
	const int DATA_LENGTH = 10;
	const int NSPLITS = 20;
	const double ALPH = 0.95;
	const double BET = 1.0;

	MatrixXd W0(DIMS, DIMS);
	W0 << 0.9,0,0,0.9;
	MatrixXd WT(DIMS, DIMS);
	WT << 1,0,0,1;
	MatrixXd VT(DIMS, DIMS);
	VT << 1,0,0,1;
	MatrixXd HT(DIMS, DIMS);
	HT << 1,0,0,1;
	MatrixXd FT(DIMS, DIMS);
	FT << 1,0,0,1;
	MatrixXd SIG(DIMS, DIMS);
	SIG << 1,0,0,1;
	VectorXd MU(DIMS);
	MU << 3,5;
	MatrixXd A_mat(DIMS, DATA_LENGTH);
		for (int i = 0; i < A_mat.size(); i++) {
    	A_mat(i) = rand()%2;
	}

	MatrixXd VTI(DIMS, DIMS);
	VTI = VT.inverse();
	MatrixXd WTI(DIMS, DIMS);
	WTI = WT.inverse();

	MatrixXd FWFT(DIMS, DIMS);
	FWFT = FT.transpose()*WT.inverse()*FT;
	MatrixXd HVHT(DIMS, DIMS);
	HVHT = HT.transpose()*VT.inverse()*HT;
	MatrixXd VHT(DIMS, DIMS);
	VHT = VTI*HT;
	MatrixXd WFT(DIMS, DIMS);
	WFT = WTI*FT.transpose();
	MatrixXd WF(DIMS, DIMS);
	WF = WTI*FT;

	// Get intial values for posterior conditional
	MatrixXd AT_0(DIMS, DIMS);
	AT_0 = W0.inverse()+FWFT;
	MatrixXd AI_0(DIMS, DIMS);
	AI_0.setZero();
	VectorXd D_0(DIMS);
	D_0 = MU.transpose()*W0;

	double VTD = VT.determinant();
	double WTD = WT.determinant();
	double LG_W0D = -log(2*M_PI*W0.determinant());
	double W0_SQ = MU.transpose()*W0.inverse()*MU;
	double LG_AD = log(2*M_PI*AT_0.determinant());
	double ADSQ_0 = D_0.transpose()*AT_0*D_0;
	double POST_INIT = 0.5*(LG_AD+LG_W0D)+0.5*(ADSQ_0+W0_SQ);
	

	Forest ftest(FOREST);
	Tree * tree = ftest.fptr;	
	vector<int> leaf_nodes_up;
	vector<int> leaf_nodes_prop;
	vector<int> int_nodes;
	VectorXd datavec(DATA_LENGTH);

 	ifstream ifs("xsim_data2.txt", ios::in);

	// Get a line of data from .txt file and make predictor names
	MatrixXd datai(NSPLITS, DATA_LENGTH);
	datai.row(0) = read_line(0, DATA_LENGTH, ifs);
	vector<string> myNames = make_names(100, 10);

	// Get y
	VectorXd yi(DIMS);
	//yi = sim_y(data1, VT, A_mat, HT);
   
	// Use a k-loop to intialise forest
	for (int k = 0; k < FOREST; k++) {
		ftest.fvec[k] = Tree(1, DIMS);
		tree = &ftest.fvec[k];
		cout << "Address of ftest.fptr->root is: " << tree->root 
				<< endl;
		cout << "Address of ftest.fptr is: " << &tree << endl;
		ftest.init_tree(tree, AT_0, AI_0, ALPH, BET, D_0, D_0, 
				POST_INIT, myNames, datai, MU, SIG);
		//cout << "Initial Trees are:\n";
		//tree->display(tree->root);
		//cout << endl;
	}

	VectorXd qvec(FOREST);
	VectorXd postvec(FOREST);
	//VectorXd max_vec(FOREST);
	double max_post;
	double qstar;
	double const_prop;

	for (int l = 1; l < LOOPS; l++) {
		int i = l%NSPLITS;
		cout << endl;
		cout << "For iter = " << i << " and loop = " << l << "\n\n";
		double num_splits = nsplit(1, NSPLITS);
		double nprule = prule(DATA_LENGTH, num_splits);
		datai.row(i) = read_line(i*DATA_LENGTH, DATA_LENGTH*(i+1), ifs);
		datavec = datai.row(i);
		//cout << "Datai is:\n" << datai << endl;
		//cout << "Datavec is:\n" << datavec.transpose() << endl;
		yi = sim_y(datavec, VT, A_mat, HT);
		//cout << "yi is:\n" << yi << endl;
		Node * found;
		for (int k = 0; k < FOREST; k++) {
			cout << "The tree number is: " << k << endl;
			tree = &ftest.fvec[k];
			leaf_nodes_up = ftest.lnodes(tree);
			leaf_nodes_prop = leaf_nodes_up;
			int nleaves = leaf_nodes_up.size();
			tree->find_int_nodes(tree->root, int_nodes);
			cout << "The number of leaves is: " << 
						nleaves << endl;

			found = ftest.ch_node(tree, datavec, myNames);
			//cout << "Found node is: " << found->nnode << endl << endl;

			// Get updated Tree
			ftest.update(tree, found,
				leaf_nodes_up, DIMS, yi, WT, VT, WTI, 
				VTI, FT, HT, HVHT, VHT, WFT, WF, FWFT, 
				VTD, WTD);

			// Copy tree
			Tree copy(1, DIMS);
			copy = *tree;
			Tree * cptr = &copy; 
			//Tree * copy = tree; 
			
			//cout << "The copy before the proposal:\n";
			//cout << "Copy Posterior is: " << cptr->tree_post << endl;
			//cout << "Copy Prior is: " << cptr->prior << endl;
			//cout << "Copy prop_new is: " << cptr->prop_new << endl;
			//cout << "Copy prop_current is: " << cptr->prop_current 
			//		<< endl << endl;

			// Update tree prior
			tree->find_prior(tree->root);

			//cout << "The Tree after the update is:\n";
			//cout << "Tree Posterior is: " << tree->tree_post << endl;
			//cout << "Tree Prior is: " << tree->prior << endl;
			//cout << "Tree prop_new is: " << tree->prop_new << endl;
			//cout << "Tree prop_current is: " << tree->prop_current 
			//		<< endl << endl;

			//cout << "The copy after the prior update is:\n";
			//cout << "Copy Posterior is: " << cptr->tree_post << endl;
			//cout << "Copy Prior is: " << cptr->prior << endl;
			//cout << "Copy prop_new is: " << cptr->prop_new << endl;
			//cout << "Copy prop_current is: " << cptr->prop_current 
			//		<< endl << endl;
			//cout << " The tree copy is:\n";
			//copy->display(copy->root);
			//cout << endl;
			//cout << endl;
			// Propose new tree
			ftest.proposal(cptr, found, leaf_nodes_prop, int_nodes,
					ALPH, BET, nprule, DIMS, myNames, datai, i, WT, VT,
					WTI, VTI, FT, HT, HVHT, VHT, WFT, WF,
					FWFT, yi, VTD, WTD);
			if (nleaves > 1) {
			cout << "The copy after the proposal:\n";
			cout << "Copy Posterior is: " << cptr->tree_post << endl;
			cout << "Copy Prior is: " << cptr->prior << endl;
			cout << "Copy Prop_new is: " << cptr->prop_new << endl;
			cout << "Copy Prop_current is: " << cptr->prop_current 
					<< endl << endl;
			}
			int newtree;
			newtree = ftest.ch_tree(tree, cptr);
			if (newtree == 2) {
				*tree = copy;
				copy.~Tree();
			} 
			else { 
				copy.~Tree();
			}
			//cout << "The tree is:\n";
			//tree->display(tree->root);
			//cout << endl;	

		//	m.change(found->left, datai, myNames, i);
		//	tree->display(tree->root);	
		//	cout << endl;
		//	m.mswap(found);
		//	tree->display(tree->root);
		//	cout << endl;	
		//	m.prune(found, DIMS, datai, found->nnode, WT, VT, WTI, VTI,
		//			FT, HT, HVHT, VHT, WFT, WF, FWFT, yi, VTD, WTD);
		//	tree->display(tree->root);
		//	cout << endl;	
			//tree->find_post(tree->root);
			//tree->find_prior(tree->root);
			qvec(k) = log(tree->prior) + tree->tree_post;
			//max_vec(k) = tree->tree_post;
			if (nleaves > 1) {
			cout << "The Final tree:\n";
			cout << "Tree Posterior is: " << tree->tree_post << endl;
			cout << "Tree Prior is: " << tree->prior << endl;
			cout << "Tree prop_new is: " << tree->prop_new << endl;
			cout << "Tree prop_current is: " << tree->prop_current 
					<< endl << endl;
			} 
			//cout << "Log Tree Prior is:\n" << log(tree->prior) << endl;
			//cout << "State estimate is:\n" << tree->root->state << endl;
			//cout << "Leaf Posterior is:\n" << found->lpost << endl;

			tree->prior = 1.0;
			tree->tree_post = 0.0;
		}
		//cout << endl;
		//cout << "qvec is:\n" << qvec << endl;
		max_post = qvec.maxCoeff();
		//cout << "max_post is: " << max_post << endl;
		qvec = qvec.array() - max_post;
		//cout << "qvec diff is:\n" << qvec << endl;
		qvec = qvec.array().exp();
		//cout << "qvec end is:\n" << qvec << endl;
		const_prop = qvec.sum();
		postvec = qvec.array()/const_prop;
		cout << "The vector of posteriors for each tree is:\n" <<
				postvec.transpose() << endl;
	}
}*/
