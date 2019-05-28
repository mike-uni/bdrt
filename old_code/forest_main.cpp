//FOREST MAIN

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
# include "forest_class.h"

using namespace std;
using namespace Eigen;

int main() {
	// Input folder name
	string infolder = "cpp/180521_cpp/oned";

	// Output folder name
	string outfolder = "output";

	// Open files to write to:
	string path1 = "/home/michael/Desktop/"+infolder+"/"+outfolder+"/state.txt";
	string pathfile1 = path1;
	ofstream file1;
	file1.open(pathfile1);
	string header1 = "iter tree ysim1 ysim2 zhat1 zhat2 qstar\n";
	file1 << header1;

	string path2 = "/home/michael/Desktop/"+infolder+"/"+outfolder+"/moves.txt";
	string pathfile2 = path2;
	ofstream file2;
	string header2 = "iter tree move accept\n";
	file2.open(pathfile2);
	file2 << header2;

	string path3 = "/home/michael/Desktop/"+infolder+"/"+outfolder+"/posterior.txt";
	string pathfile3 = path3;
	ofstream file3;
	//string header3 = "iter posterior\n";
	file3.open(pathfile3);
	//file3 << header3;

	string path4 = "/home/michael/Desktop/"+infolder+"/"+outfolder+"/treeout.txt";
	string pathfile4 = path4;
	ofstream file4;
	file4.open(pathfile4);
	
	const int LOOPS = 100;
	const int FOREST = 20;
	const int DIMS = 1;
	const int DATA_LENGTH = 10;
	const int NSPLITS = 10;
	const double ALPH = 0.95;
	const double BET = 1.0;

	MatrixXd W0(DIMS, DIMS);
	W0 << 0.999;//,0,0,0.999;
	MatrixXd WT(DIMS, DIMS);
	WT << 1;//,0,0,1;
	MatrixXd VT(DIMS, DIMS);
	VT << 2;//,0,0,2;
	MatrixXd HT(DIMS, DIMS);
	HT << 1;//,0,0,1;
	MatrixXd FT(DIMS, DIMS);
	FT << -0.5;//,0,0,1;
	MatrixXd SIG(DIMS, DIMS);
	SIG << 1;//,0,0,1;
	VectorXd MU(DIMS);
	MU << 0;//,0;

/*	MatrixXd A_mat(DIMS, DATA_LENGTH);
		for (int i = 0; i < A_mat.size(); i++) {
    	A_mat(i) = rand()%2;
	}
*/
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
	double LPOST_INIT = 0.0;
	
	// Initialise the forest
	Forest ftest;
	Tree * tree = ftest.fptr;	
	vector<int> leaf_nodes_up;
	vector<int> leaf_nodes_prop;
	vector<int> int_nodes;
	VectorXd datavec(DATA_LENGTH);

	VectorXd qvec(FOREST);
	VectorXd postvec(FOREST);
	double max_post;
	double qstar;
	double const_prop;

	// Get input stream
	string data = "/home/michael/Desktop/"+infolder+"/"+outfolder+"/xsim.txt";
 	ifstream ifs(data, ios::in);
	string yt = "/home/michael/Desktop/"+infolder+"/"+outfolder+"/ysim.txt";
 	ifstream yfs(yt, ios::in);

	// Get a line of data from .txt file and make predictor names
	MatrixXd datai(NSPLITS, DATA_LENGTH);
	datai.row(0) = read_line(0, DATA_LENGTH, ifs);
	vector<string> myNames = make_names(100, DATA_LENGTH);

   
	// Use a k-loop to intialise forest trees
	for (int k = 0; k < FOREST; k++) {
		Tree itree(1, DIMS);
		ftest.fvec.push_back(itree);
		tree = &ftest.fvec[k];
		ftest.init_tree(tree, AT_0, AI_0, ALPH, BET, D_0, D_0, 
				POST_INIT, LPOST_INIT, myNames, datai, MU, SIG);
	}
	//cout << "Initialisation over.\n"; 

	// Begin loop over trees
	for (int l = 1; l < LOOPS; l++) {
		int i = l%NSPLITS;
		double num_splits = nsplit(1, NSPLITS);
		double nprule = prule(DATA_LENGTH, num_splits);
		file3 << l << " ";
		datai.row(i) = read_line(l*DATA_LENGTH, DATA_LENGTH*(l+1), ifs);
		datavec = datai.row(i);
		//cout << datavec.transpose() << endl;

		// Get y
		VectorXd yi(DIMS);
		//yi = sim_y(datavec, VT, A_mat, HT);
		yi = read_line(l*DIMS, DIMS*(l+1), yfs);
		//cout << yi.transpose() << endl;

		Node * found;
		for (int k = 0; k < FOREST; k++) {
			file1 << l << " " << k << " ";
			file2 << l << " " << k << " ";

			tree = &ftest.fvec[k];
			leaf_nodes_up = ftest.lnodes(tree);
			leaf_nodes_prop = leaf_nodes_up;
			int nleaves = leaf_nodes_up.size();
			tree->find_int_nodes(tree->root, int_nodes);

			found = ftest.ch_node(tree, datavec, myNames);

			ftest.update(tree, found, /*datavec, myNames, */
				leaf_nodes_up, DIMS, yi, WT, VT, WTI, 
				VTI, FT, HT, HVHT, VHT, WFT, WF, FWFT, 
				VTD, WTD);

			file1 << yi.transpose() << " ";
			file1 << found->state.transpose() << " ";

			// Copy tree
			Tree copy(1, DIMS);
			copy = *tree;
			Tree * cptr = &copy; 
			
			// Update tree prior
			tree->find_prior(tree->root);

			// Propose new tree
			ftest.proposal(cptr, found, leaf_nodes_prop, int_nodes,
					ALPH, BET, nprule, DIMS, myNames, datai, i, WT, VT,
					WTI, VTI, FT, HT, HVHT, VHT, WFT, WF,
					FWFT, yi, VTD, WTD);
			file2 << cptr->move << " ";
			int newtree;
			newtree = ftest.ch_tree(tree, cptr);
			if (newtree == 2) {
				file2 << 1 << "\n";
				*tree = copy;
			//	copy.~Tree();
			} 
			else { 
				file2 << 0 << "\n";
			//	copy.~Tree();
			}
			file4 << l << " " << k+1 << "   ";
			ftest.treeout(tree->root, file4);
			file4 << endl;
			//tree->display(tree->root);

			double qstar;
			qstar = log(tree->prior) + tree->tree_post;
			qvec(k) = qstar;
			if (isnan(qstar) || isinf(qstar)) {
				qstar = 99;
			}
			file1 << qstar << "\n";
			//cout << qstar << "\n";
			tree->prior = 1.0;
			tree->tree_post = 0.0;
			int_nodes.clear();
		}
		//cout << endl;
		max_post = qvec.maxCoeff();
		qvec = qvec.array() - max_post;
		qvec = qvec.array().exp();
		const_prop = qvec.sum();
		//cout << const_prop << endl;
		postvec = qvec.array()/const_prop;
		//file3 << tree_index.transpose() << "\n"; 
		file3 << postvec.transpose() << "\n";
		//cout << "The tree posteriors at " << l << " are is:\n" <<
		//		postvec.transpose() << endl;
		cout << " " << l << " " << flush;
	}
	cout << endl;
}
