// BASIC PARALLEL EXAMPLE USING TREES

# include <iostream>
# include <Eigen/Dense>
# include <fstream>
# include <algorithm>
# include <vector>
# include <time.h>
# include <chrono>
# include <dirent.h>
# include "kalman.h"
# include "priorf.h"
# include "posterior.h"
# include "node.h"
# include "tree.h"
# include "move.h"
# include "proposal.h"
# include "forest.h"
# include "read.h"

using namespace std;
using namespace Eigen;

int main(int argc, char * argv[]) {
 
	// choose parallelisation
	istringstream s1(argv[1]);	
	istringstream s2(argv[2]);	
	istringstream s3(argv[3]);	
	int dopara; 
	int LOOPS;
	int FOREST;
	s1 >> dopara;
	s2 >> LOOPS;
	s3 >> FOREST;

	cout << "LOOPS " << LOOPS << " FOREST " << FOREST << endl;


	// Declare path names
	string outpath = "/home/michael/cpp/190426/output/testout";
	string inpath = "/home/michael/cpp/190426/output/testin";

	// Basic variables
	const int NDIMS = 2;
	const int MDIMS = 2;
	const int DATALENGTH = 5;
	Dataline D(DATALENGTH);
	const double ALPH = 0.95;
	const double BET = 1.0;

	// Matrix variables
	// Eigen Setting for output
	IOFormat linemat(StreamPrecision, DontAlignCols, ", ", ", ", "", "", "", "");

	string basemat = "basemat2d.txt";
	string bt = inpath+"/"+basemat;
	cout << "bt " << bt << endl;
	ifstream bfs(bt, ios::in);
	string inmat;

	MatrixXd W0(MDIMS, MDIMS);
	getline(bfs, inmat);
	W0 = D.read_matbyrow(inmat, MDIMS, MDIMS);
	MatrixXd WT(MDIMS, MDIMS);
	getline(bfs, inmat);
	WT = D.read_matbyrow(inmat, MDIMS, MDIMS);
	MatrixXd VT(NDIMS, NDIMS);
	getline(bfs, inmat);
	VT = D.read_matbyrow(inmat, NDIMS, NDIMS);
	MatrixXd HT(NDIMS, MDIMS);
	getline(bfs, inmat);
	HT = D.read_matbyrow(inmat, NDIMS, MDIMS);
	MatrixXd FT(MDIMS, MDIMS);
	getline(bfs, inmat);
	FT = D.read_matbyrow(inmat, MDIMS, MDIMS);
	MatrixXd GT(MDIMS, MDIMS);
	getline(bfs, inmat);
	GT = D.read_matbyrow(inmat, MDIMS, MDIMS);
	MatrixXd SIG(MDIMS, MDIMS);
	getline(bfs, inmat);
	SIG = D.read_matbyrow(inmat, MDIMS, MDIMS);
	VectorXd MU(MDIMS);
	getline(bfs, inmat);
	MU = D.read_vec(inmat, MDIMS);

	// Create addtional constant matrices from input matrices
	MatrixXd VTI(NDIMS, NDIMS);
	VTI = VT.inverse();
	MatrixXd WTI(MDIMS, MDIMS);
	WTI = WT.inverse();
	MatrixXd W0I(MDIMS, MDIMS);
	W0I = W0.inverse();

	MatrixXd FWIFT(MDIMS, MDIMS);
	FWIFT = FT.transpose()*WTI*FT;
	MatrixXd HVIHT(MDIMS, MDIMS);
	HVIHT = HT.transpose()*VTI*HT;
	MatrixXd VIHT(NDIMS, MDIMS);
	VIHT = VTI*HT;
	MatrixXd WIFT(MDIMS, MDIMS);
	WIFT = (WTI*FT).transpose();
	MatrixXd WIF(MDIMS, MDIMS);
	WIF = WTI*FT;
	
	// Vectors for input
	VectorXd UPGWI(MDIMS);
	VectorXd GUWI(MDIMS);
	VectorXd yi(NDIMS);

	// Get intial values for posterior marginal
	MatrixXd AT_0(MDIMS, MDIMS);
	AT_0 = W0.inverse()+FWIFT;
	MatrixXd AI_0(MDIMS, MDIMS);
	AI_0 = MatrixXd::Zero(MDIMS, MDIMS);
	VectorXd D_0(MDIMS);
	D_0 = MU.transpose()*W0;
	VectorXd C_0(MDIMS);
	C_0 = VectorXd::Zero(MDIMS);

	const double VTD = VT.determinant();
	const double WTD = WT.determinant();
	const double LG_W0D = log(W0.determinant());
	const double LG_AD = log(AT_0.determinant());
	const double W0_SQ = (MU.transpose()*W0I*MU);
	const double ADSQ_0 = (D_0.transpose()*AT_0.inverse()*D_0);
	const double POST_INIT = (double)0.5*(LG_AD - LG_W0D) - (double)0.5*(W0_SQ - ADSQ_0);
	const double LPOST_INIT = 0.0;

	// Get input stream
	string infile1 = "xsim.txt";
	string xt = inpath+"/"+infile1;
	cout << "xt " << xt << endl;
	ifstream xfs(xt, ios::in);

	string infile2 = "ysim2d.txt";
	string yt = inpath+"/"+infile2;
	cout << "yt " << yt << endl;
	ifstream yfs(yt, ios::in);

	string infile3 = "usim2d.txt";
	string uin = inpath+"/"+infile3;
	cout << "uin " << uin << endl;
	ifstream ufs(uin, ios::in);

	// Get output name suffix
	string suffix = "_"+to_string(NDIMS)+to_string(MDIMS)+"_"+to_string(LOOPS)+
		"_"+to_string(FOREST)+"_"+to_string(DATALENGTH)+"_"+to_string(ALPH)+
		"_"+to_string(BET)+"_"+to_string(dopara);
	suffix.erase(remove(suffix.begin(), suffix.end(), '.'), suffix.end());
	suffix = suffix+".txt";

	// Open files to write to
	string path1 = outpath+"/"+"nodes"+suffix;
	ofstream file1;
	file1.open(path1);

	string path2 = outpath+"/"+"trees"+suffix;
	ofstream file2;
	file2.open(path2);

	string path5 = outpath+"/"+"time"+suffix;
	ofstream file5;
	file5.open(path5);

	// Fill base number of observations using
	D.make_names();
	for(string s : D.names) {
		cout << s << " ";
	}
	cout << endl;

	string typestring;
	getline(xfs, typestring);
	
	string boundstring;
	getline(xfs, boundstring);

	D.read_alltypes(typestring);
	D.read_bounds(boundstring);
	for(auto s : D.alltypes) {
		cout << s.second << " ";
	}
	cout << endl;
	for(auto s : D.dbounds) {
		cout << s.second[0] << " ";
	}
	cout << endl;


	// Get this first line of input
	VectorXd ut(MDIMS);
	VectorXd utp(MDIMS);
	utp = VectorXd::Zero(MDIMS);

	// Initialise the forest
	int k;
	Forest F;
	map<int, Tree*> copytrees;
	vector<int> leaf_nodes;
	vector<int> end_leaf_nodes;
	vector<int> int_nodes;
	vector<int> leaf_nodes_prop;
//	string detail = "full";

	double slpostn = 0.0;

	const int numthreads = omp_get_num_threads();
			
	// Begin loops over number of iterations
	for (int l = 0; l < LOOPS; l++) {
		chrono::high_resolution_clock::time_point init1, end1;
		init1 = chrono::high_resolution_clock::now();

//		int i = l%NSPLITS;
		//double num_splits = nsplit(l, NSPLITS);
//		double nprule = prule(DATALENGTH, NSPLITS);
		string newdata;
		getline(xfs, newdata);
		Dataline::line fline;
		fline = D.read_line(newdata);
		//datai.row(i) = read_line(l*DATALENGTH, DATALENGTH*(l+1), xfs);
		//datavec = datai.row(i);
		
		// Get input
		string newinput;
		getline(ufs, newinput);
		ut = D.read_vec(newinput, MDIMS);

		// Create input matrices
		if (ut.isZero(0)) {
			UPGWI = VectorXd::Zero(MDIMS);
			GUWI = VectorXd::Zero(MDIMS);
		}
		else {
			UPGWI = utp.transpose()*GT.transpose()*WTI;
			GUWI = ut.transpose()*GT.transpose()*WTI;
		}

		// Get y
		string newresponse;
		getline(yfs, newresponse);
		yi = D.read_vec(newresponse, NDIMS);
		cout << "yi " << yi << endl;

		// Declare parallel varaibles
		int k, id;
		map<int, Node*>::iterator cit;

#pragma omp parallel if(dopara == 1) shared(copytrees) firstprivate(leaf_nodes, \
		int_nodes, leaf_nodes_prop, slpostn, WT, VT, HT, FT, GT, \
	       	VTI, WTI, NDIMS, MDIMS, FWIFT, HVIHT, VIHT, WIFT, WIF, ut, \
		yi, UPGWI, GUWI, VTD, WTD, l, k, id, cit)	
	{	
#pragma omp for schedule(dynamic, 1) ordered 
		for (k = 1; k <= FOREST; k++) { 
			cout << endl;
			printf("iter: %d tree: %d\n", l, k);
			id = omp_get_thread_num();
	if (l == 0) {
			F.forest[k] = new Tree(1, MDIMS);
			copytrees[k] = new Tree(1, MDIMS);
			F.init_tree(F.forest[k]->tree, AT_0, AI_0, ALPH, BET, C_0, D_0, 
				POST_INIT, LPOST_INIT, MU, SIG, DATALENGTH, D.names, D.alltypes, 
				D.allcats, D.dbounds, D.ibounds);
			F.forest[k]->tname = k;
			F.forest[k]->threadid = 0;
			F.forest[k]->numleaves = 1;
			F.forest[k]->move = 0;
			F.forest[k]->accept = 0;
			F.forest[k]->found = 0;
			F.forest[k]->cfound = 0;

		#pragma omp ordered 
		{
			F.nodeout(F.forest[k]->tree, file1, k, 0);
			F.treeout(F.forest[k], file2, 0);
		}
			printf("iter: %d tree: %d\n", 0, k);
			F.forest[k]->display(F.forest[k]->tree);
			F.forest[k]->numleaves = 0;
	}
	else {	/*	map<int, Node*> ltree;
			map<int, Node*>::iterator cit = F.forest[k]->tree.begin();
			while (cit != F.forest[k]->tree.end()) {
				ltree[cit->first] = new Node;
				*ltree[cit->first] = *cit->second;
				++cit;
			}
		*/
			F.forest[k]->find_leaf_nodes(F.forest[k]->tree, leaf_nodes, F.forest[k]->numleaves);
			cout << "nleaves: " << endl;
			for(auto s : leaf_nodes){cout << s << " ";}
			cout << endl;
			F.forest[k]->threadid = id;
			int nint = 0;
			F.forest[k]->find_int_nodes(F.forest[k]->tree, int_nodes, nint);
			cout << "nint: " << endl;
			for(auto s : int_nodes){cout << s << " ";}
			cout << endl;
			
			//STAGE 1. CHOOSE A LEAF USING CURRENT STATE SPACE
			// AND CURRENT DATA. PERFORM INTERMITTENT UPDATE ON LEAVES
			F.forest[k]->find_uleaf(F.forest[k]->tree, fline, F.forest[k]->found);
			cout << "found: " << F.forest[k]->found << endl;

			// KALMAN UPDATE 
			F.kupdate(F.forest[k]->tree, F.forest[k]->found, leaf_nodes, NDIMS, MDIMS, yi, WT, VT, 
					FT, HT, GT, ut);
			
			// COPY THE CURRENT TREE
			//map<int, Node*> ctree;
			cit = F.forest[k]->tree.begin();
			while (cit != F.forest[k]->tree.end()) {
				copytrees[k]->tree[cit->first] = new Node;
				*copytrees[k]->tree[cit->first] = *cit->second;
				++cit;
			}

			// STAGE 2. CHOOSE A RANDOM LEAF FROM THE CURRENT (COPIED) TREE
			// AND PROPOSE A MOVE.
			F.proposal(F.forest[k], copytrees[k]->tree, leaf_nodes, 
				int_nodes, ALPH, BET, MDIMS, DATALENGTH, D.names, D.alltypes,
				D.allcats, D.dbounds, D.ibounds);
			cout << "ctree: \n";
			F.forest[k]->display(copytrees[k]->tree);
			cout << "ltree: \n";
			F.forest[k]->display(F.forest[k]->tree);
			cout << "move: " << F.forest[k]->move << endl;

			// STAGE 3. RESELECT LEAF USING CURRENT DATA AND PROPOSED TREE
			// AND UPDATE POSTERIOR DATA FOR BOTH TREES
			F.forest[k]->find_uleaf(copytrees[k]->tree, fline, F.forest[k]->cfound);
			cout << "cfound: " << F.forest[k]->cfound << endl;

			// 3a. Update current tree (ltree).
			F.lupdate(F.forest[k]->tree, F.forest[k]->found, leaf_nodes, yi,
					WTI, VTI, HVIHT, VIHT, WIFT, WIF, FWIFT, 
					UPGWI, GUWI, VTD, WTD, slpostn);

			// Update ltree prior
			F.forest[k]->find_prior(F.forest[k]->tree);

			// Update ltree logtreepost
			F.forest[k]->find_post(F.forest[k]->tree);

			// 3b. Update proposed tree (ltree).
			F.forest[k]->find_leaf_nodes(copytrees[k]->tree, leaf_nodes_prop, F.forest[k]->numcopyleaves);
			F.lupdate(copytrees[k]->tree, F.forest[k]->cfound, leaf_nodes_prop, yi,  
				WTI, VTI, HVIHT, VIHT, WIFT, WIF, FWIFT, 
				UPGWI, GUWI, VTD, WTD, slpostn);

			// Update copy prior
			copytrees[k]->find_prior(copytrees[k]->tree);

			// Update copy logtreepost
			copytrees[k]->find_post(copytrees[k]->tree);

			// STAGE 4. CHOOSE TREE AND CALCULATE POSTERIORS
			F.ch_tree(F.forest[k], F.forest[k]->prior, F.forest[k]->logtreepost,
				copytrees[k]->prior, copytrees[k]->logtreepost, F.forest[k]->accept);
			if (F.forest[k]->accept == 1) {
				F.forest[k]->tree.clear();
				cit = copytrees[k]->tree.begin();
				while (cit != copytrees[k]->tree.end()) {
					F.forest[k]->tree[cit->first] = new Node;
					*F.forest[k]->tree[cit->first] = *cit->second;
					++cit;
				}
			} 
			else { 
				copytrees[k]->tree.clear();
			}

/*			F.forest[k]->tree.clear();
			cit = F.forest[k].begin();
			while (cit != F.forest[k].end()) {
				F.forest[k]->tree[cit->first] = new Node;
				*F.forest[k]->tree[cit->first] = *cit->second;
				++cit;
			}
*/			
			F.forest[k]->find_leaf_nodes(F.forest[k]->tree, end_leaf_nodes, F.forest[k]->numendleaves);
#pragma omp ordered 
			{
			F.nodeout(F.forest[k]->tree, file1, k, l);
			F.treeout(F.forest[k], file2, l);
			}
			F.forest[k]->prior = 1.0;
			F.forest[k]->logtreepost = 0.0; // Updated in each iteration
			slpostn = 0.0; // reset sum of posterior of leaves not updated.
			int_nodes.clear();
			leaf_nodes.clear();
			leaf_nodes_prop.clear();
			end_leaf_nodes.clear();
			F.forest[k]->display(F.forest[k]->tree);
			//cout << "endleaves: " << F.forest[k]->numendleaves << endl;
			F.forest[k]->numleaves = 0;
			F.forest[k]->numendleaves = 0;
			F.forest[k]->numcopyleaves = 0;
			F.forest[k]->numintnodes = 0;
			F.forest[k]->found = 0;
			F.forest[k]->cfound = 0;
		}
	}
	}
	utp = ut;
	end1 = chrono::high_resolution_clock::now();
	file5 << chrono::duration_cast<chrono::microseconds>(end1 - init1).count() 
	<< "\n";
	}
}
