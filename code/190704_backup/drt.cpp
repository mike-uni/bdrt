// PARALLEL EXAMPLE USING TREES

# include <iostream>
# include <Eigen/Dense>
# include <fstream>
# include <algorithm>
# include <vector>
# include <time.h>
# include <chrono>
# include <dirent.h>
# include "kalman.h"
# include "prior.h"
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
	string inpath;
	string outpath;
	string basemat;
	string wmat;
	string xsim;
	string ysim;
	string usim;
	string bounds;
	string sttype;
	string general;
	int LOOPS;
	int FOREST;
	int TREESIZE;
	int UPITER;
	double NUPPROB;
	int NDIMS;
	int MDIMS;
	int DATALENGTH;
	double ALPH;
	double BET;
	int nlevels;
	int startlevel;
	double c_0;
	double n_0;
	int dopara;
	int dosimtemp;
	vector<double> tempvec;
 
	// Input from file
	if (argc > 1) {
		ifstream fs(argv[1]);
		while(fs.peek() != EOF) {
			getline(fs, inpath);
			getline(fs, outpath);
			getline(fs, basemat);
			getline(fs, wmat);
			getline(fs, xsim);
			getline(fs, ysim);
			getline(fs, usim);
			getline(fs, bounds);
			getline(fs, general);
			LOOPS = stoi(general);
			getline(fs, general);
			FOREST = stoi(general);
			getline(fs, general);
			TREESIZE = stoi(general);
			getline(fs, general);
			UPITER = stoi(general);
			getline(fs, general);
			NUPPROB = stod(general);
			getline(fs, general);
			NDIMS = stoi(general);
			getline(fs, general);
			MDIMS = stoi(general);
			getline(fs, general);
			DATALENGTH = stoi(general);
			getline(fs, general);
			ALPH = stod(general);
			getline(fs, general);
			BET = stod(general);
			getline(fs, general);
			nlevels = stoi(general);
			getline(fs, general);
			startlevel = stoi(general);
			getline(fs, general);
			c_0 = stod(general);
			getline(fs, general);
			n_0 = stod(general);
			getline(fs, general);
			dopara = stoi(general);
			getline(fs, general);
			dosimtemp = stoi(general);
			getline(fs, sttype);
			getline(fs,general);
			istringstream vs(general);
			double tempval;
			while(vs >> tempval) {
				tempvec.push_back(tempval);
			}
		}
	}
	else {
		cout << "Forgotten File" << endl;
	}
	cout << inpath << " " << outpath << " " << basemat << " " << wmat << " " << xsim << " " << ysim << " " << usim << " " << 
		bounds << " " << LOOPS << " " << FOREST << " " << TREESIZE << " " << NDIMS << " " << MDIMS << " " << DATALENGTH << " " << 
		ALPH << " " << BET << " " << nlevels << " " << c_0 << " " << n_0 << " " << dopara << " " << dosimtemp <<
		" " << sttype << endl;

	Dataline D(DATALENGTH);
	vector<double> pseudoprior(nlevels);
	fill(pseudoprior.begin(), pseudoprior.end(), 1/double(nlevels));

	// Matrix variables
	// Eigen Setting for output
	IOFormat linemat(StreamPrecision, DontAlignCols, ", ", ", ", "", "", "", "");

	string bt = inpath+"/"+basemat;
	ifstream bfs(bt, ios::in);
	string inmat;

	MatrixXd W0(MDIMS, MDIMS);
	getline(bfs, inmat);
	W0 = D.read_matbyrow(inmat, MDIMS, MDIMS);
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

	// Get a different WT for each tree
	string wt = inpath+"/"+wmat;
	ifstream wfs(wt, ios::in);
	int inwmats = count(istreambuf_iterator<char>(wfs), istreambuf_iterator<char>(), '\n');
	wfs.clear();
	cout << "inwmats: " << inwmats << endl;
	wfs.seekg(0, ios::beg);
	vector<string> wtlines(inwmats);
	for (int i = 0; i < inwmats; ++i) {
		getline(wfs, wtlines[i]);
		cout << wtlines[i] << " ";
	}
	wfs.close();

	// Create addtional constant matrices from input matrices
	MatrixXd VTI(NDIMS, NDIMS);
	VTI = VT.inverse();
	MatrixXd W0I(MDIMS, MDIMS);
	W0I = W0.inverse();

	MatrixXd HVIHT(MDIMS, MDIMS);
	HVIHT = HT.transpose()*VTI*HT;
	MatrixXd HTVI(MDIMS, NDIMS);
	HTVI = HT.transpose()*VTI;
	
	const double W0_SQ = (MU.transpose()*W0I*MU);
	const double LG_W0D = log(2*M_PI*W0.determinant());
	
	// Vectors for input
	VectorXd GU(MDIMS);
	VectorXd yi(NDIMS);
	VectorXd ut(MDIMS);

	// Get input stream
	string xt = inpath+"/"+xsim;
	ifstream xfs(xt, ios::in);

	string yt = inpath+"/"+ysim;
	ifstream yfs(yt, ios::in);

	string uin = inpath+"/"+usim;
	ifstream ufs(uin, ios::in);


	// Get output name suffix
	string suffix = "_"+to_string(NDIMS)+to_string(MDIMS)+"_"+to_string(LOOPS)+
		"_"+to_string(FOREST)+"_"+to_string(DATALENGTH)+"_"+to_string(ALPH)+
		"_"+to_string(BET)+"_"+to_string(dopara)+"_"+to_string(dosimtemp)+"_"+sttype;
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
	for(auto s : D.names) {cout << s << " ";}
	cout << endl;

	string bin = inpath+"/"+bounds;
	ifstream bdfs(bin, ios::in);

	string typestring;
	getline(bdfs, typestring);
	
	string boundstring;
	getline(bdfs, boundstring);

	bfs.close();

	D.read_alltypes(typestring);
	D.read_bounds(boundstring);

	// Initialise the forest
	int k;
	Forest F;
	map<int, Tree*> copytrees;
	vector<int> leaf_nodes;
	vector<int> end_leaf_nodes;
	vector<int> int_nodes;
	vector<int> leaf_nodes_prop;
	double pgmove = 1/double(8);
	double ppmove = 1/double(8);
	double pfmove = 1/double(8);
	double prmove = 1/double(8);

	const int numthreads = omp_get_num_threads();
			
	// Begin loops over number of iterations
	for (int l = 0; l < LOOPS; l++) {
		chrono::high_resolution_clock::time_point init1, end1;
		init1 = chrono::high_resolution_clock::now();

		string newdata;
		Dataline::line fline;
		string newinput;
		string newresponse;
		if (l > 0) {	
			// Get new predictors
			getline(xfs, newdata);
			fline = D.read_line(newdata);
			
			// Get input
			getline(ufs, newinput);
			ut = D.read_vec(newinput, MDIMS);

			// Create input matrices
			if (ut.isZero(0)) {
				GU = VectorXd::Zero(MDIMS);
			}
			else {
				GU = GT*ut;
			}

			// Get y
			getline(yfs, newresponse);
			yi = D.read_vec(newresponse, NDIMS);
		}
		// Declare parallel varaibles
		int k, id;
		map<int, Node*>::iterator cit;
		vector<string>::const_iterator vit;
		vit = wtlines.begin();

#pragma omp parallel if(dopara == 1) shared(copytrees) firstprivate(leaf_nodes, \
		int_nodes, leaf_nodes_prop, VT, HT, FT, GT, LG_W0D, W0_SQ, W0, MU, \
	       	VTI, NDIMS, MDIMS, HVIHT, HTVI, ut, pgmove, ppmove, pfmove, prmove,\
		yi, GU, l, k, id, vit, cit, D, tempvec, c_0, n_0, \
		nlevels, pseudoprior, wtlines)	
	{	
#pragma omp for schedule(dynamic, 1) ordered 
		for (k = 1; k <= FOREST; k++) { 
			cout << endl;
			printf("iter: %d tree: %d\n", l, k);
			id = omp_get_thread_num();
			if (l == 0) {
				F.forest[k] = new Tree(1, MDIMS, NDIMS);
				copytrees[k] = new Tree(1, MDIMS, NDIMS);

				F.init_tree(F.forest[k]->tree, ALPH, BET, MU, SIG, DATALENGTH, D.names, D.alltypes, D.allcats, 
						D.dbounds, D.ibounds);
				F.forest[k]->tname = k;
				F.forest[k]->threadid = 0;
				F.forest[k]->numleaves = 1;
				F.forest[k]->move = 0;
				F.forest[k]->level = startlevel;
				F.forest[k]->accept = 0;
				F.forest[k]->found = 0;

			#pragma omp ordered 
				{	
					size_t nwtl = wtlines.size();
					F.forest[k]->WT = D.read_matbyrow(wtlines[(k-1)%nwtl], MDIMS, MDIMS);
					cout << "WT " << F.forest[k]->WT << endl;
					F.forest[k]->WI = F.forest[k]->WT.inverse();
					F.forest[k]->FTWIF = FT.transpose()*F.forest[k]->WI*FT;
					F.forest[k]->FTWI = FT.transpose()*F.forest[k]->WI;
					F.forest[k]->WIF = F.forest[k]->WI*FT;
					F.forest[k]->vdet = (2*M_PI*VT).determinant();
					F.forest[k]->wdet = (2*M_PI*F.forest[k]->WT).determinant();

					// Get intial values for posterior marginal
					map<int, Node*>::const_iterator tit;
					tit = F.forest[k]->tree.find(1);
					tit->second->atmat = W0.inverse()+F.forest[k]->FTWIF;
					tit->second->aimat = MatrixXd::Zero(MDIMS, MDIMS);
					tit->second->dvec = MU.transpose()*W0;
					tit->second->HTVYp = VectorXd::Zero(MDIMS);

					tit->second->adetp = log(2*M_PI*tit->second->atmat.determinant());
					tit->second->adsqp = (tit->second->dvec.transpose()*tit->second->atmat.inverse()*tit->second->dvec);
					tit->second->lpost = 0.5*(tit->second->adetp-LG_W0D-W0_SQ+tit->second->adsqp);
					
					// Output to file
					F.nodeout(F.forest[k]->tree, file1, k, 0);
					F.treeout(F.forest[k], file2, 0);
				}
			}
			else {		

			F.forest[k]->find_leaf_nodes(F.forest[k]->tree, leaf_nodes, F.forest[k]->numleaves);
			cout << "Main ln: ";
			for(auto s: leaf_nodes){cout << s << " ";}
			cout << endl;
			F.forest[k]->threadid = id;
			int nint = 0;
			F.forest[k]->find_int_nodes(F.forest[k]->tree, int_nodes, nint);
			size_t treesize = F.forest[k]->numleaves;
			
			// STAGE 1. COPY THE CURRENT TREE
			cit = F.forest[k]->tree.begin();
			while (cit != F.forest[k]->tree.end()) {
				copytrees[k]->tree[cit->first] = new Node;
				*copytrees[k]->tree[cit->first] = *cit->second;
				++cit;
			}

			// STAGE 2. CHOOSE A RANDOM LEAF FROM THE COPIED TREE
			// AND PROPOSE A MOVE.
			F.proposal(F.forest[k], copytrees[k]->tree, leaf_nodes, 
				int_nodes, ALPH, BET, MDIMS, NDIMS, DATALENGTH, D.names, D.alltypes,
				D.allcats, D.dbounds, D.ibounds, nlevels, TREESIZE, treesize, pgmove, 
				ppmove, NUPPROB, pfmove, prmove);
			cout << "Main: move " << F.forest[k]->move << endl;

			// STAGE 3. SELECT LEAF FROM BOTH TREES AND UPDATE THE POSTERIOR OF EACH
			cout << "Main display ltree: ";
			Tree T;
			T.display(F.forest[k]->tree);
			cout << endl;
			F.forest[k]->find_uleaf(F.forest[k]->tree, fline, F.forest[k]->found);
			cout << "Main display copytree: ";
			T.display(copytrees[k]->tree);
			cout << endl;
			F.forest[k]->find_uleaf(copytrees[k]->tree, fline, copytrees[k]->found);

			// 3.1 Update current tree (ltree).
			F.lupdate(F.forest[k]->tree, F.forest[k]->found, leaf_nodes, yi, F.forest[k]->WI, VTI, 
					HVIHT, HTVI, F.forest[k]->FTWI, F.forest[k]->WIF, F.forest[k]->FTWIF, 
					GU, F.forest[k]->vdet, F.forest[k]->wdet, F.forest[k]->slpostn, UPITER);

			// Update ltree prior
			F.forest[k]->find_prior(F.forest[k]->tree);

			// Update ltree logtreepost
			F.forest[k]->find_post(F.forest[k]->tree);

			// Calculate ltree zk loglikelihood
			double ltreezk = F.logzk(F.forest[k]->tree);

			// 3.2 Update proposed tree (copytree).
			F.forest[k]->find_leaf_nodes(copytrees[k]->tree, leaf_nodes_prop, copytrees[k]->numleaves);
			cout << "Main lnprop: ";
			for(auto s: leaf_nodes_prop){cout << s << " ";}
			cout << endl;
			F.lupdate(copytrees[k]->tree, copytrees[k]->found, leaf_nodes_prop, yi, F.forest[k]->WI, VTI, 
					HVIHT, HTVI, F.forest[k]->FTWI, F.forest[k]->WIF, F.forest[k]->FTWIF, 
					GU, F.forest[k]->vdet, F.forest[k]->wdet, copytrees[k]->slpostn, UPITER);

			// Update copy prior
			copytrees[k]->find_prior(copytrees[k]->tree);

			// Update copy logtreepost
			copytrees[k]->find_post(copytrees[k]->tree);

			// Calculate copy zk loglikelihood
			double copyzk = F.logzk(copytrees[k]->tree);

			// STAGE 4. CHOOSE TREE
			F.ch_tree(F.forest[k], F.forest[k]->prior, F.forest[k]->logtreepost,
				copytrees[k]->prior, copytrees[k]->logtreepost, F.forest[k]->accept,
				ltreezk, copyzk);

			if (F.forest[k]->accept == 1) {
				F.forest[k]->prior = copytrees[k]->prior;
				F.forest[k]->numleaves = copytrees[k]->numleaves;
				F.forest[k]->found = copytrees[k]->found;
				F.forest[k]->logtreepost = copytrees[k]->logtreepost;
				F.forest[k]->slpostn = copytrees[k]->slpostn;
				F.forest[k]->qstar = 0;
				F.forest[k]->tree.clear();
				cit = copytrees[k]->tree.begin();
				while (cit != copytrees[k]->tree.end()) {
					F.forest[k]->tree[cit->first] = new Node;
					*F.forest[k]->tree[cit->first] = *cit->second;
					++cit;
				}
				copytrees[k]->tree.clear();
				copytrees[k]->prior = 1.0;
				copytrees[k]->tname = 0;
				copytrees[k]->threadid = 0;
				copytrees[k]->numleaves = 0;
				copytrees[k]->move = 0;
				copytrees[k]->level = 0;
				copytrees[k]->accept= 0;
				copytrees[k]->found = 0;
				copytrees[k]->lfratio = 0;
				copytrees[k]->lrratio = 0;
				copytrees[k]->logtreepost = 0;
				copytrees[k]->slpostn = 0;
				copytrees[k]->qstar = 0;
			} 
			else { 
				copytrees[k]->tree.clear();
				copytrees[k]->prior = 1.0;
				copytrees[k]->tname = 0;
				copytrees[k]->threadid = 0;
				copytrees[k]->numleaves = 0;
				copytrees[k]->move = 0;
				copytrees[k]->level = 0;
				copytrees[k]->accept= 0;
				copytrees[k]->found = 0;
				copytrees[k]->lfratio = 0;
				copytrees[k]->lrratio = 0;
				copytrees[k]->logtreepost = 0;
				copytrees[k]->slpostn = 0;
				copytrees[k]->qstar = 0;
			}

			// Choose new temperature level
			int clevel = F.forest[k]->level;
			double cpseudo = pseudoprior[clevel-1];
			double npseudo;
			cpseudo = exp(log(cpseudo)-(c_0/(double(k+n_0))));
			for(int i = 0; i < nlevels; ++i) {
				npseudo = pseudoprior[i];
				npseudo = exp(log(npseudo)+(c_0/(nlevels*double(k+n_0))));
				pseudoprior[i] = npseudo;
			}
			pseudoprior[clevel-1] = cpseudo;
			F.ch_temp(F.forest[k], F.forest[k]->prior, F.forest[k]->logtreepost, tempvec, pseudoprior);

			// Add qstar to tree
			F.forest[k]->qstar = log(F.forest[k]->prior) + F.forest[k]->logtreepost
				+ F.forest[k]->slpostn;

			//STAGE 5. UPDATE THE KALMAN FILTERS OF THE SELECTED TREE
			leaf_nodes.clear();
			F.forest[k]->numleaves = 0;
			F.forest[k]->find_leaf_nodes(F.forest[k]->tree, leaf_nodes, F.forest[k]->numleaves);
			F.kupdate(F.forest[k]->tree, F.forest[k]->found, leaf_nodes, NDIMS, MDIMS, yi, F.forest[k]->WT, VT, 
					FT, HT, GT, ut);
#pragma omp ordered 
			{
			F.nodeout(F.forest[k]->tree, file1, k, l);
			F.treeout(F.forest[k], file2, l);
			}
			F.forest[k]->prior = 1.0;
			F.forest[k]->logtreepost = 0.0; // Updated in each iteration
			F.forest[k]->qstar = 0.0; // Updated in each iteration
			F.forest[k]->slpostn = 0.0; // reset sum of posterior of leaves not updated.
			int_nodes.clear();
			leaf_nodes.clear();
			leaf_nodes_prop.clear();
			F.forest[k]->numleaves = 0;
			F.forest[k]->found = 0;
		}
	}
	}
	end1 = chrono::high_resolution_clock::now();
	file5 << chrono::duration_cast<chrono::microseconds>(end1 - init1).count() << "\n";
	}
}
