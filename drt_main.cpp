// THIS IS THE MAIN FILE FOR BAYESIAN DYNAMIC REGRESSION TREES

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
# include "settings_class.h"
# include "stree_class.h"
# include "snode_class.h"
# include "sim_class.h"

using namespace std;
using namespace Eigen;
using namespace libconfig;

int main() {	
	Config cfg;
	Settings set;
	Dist d;
	Simtree simtree;
	
	// Introduction
	cout << "Welcome to the BDRT!\n"
		"There are two main types of experiment:\n"
		"\t1. Those where model parameters are varied, for instance,\n"
		"the number of trees in the ensemble changes, but the\n"
		"'known' parameters (initial conditions) are not varied.\n"
		"\t2. Those where the input parameters are varied, for instance,\n"
		"the signal to noise ratio (W/V) changes, but most of the model\n"
		"parameters are kept constant\n";
	cout << "In what follows, option 1 above is referred to as the 'basic_settings'\n"
		"and the second option is the 'matrix_settings'. To limit the number\n"
		"of factors each of these will have a separate config file and a\n"
		"type of experiment must be chosen.\n";

	// Display current cfg files
	cout << "These are the custom config files stored in the current\n"
		"working directory:\n";

	// Get the current working directory and place in
	// settings class
	string cwd = get_wd();
	set.outpath = cwd;

	// Declare path names
	string outpath;
	string inpath;

	// Show .cfg files in cwd	
	get_exts(cwd, ".cfg");

	cout << "If you wish to use one of these config files please enter the name, \n"
		"if you wish to modify an existing one, please type 'mod'\n"
		"or if you wish to create a new one, please type 'new'\n";
	string exname;
	while (cin.peek() != '\n') {
		cin >> exname;
		if (cin.peek() == '\b') {
			cin.clear();
		}
	}
	cin.ignore();
	if (exname == "new") {
		// Choose type of experiment
		cout << "As you have selected to create a new config file you will \n"
			"need to choose the type of experiment you wish to conduct.\n"
			"Please enter either either 1 or 2 corresponding to the options\n"
			"described in the introduction and press enter.\n";
		int opt;
		while (cin.peek() != '\n') {
			cin >> opt;
			if (cin.peek() == '\b') {
				cin.clear();
			}
		}
		cin.ignore();

		// Get input folder name for this run of experiments
		cout << "Please enter the name of the separate folders in which\n" 
"the predictors and responses are currently being stored. " 
"If they do not yet exist, they will be created.\n";
		cout << "First the input folder:\n";
		string infoldername;
		while (cin.peek() != '\n') {
			cin >> infoldername;
			if (cin.peek() == '\b') {
				cin.clear();
			}
		}
		cin.ignore();
		
		set.infolder = infoldername;

		// Make new folder path 
		inpath = cwd+"/"+infoldername;

		// Make new folder
		set.make_folder(set.infolder);

		// Get folder name for output
		cout << "Now the output folder.\n";
		string foutname;
		while (cin.peek() != '\n') {
			cin >> foutname;
			if (cin.peek() == '\b') {
				cin.clear();
			}
		}	
		cin.ignore();
		set.outfolder = foutname;
	
		// Make new folder
		set.make_folder(set.outfolder);

		// Make new folder path
		outpath = cwd+"/"+foutname;

		// Set the config file name
		cout << "Please enter the name of the new .cfg file (without extension): "; 
		string newname;
		while (cin.peek() != '\n') {
			cin >> newname;
			if (cin.peek() == '\b') {
				cin.clear();
			}
		}	
		cin.ignore();
		set.outname = newname+".cfg";
		cout << endl;
	
		// Make cfg file based on type of experiment
		if (opt == 1) {
			// Set the number of matrix examples
			set.nmat = 1;

			// Show the example output file
			set.example();
			cout << endl;

			// Make the cfg file
			set.make_file(cfg);
			cin.ignore();
		}
		else if (opt == 2) {
			// Set the number of matrix examples
			cout << "Please enter the number of matrix examples: \n";
			int nmats;
			cin >> nmats;
			set.nmat = nmats;

			// Show the example output file
			set.examplemat();
			cout << endl;

			// Make the cfg file
			set.make_file(cfg);
			cin.ignore();
		}

		// Simulate Data
		cout << "You will now simulate data based on the configuration\n"
			"you have created. This will be used in the BDRT model below.\n";
		simtree.simcfg(cfg, set, d, exname, set.outname);
	}
	else if (exname == "mod") {
		cout << "Please enter the name of the file you wish to modify\n"
			"from those listed above:\n";
		string modname;
		while (cin.peek() != '\n') {
			cin >> modname;
			if (cin.peek() == '\b') {
				cin.clear();
			}
		}	
		cin.ignore();
		modname = modname+".cfg";

		try {
			cfg.readFile(modname.c_str());
		}
		catch(const FileIOException &fioex) {
			cerr << "I/O error while reading file." << endl;
			return (EXIT_FAILURE);
		}
		catch(const ParseException &pex) {
			cerr << "Parse error at " << pex.getFile() << ":"
				<< pex.getLine() << " - " << pex.getError() << endl;
			return(EXIT_FAILURE);
		}

		// Get the base of the configuration file
		Setting& modbase = cfg.getRoot();

		try {
			Setting &setname = modbase["setting_name"];
			int namel = setname.getLength();
			cout << "The setting names are:\n";
			for (int sn = 0; sn < namel; sn++) {
				string sname = setname[sn];
				cout << sname << endl;
			}
		}
		catch(const SettingNotFoundException &nfex) {
			cerr << "No 'setting_name' setting in config file." << endl;
		}
		
		cout <<"Please enter the setting you would like to modify:\n";
		string modset;
		bool right = true;
		while (right) {
			while (cin.peek() != '\n') {
				cin >> modset;
				if (cin.peek() == '\b') {
					cin.clear();
				}
			}
			Setting &upmod = modbase[modset];
			if (modset == "basic_setup") {
				try {
					int blen = upmod.getLength();
					for (int i = 0; i<blen; i++) {
						upmod.remove(i);
					}
					set.bsetup(1);

					Setting& curinfo = modbase["setting_info"];
					set.nmat = curinfo["nmat"];
				
					modbase.remove("setting_info");
				
					Setting& setinfo = modbase.add("setting_info", 
										Setting::TypeGroup);
					setinfo.add("dim", Setting::TypeInt) = set.dim;
					setinfo.add("nloops", Setting::TypeInt) = set.nloops;
					setinfo.add("nforests", Setting::TypeInt) = set.nforest;
					setinfo.add("nrest", Setting::TypeInt) = set.nrest;
					setinfo.add("nmat", Setting::TypeInt) = set.nmat;

					set.bset_fill(upmod);

					right = false;
				}
				catch(const SettingNotFoundException &nfex) {
				//	cerr << "Setting not found. Please try again."<< endl;
				}
			}
			else if (modset == "matrix_setup") {
				try{
					int blen = upmod.getLength();
					for (int i = 0; i<blen; i++) {
						upmod.remove(i);
					}
					cout << "Please enter the number of matrix examples: \n";
					int nmats;
					cin >> nmats;
					set.nmat = nmats;

					set.mat_setup(upmod, set.nmat);

					right = false;
				}
				catch(const SettingNotFoundException &nfex) {
				//	cerr << "Setting not found. Please try again."<< endl;
				}
			}
			else {
				cout << "This setting is unavailable. Please choose another: ";		
				while (cin.peek() != '\n') {
					cin >> modset;
					if (cin.peek() == '\b') {
						cin.clear();
					}
				}
			}
		}
		const char * nout = modname.c_str();
		// Write out the new configuration.
		try {
			cfg.writeFile(nout);
			cerr << "New config file written to: " << modname << endl;
		}
		catch (const FileIOException &fioex) {
			cerr << "I/O error while writing file: " << modname << endl;
		}
		Setting &filemod = modbase["file_path"];

		// Make new input folder path 
		string infolder = filemod.lookup("infolder");
		set.infolder = infolder;
		inpath = cwd+"/"+set.infolder;

		// Make new output folder path 
		string outfolder = filemod.lookup("outfolder");
		set.outfolder = outfolder;
		outpath = cwd+"/"+set.outfolder;

		set.outname = modname;
		// Re-simulate Data
		cout << "You will now simulate data based on the configuration\n"
			"you have modified. This will be used in the BDRT model below.\n";
		cin.ignore();
		simtree.simcfg(cfg, set, d, exname, set.outname);
	}
	else {
		set.outname = exname+".cfg";
		//cin.ignore();
	}
	cin.ignore();
	cout << "The current config file is "<< set.outname <<"." << endl;
	cout << "The config is now complete. To continue onto the BDRT please type 'yes', else 'exit' will exit the program:\n";

	string cont;
	while (cin.peek() != '\n') {
		cin >> cont;
		if (cont == "yes") break;
		else return(EXIT_SUCCESS);
	}
	cin.ignore();		
	// Read the config file. If there is an error, 
	// report it and exit;
	try {
		cfg.readFile(set.outname.c_str());
	}
	catch(const FileIOException &fioex) {
		cerr << "I/O error while reading file." << endl;
		return (EXIT_FAILURE);
	}
	catch(const ParseException &pex) {
		cerr << "Parse error at " << pex.getFile() << ":"
			<< pex.getLine() << " - " << pex.getError() << endl;
		return(EXIT_FAILURE);
	}

	// Get program intro.
	try {
		string model_info = cfg.lookup("model_info");
		cout << model_info << endl;
	}
	catch(const SettingNotFoundException &nfex) {
		cerr << "No 'model_info' setting in config file." << endl;
	}

	// Get the base of the configuration file
	const Setting& base = cfg.getRoot();

	// Get paths
	const Setting& fileinfo = base["file_path"];
	string outfolder = fileinfo.lookup("outfolder");
	string infolder = fileinfo.lookup("infolder");
	inpath = cwd+"/"+infolder;
	outpath = cwd+"/"+outfolder;

	// Numbers of iterations of different settings
	const Setting& setinfo = base["setting_info"];
	int DIMS = setinfo["dim"];
	int nLOOPS = setinfo["nloops"];
	int nFORESTS = setinfo["nforests"];
	int nDATA = setinfo["nrest"];
	int nMATS = setinfo["nmat"];
	
	int NUMEXP = nLOOPS*nFORESTS*nDATA;

	// Basic variables
	const Setting& basic_setup = base["basic_setup"];
	int LOOPS;
	int FOREST;
	int DATA_LENGTH;
	int pDATA_LENGTH = basic_setup[0][3];
	int NSPLITS;
	double ALPH;
	double BET;

	// Matrix variables
	const Setting &matrix_list = base["matrix_setup"];
	double mymatarr[DIMS*DIMS];
	double myvecarr[DIMS];

	MatrixXd W0(DIMS, DIMS);
	MatrixXd WT(DIMS, DIMS);
	MatrixXd VT(DIMS, DIMS);
	MatrixXd HT(DIMS, DIMS);
	MatrixXd FT(DIMS, DIMS);
	MatrixXd SIG(DIMS, DIMS);
	VectorXd MU(DIMS);

	int nd = 0;
	int step = nLOOPS*nFORESTS;
	for (int rn = 0; rn < NUMEXP; rn++) {

		// Get iterations of data variable xor number of matrices
		if (rn == step) {
			nd++;
			step = step+(nLOOPS*nFORESTS);
		}
		cout << "nd is: " << nd << " and step is: " << step 
			<< " and nMATS is: " << nMATS << endl;		
		
		for (int rm = 0; rm < nMATS; rm++) {
			// Set up matix inputs for type of matrix.
			const Setting &matrix_setup = matrix_list[rm];

			// Matrix Variables
			set.getdArray(mymatarr, matrix_setup[0], 0, DIMS*DIMS);
			W0 = Map<MatrixXd>(mymatarr, DIMS, DIMS);
			//cout << W0 << endl;

			set.getdArray(mymatarr, matrix_setup[1], 0, DIMS*DIMS);
			WT = Map<MatrixXd>(mymatarr, DIMS, DIMS);
			//cout << WT << endl;

			set.getdArray(mymatarr, matrix_setup[2], 0, DIMS*DIMS);
			VT = Map<MatrixXd>(mymatarr, DIMS, DIMS);
			//cout << VT << endl;

			set.getdArray(mymatarr, matrix_setup[3], 0, DIMS*DIMS);
			HT = Map<MatrixXd>(mymatarr, DIMS, DIMS);
			//cout << HT << endl;

			set.getdArray(mymatarr, matrix_setup[4], 0, DIMS*DIMS);
			FT = Map<MatrixXd>(mymatarr, DIMS, DIMS);
			//cout << FT << endl;

			set.getdArray(mymatarr, matrix_setup[5], 0, DIMS*DIMS);
			SIG = Map<MatrixXd>(mymatarr, DIMS, DIMS);
			//cout << SIG << endl;

			set.getdArray(myvecarr, matrix_setup[6], 0, DIMS);
			MU = Map<VectorXd>(myvecarr, DIMS);
			//cout << MU << endl;

			// Create addtional constant matrices from input matrices
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

			//cout << "LPOST_INIT: " << LPOST_INIT << " POST_INIT: " << POST_INIT
			//	<< endl;

			// Initialise Base Experiment varaibles
			LOOPS = basic_setup[rn][1];
			FOREST = basic_setup[rn][2];
			DATA_LENGTH = basic_setup[rn][3];
			NSPLITS = basic_setup[rn][4];
			ALPH = basic_setup[rn][5];
			BET = basic_setup[rn][6];

			cout << "LOOPS: " << LOOPS << " FOREST: " << FOREST
			 << " DATA_LENGTH: " << DATA_LENGTH << " NSPLITS: " << NSPLITS
			 << " ALPHA: " << ALPH << " BETA: " << BET
			 << endl;

			if (nMATS > 1 || pDATA_LENGTH == DATA_LENGTH) {
				nd = rm;
			}

			cout << "If nmat = 1, or pdl == dl, nd is: " << nd 
				<< " and rm is: " << rm << endl; 

			// Get input stream
			const Setting &xfiles = base["file_path"]["infilesx"];
			string infile1 = xfiles[nd];
			string indata = inpath+"/"+infile1; //xsim.txt
			cout << indata << endl;
			ifstream ifs(indata, ios::in);
		
			const Setting &yfiles = base["file_path"]["infilesy"];
			string infile2 = yfiles[nd];
			string yt = inpath+"/"+infile2; //ysim.txt
			cout << yt << endl;
			ifstream yfs(yt, ios::in);
		
			// Get output name suffix
			string suffix = "_"+to_string(LOOPS)+"_"+to_string(FOREST)+
							"_"+to_string(DATA_LENGTH)+"_"+to_string(NSPLITS)+
							"_"+to_string(ALPH)+"_"+to_string(BET)+
							"_"+to_string(rm);
			suffix.erase(remove(suffix.begin(), suffix.end(), '.'), suffix.end());
			suffix = suffix+".txt";
			//string delim = "_";
			//size_t start = infile1.find(delim);
			//string suffix = infile1.substr(start, infile1.length());
				
			// Open files to write to:
			const Setting &outfiles = base["file_path"]["outfiles"];

			string path1 = outpath+"/"+"state"+suffix;
		//	cout << path1 << endl;
			ofstream file1;
			file1.open(path1);
			//string header1 = "iter tree ysim1 ysim2 zhat1 zhat2 qstar\n";
			//file1 << header1;

			string path2 = outpath+"/"+"moves"+suffix;
			ofstream file2;
			file2.open(path2);
			//string header2 = "iter tree move accept\n";
			//file2 << header2;

			string path3 = outpath+"/"+"posterior"+suffix;
			ofstream file3;
			file3.open(path3);

			string path4 = outpath+"/"+"tree"+suffix;
			ofstream file4;
			file4.open(path4);

			// Get a line of data from .txt file and make predictor names
			MatrixXd datai(NSPLITS, DATA_LENGTH);
			datai.row(0) = read_line(0, DATA_LENGTH, ifs);
			vector<string> myNames = make_names(100, DATA_LENGTH);
			//cout << datai.row(0) << endl;

			// Initialise the forest
			//Forest forest(FOREST);
			Forest forest;
			Tree * tree = forest.fptr;	
			vector<int> leaf_nodes_up;
			vector<int> leaf_nodes_prop;
			vector<int> int_nodes;
			VectorXd datavec(DATA_LENGTH);

			VectorXd qvec(FOREST);
			VectorXd postvec(FOREST);
			double max_post;
			double qstar;
			double const_prop;
		
			// Use a k-loop to intialise forest trees
			for (int k = 0; k < FOREST; k++) {
				//cout << "k is: " << k << endl;
				Tree itree(1, DIMS);
				//itree.display(itree.root);
				forest.fvec.push_back(itree);
				//cout << &forest.fvec[k] << endl;
				tree = &forest.fvec[k];
				//cout << "Address of Tree is " << tree << endl;
				forest.init_tree(tree, AT_0, AI_0, ALPH, BET, D_0, D_0, 
					POST_INIT, LPOST_INIT, myNames, datai, MU, SIG);
			//	cout << "For Tree " << k << " lpost is: " << 
			//	tree->root->lpost <<" and ppost is:" << tree->root->ppost
			//	<< " and tree->tree_post is: " << tree->tree_post 
			//	<< endl;
			//	cout << "For Tree " << k << " aimat is: " << 
			//	tree->root->aimat <<" and atmat is:" << tree->root->atmat
			//	<< " and tree->tree_post is: " << tree->tree_post 
			//	<< endl;
			}
			//cout << "Initialisation over.\n"; 

			// Begin loop over trees
			for (int l = 1; l < LOOPS; l++) {
			//	cout << "ITERATION NUMBER IS: " << l << endl << endl;
				int i = l%NSPLITS;
				double num_splits = nsplit(1, NSPLITS);
				double nprule = prule(DATA_LENGTH, num_splits);
				file3 << l << " ";
				datai.row(i) = read_line(l*DATA_LENGTH, 
								DATA_LENGTH*(l+1), ifs);
				datavec = datai.row(i);
			//	cout << datavec.transpose() << endl;

				// Get y
				VectorXd yi(DIMS);
				//yi = sim_y(datavec, VT, A_mat, HT);
				yi = read_line(l*DIMS, DIMS*(l+1), yfs);
				//cout << yi.transpose() << endl;

				Node * found;
				for (int k = 0; k < FOREST; k++) {
					file1 << l << " " << k << " ";
					file2 << l << " " << k << " ";
				//	cout << "TREE NUMBER IS: " << k << endl << endl;

					tree = &forest.fvec[k];
					leaf_nodes_up = forest.lnodes(tree);
					leaf_nodes_prop = leaf_nodes_up;
					int nleaves = leaf_nodes_up.size();
					tree->find_int_nodes(tree->root, int_nodes);

					found = forest.ch_node(tree, datavec, myNames);
				//	cout << "Found ppost is: " << found->ppost << endl;
				//	cout << "Found lpost is: " << found->lpost << endl;

					forest.update(tree, found,
					leaf_nodes_up, DIMS, yi, WT, VT, WTI, 
					VTI, FT, HT, HVHT, VHT, WFT, WF, FWFT, 
					VTD, WTD);

					file1 << yi.transpose() << " ";
					file1 << found->state.transpose() << " ";

					// Copy tree
					Tree copy(1, DIMS);
					//cout << "init tree copy:\n";
					//copy.display(copy.root);
					copy = *tree;
					//cout << "tree copy of tree:\n";
					//copy.display(copy.root);
					Tree * cptr = &copy; 
				
					// Update tree prior
					tree->find_prior(tree->root);

					// Propose new tree
					forest.proposal(cptr, found, leaf_nodes_prop, int_nodes,
						ALPH, BET, nprule, DIMS, myNames, datai, i, WT, VT,
						WTI, VTI, FT, HT, HVHT, VHT, WFT, WF,
						FWFT, yi, VTD, WTD);
					file2 << cptr->move << " ";
					cout << cptr->move << " ";
					int newtree;
					newtree = forest.ch_tree(tree, cptr);
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
					forest.treeout(tree->root, file4);
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
			//	cout << " " << l << " " << flush;
			}
			cout << endl;
		}
	}
}
