// THIS IS THE MAIN FILE FOR BAYESIAN DYNAMIC REGRESSION TREES

# include <iostream>
# include <Eigen/Dense>
# include <fstream>
# include <algorithm>
# include <vector>
# include <time.h>
# include <chrono>
# include <dirent.h>
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
# include "simfixed_class.h"

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
	string pwd = cwd.substr(0, cwd.size()-3);
	set.outpath = pwd;

	// Declare path names
	string outpath;
	string inpath;

	// Show .cfg files in cwd	
	get_exts(cwd, ".cfg");

	cout << "If you wish to use one of these config files please enter the name, \n"
		"if you wish to modify an existing one, please type 'mod'\n"
		"or if you wish to create a new one, please type 'new'\n";
	string exname;
	getline(cin, exname);
	if (exname == "new") {
		// Choose type of experiment
		cout << "As you have selected to create a new config file you will \n"
			"need to choose the type of experiment you wish to conduct.\n"
			"Please enter either either 1 or 2 corresponding to the options\n"
			"described in the introduction and press enter.\n";
		int opt = set.get_int();

		// Get input folder name for this run of experiments
		cout << "Please enter the name of the separate folders in which\n" 
"the predictors and responses are currently being stored. " 
"If they do not yet exist, they will be created.\n";
		cout << "First the input folder:\n";
		string infoldername;
		getline(cin, infoldername);
		
		set.infolder = infoldername;

		// Make new folder path 
		inpath = pwd+infoldername;
		cout << inpath << endl;

		// Make new folder
		set.make_folder(inpath);

		// Get folder name for output
		cout << "Now the output folder.\n";
		string foutname;
		getline(cin, foutname);
		set.outfolder = foutname;

		// Make new folder path
		outpath = pwd+foutname;
		cout << outpath << endl;
	
		// Make new folder
		set.make_folder(outpath);

		// Set the config file name
		cout << "Please enter the name of the new .cfg file (without extension): "; 
		string newname;
		getline(cin, newname);
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
			//cin.ignore();
		}
		else if (opt == 2) {
			// Set the number of matrix examples
			cout << "Please enter the number of matrix examples: \n";
			int nmats = set.get_int();
			set.nmat = nmats;

			// Show the example output file
			set.examplemat();
			cout << endl;

			// Make the cfg file
			set.make_file(cfg);
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
		getline(cin, modname);
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

		Setting& fileinfo = modbase["file_path"];
		fileinfo.remove("outpath");
		fileinfo.add("outpath", Setting::TypeString) = set.outpath;
		
		cout <<"Please enter the setting you would like to modify:\n";
		string modset;
		bool right = true;
		while (right) {
			getline(cin, modset);
			if (modset == "basic_setup") {
				try {
					modbase.remove("basic_setup");
					modbase.add("basic_setup", Setting::TypeList);
					
					Setting& curinfo = modbase["setting_info"];
					set.nmat = curinfo["nmat"];
				
					modbase.remove("setting_info");
					Setting& setinfo = modbase.add("setting_info", 
										Setting::TypeGroup);
				
					Setting &upmod = modbase[modset];
					set.bsetup(set.nmat);

					setinfo.add("dim", Setting::TypeInt) = set.dim;
					setinfo.add("nloops", Setting::TypeInt) = set.nloops;
					setinfo.add("nforests", Setting::TypeInt) = set.nforest;
					setinfo.add("nrest", Setting::TypeInt) = set.nrest;
					setinfo.add("nmat", Setting::TypeInt) = set.nmat;

					set.bset_fill(upmod);

					right = false;
				}
				catch(const SettingNotFoundException &nfex) {
					cerr << "Setting not found. Please try again."<< endl;
				}
			}
			else if (modset == "matrix_setup") {
				try{
					modbase.remove("matrix_setup");
					modbase.add("matrix_setup", Setting::TypeList);

					cout << "Please enter the number of matrix examples: \n";
					int nmats = set.get_int();
					set.nmat = nmats;
					
					Setting& curinfo = modbase["setting_info"];
					set.dim = curinfo["dim"];

					Setting &upmod = modbase[modset];
					set.mat_setup(upmod, set.nmat);

					right = false;
				}
				catch(const SettingNotFoundException &nfex) {
					cerr << "Setting not found. Please try again."<< endl;
				}
			}
			else {
				cout << "This setting is unavailable. Please choose another: ";		
				getline(cin, modset);
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
		inpath = cwd+set.infolder;

		// Make new output folder path 
		string outfolder = filemod.lookup("outfolder");
		set.outfolder = outfolder;
		outpath = cwd+set.outfolder;

		set.outname = modname;
		// Re-simulate Data
		cout << "You will now simulate data based on the configuration\n"
			"you have modified. This will be used in the BDRT model below.\n";
		//cin.ignore();
		simtree.simcfg(cfg, set, d, exname, set.outname);
	}
	else {
		set.outname = exname+".cfg";
		cout << "Would you like to re-simulate the data? (yes/no)" << endl;
		string resim;
		getline(cin, resim);
		if (resim == "yes") {
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
			// Get the base of the configuration file
			const Setting& base = cfg.getRoot();

			// Get paths
			Setting& fileinfo = base["file_path"];
			fileinfo.remove("infilesx");
			fileinfo.remove("infilesy");
			fileinfo.remove("infilesz");
			fileinfo.remove("infilest");
			fileinfo.add("infilesx", Setting::TypeList);
			fileinfo.add("infilesy", Setting::TypeList);
			fileinfo.add("infilesz", Setting::TypeList);
			fileinfo.add("infilest", Setting::TypeList);
			fileinfo.remove("outpath");
			fileinfo.add("outpath", Setting::TypeString) = pwd;

			// Write out the new configuration.
			const char * reout = set.outname.c_str();
			try {
				cfg.writeFile(reout);
				//cerr << "New config file written to: " << modname << endl;
			}
			catch (const FileIOException &fioex) {
				cerr << "I/O error while writing file: " << set.outname << endl;
			}
			// Delete files in folders
			// Get paths
			string outfolder = fileinfo.lookup("outfolder");
			string infolder = fileinfo.lookup("infolder");
			inpath = pwd+infolder;
			outpath = pwd+outfolder;

			DIR * dir = opendir(outpath.c_str());
			if (dir) { closedir(dir); }
			else if (ENOENT == errno) {
				set.make_folder(outpath);
				set.make_folder(inpath);
			}

			del_files(inpath);
			del_files(outpath);

			// Resimulate data
			simtree.simcfg(cfg, set, d, exname, set.outname);
		}
		//cin.ignore();
	}
	cout << "The current config file is "<< set.outname <<"." << endl;
	cout << "The config is now complete. To continue onto the BDRT please type 'yes', else 'exit' will exit the program:\n";

	string cont;
	while (cin.peek() != '\n') {
		getline(cin, cont);
		if (cont == "yes") break;
		else return(EXIT_SUCCESS);
	}
	//cin.ignore();		
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
	inpath = pwd+infolder;
	outpath = pwd+outfolder;

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

		// Get iterations of data variable or number of matrices
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
			MatrixXd W0I(DIMS, DIMS);
			W0I = W0.inverse();

			MatrixXd FWIFT(DIMS, DIMS);
			FWIFT = FT.transpose()*WTI*FT;
			MatrixXd HVIHT(DIMS, DIMS);
			HVIHT = HT.transpose()*VTI*HT;
			MatrixXd VIHT(DIMS, DIMS);
			VIHT = VTI*HT;
			MatrixXd WIFT(DIMS, DIMS);
			WIFT = (WTI*FT).transpose();
			MatrixXd WIF(DIMS, DIMS);
			WIF = WTI*FT;

			// Get intial values for posterior conditional
			MatrixXd AT_0(DIMS, DIMS);
			AT_0 = W0.inverse()+FWIFT;
			MatrixXd AI_0(DIMS, DIMS);
			AI_0.setZero(DIMS, DIMS);
			VectorXd D_0(DIMS);
			D_0 = MU.transpose()*W0;

			double VTD = VT.determinant();
			double WTD = WT.determinant();
			double LG_W0D = log(W0.determinant());
			double LG_AD = log(AT_0.determinant());
			double W0_SQ = (MU.transpose()*W0I*MU);
			double ADSQ_0 = (D_0.transpose()*AT_0.inverse()*D_0);
			double POST_INIT = -0.5*(LG_W0D - LG_AD + W0_SQ - (0.25)*ADSQ_0);
			double LPOST_INIT = 0.0;

			cout << "LG_W0D: " << LG_W0D << " W0_SQ: " << W0_SQ 
				 << " LG_AD: " << LG_AD << " ADSQ_0: " << ADSQ_0 
				 << " LPOST_INIT: " << LPOST_INIT << " POST_INIT: " << POST_INIT
				<< endl;

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
							"_"+to_string(ALPH)+"_"+to_string(BET)+"_"+to_string(rn)+
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

			string path5 = outpath+"/"+"time"+suffix;
			ofstream file5;
			file5.open(path5);

			// Get a line of data from .txt file and make predictor names
			MatrixXd datai(NSPLITS, DATA_LENGTH);
			datai.row(0) = read_line(0, DATA_LENGTH, ifs);
			vector<string> myNames = make_names(100, DATA_LENGTH);
			//cout << datai.row(0) << endl;

			// Initialise the forest
			//Forest forest(FOREST);
			Forest forest;
			Tree * tree = forest.fptr;	
			vector<int> leaf_nodes;
			vector<int> int_nodes;
			vector<int> leaf_nodes_prop;
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
			}

			for (int k = 0; k < FOREST; k++) {
				tree = &forest.fvec[k];
				//cout << "Address of Tree is " << tree << endl;
				forest.init_tree(tree, AT_0, AI_0, ALPH, BET, D_0, D_0, 
					POST_INIT, LPOST_INIT, k, myNames, datai, MU, SIG);
				cout << "For Tree " << tree->tname << " lpost is: " << 
				tree->root->lpost <<" and ppost is:" << tree->root->ppost
				<< " and tree->tree_post is: " << tree->tree_post 
				<< endl;
				cout << "For Tree " << k << " aimat is: " << 
				tree->root->aimat <<" and atmat is:" << tree->root->atmat
				<< " and tree->tree_post is: " << tree->tree_post 
				<< endl;
			}
			//cout << "Initialisation over.\n";
			cout << "The names of initalised trees" << endl;	 
			for (int k = 0; k < FOREST; k++) {
				cout << forest.fvec[k].tname << " ";
				cout << forest.fvec[k].root->ppost << " ";
				cout << forest.fvec[k].tree_post << " ";
			}
			cout << endl;

			Tree * ltree = forest.fptr;	

			// Begin loop over trees
			for (int l = 1; l < LOOPS; l++) {
				chrono::high_resolution_clock::time_point init1, init2, end2, end1;
				init1 = chrono::high_resolution_clock::now();
				file5 << l << " ";
				cout << "ITERATION NUMBER IS: " << l << endl << endl;
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
				Node * cfound;
				for (int k = 0; k < FOREST; k++) {
					file1 << l << " " << k << " ";
					file2 << l << " " << k << " ";
					file5 << k << " ";
					init2 = chrono::high_resolution_clock::now();
				//	cout << "TREE NUMBER IS: " << k << endl << endl;

					ltree = &forest.fvec[k];
					cout << "TREE NUMBER IS: " << ltree->tname << endl;
					cout << "The initial tree display is:\n";
					ltree->display(ltree->root);
					cout << endl;
					cout << "&forest.fvec[k] " << &forest.fvec[k] << endl;

					leaf_nodes = forest.lnodes(ltree);
					int nleaves = leaf_nodes.size();
					ltree->find_int_nodes(ltree->root, int_nodes);
					cout << "The int_nodes are:\n";
					for (i = 0; i < int_nodes.size(); i++) {
						cout << int_nodes[i] << " ";
					}
					cout << endl;

					//STAGE 1. CHOOSE A LEAF USING CURRENT STATE SPACE
					// AND CURRENT DATA.
					found = forest.ch_node(ltree, datavec, myNames);
					cout << "Found nnode is: " << found->nnode << endl;
					cout << "Found ppost is: " << found->ppost << endl;
					cout << "Found lpost is: " << found->lpost << endl;

					forest.kupdate(ltree, found, leaf_nodes, DIMS, yi, 
									WT, VT, FT, HT);

					file1 << yi.transpose() << " ";
					file1 << found->state.transpose() << " ";

					// COPY THE CURRENT TREE
					Tree copy(1, DIMS);
					//cout << "init ltree copy:\n";
					//copy.display(copy.root);
					copy = *ltree;
					//cout << "ltree copy of ltree:\n";
					//copy.display(copy.root);
					Tree * cptr = &copy; 

					// STAGE 2. CHOOSE A RANDOM LEAF FROM THE CURRENT TREE
					// AND PROPOSE A MOVE.
					forest.proposal(cptr, found, leaf_nodes, int_nodes,
						ALPH, BET, nprule, DIMS, myNames, datai, i);
					file2 << cptr->move << " ";
					cout << "The ltree move is: " << cptr->move << "\n";
					cout << "The int_nodes are:\n";
					for (i = 0; i < int_nodes.size(); i++) {
						cout << int_nodes[i] << " ";
					}
					cout << endl;

					// STAGE 3. RESELECT LEAF USING CURRENT DATA AND PROPOSED TREE
					// AND UPDATE POSTERIOR DATA
					cfound = forest.ch_node(cptr, datavec, myNames);

					// 3a. Update current tree (ltree).
					cout << "Update current tree:\n";
					forest.lupdate(ltree, found, leaf_nodes, yi, WTI, VTI, HVIHT,
									VIHT, WIFT, WIF, FWIFT, VTD, WTD);

					// Update ltree prior
					cout << "Update current tree Prior:\n";
					ltree->find_prior(ltree->root);
					cout << ltree->prior << endl << endl;

					// Update ltree tree_post
					cout << "Update current tree Posterior:\n";
					ltree->find_post(ltree->root);
					cout << ltree->tree_post << endl << endl;

					// 3b. Update proposed tree (ltree).
					cout << "Update copy tree:\n";
					leaf_nodes_prop = forest.lnodes(cptr);
					forest.lupdate(cptr, cfound, leaf_nodes_prop, yi, WTI, VTI, HVIHT,
									VIHT, WIFT, WIF, FWIFT, VTD, WTD);

					// Update copy prior
					cout << "Update copy tree Prior:\n";
					cptr->find_prior(cptr->root);
					cout << cptr->prior << endl << endl;

					// Update copy tree_post
					cout << "Update copy tree Poserior:\n";
					cptr->find_post(cptr->root);
					cout << cptr->tree_post << endl << endl;

					// STAGE 4. CHOOSE TREE AND CALCULATE POSTERIORS
					cout << "Choose Tree:\n";
					int newtree;
					newtree = forest.ch_tree(ltree, cptr);
					cout << newtree << endl << endl;
					if (newtree == 2) {
						file2 << 1 << "\n";
						*ltree = copy;
					//	copy.~Tree();
					} 
					else { 
						file2 << 0 << "\n";
					//	copy.~Tree();
					}
					file4 << l << " " << k+1 << "   ";
					forest.treeout(ltree->root, file4);
					file4 << endl;
					cout << "Display Tree " << ltree->tname << " :\n";
					ltree->display(ltree->root);
					cout << endl;

					double qstar;
					qstar = log(ltree->prior) + ltree->tree_post;
					cout << "log(ltree->prior) is: " << log(ltree->prior) << " and " <<
						    "ltree->ltree_post is: " << ltree->tree_post << endl;
					qvec(k) = qstar;
					if (isnan(qstar) || isinf(qstar)) {
						cout << "Tree_post is nan or inf." << endl << endl;
						break;
					}
					file1 << qstar << "\n";
					cout << "qstar is: " << qstar << "\n\n";
					ltree->prior = 1.0;
					ltree->tree_post = 0.0; // Updated in each iteration
					int_nodes.clear();
					end2 = chrono::high_resolution_clock::now();
					file5 << chrono::duration_cast<chrono::microseconds>(end2 - init2).count() << " ";
				}
				//cout << endl;
				max_post = qvec.maxCoeff();
				qvec = qvec.array() - max_post;
				qvec = qvec.array().exp();
				cout << "qvec is: " << qvec.transpose() << endl;
				const_prop = qvec.sum();
				cout << "const_prop is: " << const_prop << endl;
				if (isnan(const_prop) || isinf(const_prop)) {
					cout << "Tree_post is nan or inf." << endl << endl;
					break;
				}
				postvec = qvec.array()/const_prop;
				//file3 << tree_index.transpose() << "\n"; 
				file3 << postvec.transpose() << "\n";
				cout << "The ltree posteriors at " << l << " are:\n" <<
						postvec.transpose() << endl << endl;
				//cout << " " << l << " " << flush;
				end1 = chrono::high_resolution_clock::now();
				file5 << chrono::duration_cast<chrono::microseconds>(end1 - init1).count() << "\n";
			}
			cout << endl;
		}
	}
}
