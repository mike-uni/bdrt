// A SIMULATED TREE

# pragma once
# include <iostream>
# include <vector>
# include <string>
# include <fstream>
# include <Eigen/Dense>
# include <iomanip>
# include <cstdlib>
# include <libconfig.h++>
# include "read_data.h"
# include "stree_class.h"
# include "snode_class.h"
# include "settings_class.h"

using namespace std;
using namespace libconfig;

class Simtree {
	public:
	Simtree(){};
	void simcfg(Config& simcfg, Settings& simset, Dist& dist, string simtype,
				string cfgname);
};

void Simtree:: simcfg(Config& simcfg, Settings& simset, Dist& dist, string simtype,
				string cfgname) {	
/*	cout << "Welcome to the data simulation for BDRT. Please enter\n"
"the name of the configuration file (without extension):\n";
	string cfgname;
	while(cin.peek() != '\n'){
		cin >> cfgname;
		if (cin.peek() == '\b') cin.ignore();
	}
	cin.ignore();
	cfgname = cfgname+".cfg";
	cout << cfgname << endl;
	// Read the file. If there is an error, report it and exit;
*/	
	cout << "The simulation has begun: " << endl;
	cout << "The config file is: " << cfgname << endl;
	try {
		simcfg.readFile(cfgname.c_str());
	}
	catch(const FileIOException &fioex) {
		cerr << "I/O error while reading file." << endl;
	}
	catch(const ParseException &pex) {
		cerr << "Parse error at " << pex.getFile() << ":"
			<< pex.getLine() << " - " << pex.getError() << endl;
	}
			
	// Set config base
	const Setting& base = simcfg.getRoot();
			
	// Input path
	Setting& filepath = base["file_path"];
	string outpath = filepath["outpath"];
	string infolder = filepath["infolder"];
	string outfolder = filepath["outfolder"];
 
	const Setting& basicsetup = base["basic_setup"];
	int ndset = basicsetup.getLength();
	
	int DIMS = basicsetup[0][0];	

	vector<int> iloops;
	int LOOPS;
	for (int j=0; j<ndset; j++) {
		const Setting& lsets = basicsetup[j];
		iloops.push_back(lsets[1]);
	}
	LOOPS = *max_element(iloops.begin(), iloops.end());

	int DL2 = 0;
	int DATA_LENGTH;
	bool keepfollowing = false;
	for (int i = 0; i < ndset; i++) {
			const Setting& sets = basicsetup[i];
			DATA_LENGTH = sets[3];
			
			//cout << "DL1 = " << DATA_LENGTH << " DL2 = " << DL2 << endl;
			if (DATA_LENGTH == 0) {
				break;
			}
			else if (DL2 != DATA_LENGTH) {
				DL2 = DATA_LENGTH;

				MatrixXd WT(DIMS, DIMS);
				MatrixXd VT(DIMS, DIMS);
				MatrixXd HT(DIMS, DIMS);
				MatrixXd FT(DIMS, DIMS);
				VectorXd MU(DIMS);

				const Setting& matrix_list = base["matrix_setup"];
				int matset = matrix_list.getLength();
				double mymatarr[DIMS*DIMS];
				double myvecarr[DIMS];
				
				for (int k = 0; k<matset; k++) {

					const Setting& matrix_setup = matrix_list[k];

					// Matrix Variables
					simset.getdArray(mymatarr, matrix_setup[1], 0, DIMS*DIMS);
					WT = Map<MatrixXd>(mymatarr, DIMS, DIMS);
					cout << WT << endl;

					simset.getdArray(mymatarr, matrix_setup[2], 0, DIMS*DIMS);
					VT = Map<MatrixXd>(mymatarr, DIMS, DIMS);
					cout << VT << endl;

					simset.getdArray(mymatarr, matrix_setup[3], 0, DIMS*DIMS);
					HT = Map<MatrixXd>(mymatarr, DIMS, DIMS);
					cout << HT << endl;

					simset.getdArray(mymatarr, matrix_setup[4], 0, DIMS*DIMS);
					FT = Map<MatrixXd>(mymatarr, DIMS, DIMS);
					cout << FT << endl;

					simset.getdArray(myvecarr, matrix_setup[6], 0, DIMS);
					MU = Map<VectorXd>(myvecarr, DIMS);
					cout << MU << endl;

					vector<string> myNames = make_names(100, DATA_LENGTH);
					vector<int> nodevec;
					string mypred;
					double mysplit;
					string mypredl;
					double mysplitl;
					string mypredr;
					double mysplitr;

					sTree sim(1, DIMS);
					nodevec = {1,2,3,4,6,12};
					int num_nodes = nodevec.size();
					sNode * temp;

					if (keepfollowing == false && simtype == "mod") {
						filepath.remove("infilesx");
						filepath.remove("infilesy");
						filepath.remove("infilesz");
						filepath.remove("infilest");
						filepath.add("infilesx", Setting::TypeList);
						filepath.add("infilesy", Setting::TypeList);
						filepath.add("infilesz", Setting::TypeList);
						filepath.add("infilest", Setting::TypeList);
						keepfollowing = true;
					}
					else { 
						keepfollowing = true;
					}

					Setting &xfiles = filepath["infilesx"];
					//int xlen = xfiles.getLength();
					//string x1 = xfiles[0];
					//cout << x1 << endl;
					Setting &yfiles = filepath["infilesy"];
					//int ylen = yfiles.getLength();
					//string y1 = yfiles[0];
					//cout << y1 << endl;
					Setting &zfiles = filepath["infilesz"];
					//int zlen = zfiles.getLength();
					//string z1 = zfiles[0];
					//cout << z1 << endl;
					Setting &tfiles = filepath["infilest"];
					//int tlen = tfiles.getLength();
					//string t1 = tfiles[0];
					//cout << t1 << endl;
					
					string filesuffix = "_"+to_string(DIMS)
							+to_string(LOOPS)+to_string(DATA_LENGTH)
							+to_string(k);
					cout << filesuffix << endl;
				
					string path4=outpath+"/"+infolder+"/"+"tree"+
									filesuffix+".txt";
					string s1 = "tree"+filesuffix+".txt";
					tfiles.add(Setting::TypeString) = s1;
					ofstream file4;
					file4.open(path4);

					for (int i = 0; i < DATA_LENGTH; i++) {
						if (sim.root->split == 0) {
							int rand0 = dist.runif(1, DATA_LENGTH);
							mypred = ("X" + to_string(rand0));
							mysplit = dist.runif(1, rand0*3);
							sim.root->pred = mypred;
							sim.root->split = mysplit;
						}
						else {
							file4 << sim.root->nnode << " " << mypred << 
								" " << mysplit << endl;
							break;
						}
					}	

					for (int i = 0; i < num_nodes; i++) {
						int nnum = nodevec[i];
						int rand1 = dist.runif(1, DATA_LENGTH);
						mypredl = ("X" + to_string(rand1));
						mysplitl = dist.runif(1, rand1*3);
						int rand2 = dist.runif(1, DATA_LENGTH);
						mypredr = ("X" + to_string(rand()%DATA_LENGTH+1));
						mysplitr = dist.runif(1, rand2*3);
						temp = sim.find_node(sim.root, nnum);
						sim.sgrow(temp, DIMS, mypredl, mysplitl, 
									mypredr, mysplitr);
						file4 << 2*temp->nnode << " " << mypredl << 
							" " << mysplitl << " " << 2*temp->nnode+1 
							<< " " << mypredr << " " << mysplitr 
							<< endl;
					}

					//sim.display(sim.root);
				
					// Open files to write to:
					string path1=outpath+"/"+infolder+"/ysim"+
									filesuffix+".txt";
					yfiles.add(Setting::TypeString)="ysim"+filesuffix+".txt";
					ofstream file1;
					file1.open(path1);

					string path2=outpath+"/"+infolder+"/zsim"+
									filesuffix+".txt";
					zfiles.add(Setting::TypeString)="zsim"+filesuffix+".txt";
					ofstream file2;
					file2.open(path2);

					string path3 = outpath+"/"+outfolder+"/xline.txt";

					string path5 = outpath+"/"+infolder+"/xsim"+
									filesuffix+".txt";
					xfiles.add(Setting::TypeString)="xsim"+filesuffix+".txt";
					ofstream file5;
					file5.open(path5);

					VectorXd inData(DATA_LENGTH);
					sNode * leaf;

					for (int j = 0; j < LOOPS; j++) {
						sim_x(DATA_LENGTH, path3, 1);
						ifstream ifs2(path3, ios::in);
						inData = read_line(0, DATA_LENGTH, ifs2);
						file5 << inData.transpose() << endl;
						leaf = sim.find_leaf(sim.root, inData, myNames);
						//cout << leaf->nnode << " ";
						if (leaf->zprev.isZero()) {
							leaf->zprev = MU;
							leaf->wvar = WT;
							leaf->vvar = VT;
						}
						leaf->zcur = dist.mvn(1, FT*leaf->zprev, leaf->wvar);
						//cout << "leaf->zcur is:\n" << leaf->zcur << endl;
						leaf->zprev = leaf->zcur;
						leaf->ycur = dist.mvn(1, HT*leaf->zcur, leaf->vvar);
						//cout << "leaf->ycur is:\n" << leaf->ycur << endl;
						file1 << leaf->ycur.transpose() << endl;
						file2 << leaf->zcur.transpose() << endl;
					}
				}
			}
			else {
				continue;
			}
		}

  	// Write out the updated configuration.
  	try {
		simcfg.writeFile(cfgname.c_str());
		cerr << "Updated config written to: " << cfgname << endl;
  	}
  	catch(const FileIOException &fioex) {
		cerr << "I/O error while writing file: "<< cfgname << endl;
  	}
	cout << "The simulation is over." << endl;
}
/*
int main() {
	Config cfg;
	Settings set;
	Dist d;
	Simtree simtree;
	simtree.simcfg(cfg, set, d);
}*/
