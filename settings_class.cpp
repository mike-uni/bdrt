// GENERATE CFG FILE.

# pragma once
# include <iostream>
# include <fstream>
# include <string>
# include <vector>
# include <unistd.h>
# include <sys/stat.h>
# include <sys/types.h>
# include <iomanip>
# include <cstdlib>
# include <libconfig.h++>

using namespace std;
using namespace libconfig;

class Settings {
	public:
	string outname;
	string outpath;
	string outfolder;
	string infolder;
	int dim;
	int nloops;
	int nforest;
	int nrest;
	int nmat;
	vector <int> loops;
	vector <int> forest;
	vector <int> dl;
	vector <int> ns;
	vector <float> alpha;
	vector <float> beta;

	Settings(){};
	//Settings(int dims, int nloops, int nforest, int n) : 
	//		dim(dims), loops(nloops, 0), forest(nforest, 0), 
	//		dl(n, 0), ns(n, 0), alpha(n, 0), beta(n, 0) {};

	void get_wd(); 
	void example();
	void examplemat();
	float get_float();
	int get_int();
	void add_mat(Setting& set);
	void add_vec(Setting& set);
	vector<int> get_vint(vector<int>& vec);
	vector<float> get_vfloat(vector<float>& vec);
	int get_dim();
	void set_file();
	void make_file(Config& cfg);
	void make_folder(string foldername);
	void set_vint(vector<int>& vec, const int vsize);
	void set_vfloat(vector<float>& vec, const int vsize);
	void bsetup(int nmats);
	void mat_setup(Setting& matlist, int nmats);
	void bset_fill(Setting& basiclist);
	void getintArray(int marr[], const Setting& mat, 
						int start, int end);
	void getdArray(double marr[], const Setting& mat, 
						int start, int end);
//	void make_matfile();
};

void Settings :: get_wd() {
	char nbuf[1024];
	outpath = getcwd(nbuf, 1024); 
}

void Settings :: make_folder(string foldername) {
	mode_t mmode = 0777;
	int ret = mkdir(foldername.c_str(), mmode);
}

void Settings :: example() {
	cout << "This is an example of the top few lines of a basic\n" 
"settings file in which there is only one matrix setup.\n"
"Note that there is no header so the entries,\n"
"from left to right, are:\n"
"DIM, LOOPS(N), FOREST, DATA_LENGTH(P), NSPLITS(DEPTH), ALPHA, BETA.\n";
	cout << "\t1 1000 10 10 10 0.95 1.0 \n";
	cout << "\t1 1000 50 10 10 0.95 1.0 \n";
	cout << "\t1 1000 100 10 10 0.95 1.0 \n";
	cout << "\t1 1000 500 10 10 0.95 1.0 \n";
	cout << "\t...\n";
	cout << "\t1 2000 10 10 10 0.95 1.5 \n";
	cout << "\t1 2000 50 10 10 0.95 1.5 \n";
	cout << "\t ...\n";
	cout << "\t1 5000 100 100 50 0.95 2.0 \n";
	cout << "\t1 5000 500 100 50 0.95 2.0 \n\n";
	cout << "NOTE: DIM is fixed for a particular settings file.\n"
"Numbers of LOOPS and FOREST can be any positive integer.\n"
"Numbers of DATA_LENGTH, NSPLITS, ALPHA and BETA examples will need\n"
"to be the same. For instance:\n"
"BETA = {1.0, 1.25, 1.5, 1.75, 2.0} then\n"
"ALPHA = {0.85, 0.85, 0.95, 0.95, 0.95}\n"
"and similarly for the others.\n";
}

void Settings :: examplemat() {
	cout << "This is an example of the top few lines of a file with\n"
"different matrices but with DATA_LENGTH, NSPLITS, ALPHA and BETA being\n"
"constant. Note that there is no header so the entries,\n"
"from left to right, are:\n"
"DIM, LOOPS(N), FOREST, W0, W, V, H, F, SIG, MU\n";
	cout << "\t1 1000 10 0.99 1 1 1 1 1 0  \n";
	cout << "\t1 1000 10 0.99 1 2 1 1 1 0  \n";
	cout << "\t1 1000 10 0.99 1 4 1 1 1 0  \n";
	cout << "\t...\n";
	cout << "\t1 1000 10 0.99 1 1 0.5 1 0  \n";
	cout << "\t1 1000 10 0.99 1 1 0.75 1 0  \n";
	cout << "\t1 1000 10 0.99 1 1 0.95 1 0  \n";
	cout << "\t ...\n";
	cout << "\t1 1000 10 0.99 1 2 0.75 1 0  \n";
	cout << "\t1 1000 10 0.99 1 2 0.75 1 0  \n";
	cout << "\t1 1000 10 0.99 1 2 0.75 1 0  \n";
	cout << "NOTE: DIM is fixed for a particular settings file.\n"
"Numbers of LOOPS and FOREST can be any positive integer.\n"
"The entries for the matrices can be any double real. W0 CANNOT be the\n"
"same as W and F is typically < 1 for 1d or if d > 1 then the eigenvalues\n"
"of F must have negative real values.\n. Also note that each heading will\n"
"not necessarily correspond to a column. For instance, if DIM = 2 then:\n"
"DIM, LOOPS(N), FOREST, W0, W, V, H, F, SIG, MU\n";
	cout << "\t1 1000 10 0.99 0 0 0.99 1 0 0 1 2 0 0 2 1 0 0 1 -0.5 0 0 -0.5 1 0 0 1 0 0\n";
	cout << "\t1 1000 100 0.99 0 0 0.99 1 0 0 1 2 0 0 2 1 0 0 1 -0.5 0 0 -0.5 1 0 0 1 0 0\n";
	cout << "\t1 1000 150 0.99 0 0 0.99 1 0 0 1 2 0 0 2 1 0 0 1 -0.5 0 0 -0.5 1 0 0 1 0 0\n";
	cout << "\t...\n";
}

vector<int> Settings :: get_vint(vector<int>& vec) {
	return vec;
}
vector<float> Settings :: get_vfloat(vector<float>& vec) {
	return vec;
}

int Settings :: get_dim() {
	return dim;
}

float Settings :: get_float() {
	string input;
	float entry;
	while (true) {
		cout << "Please input a float value: ";
		getline(cin, input);
		stringstream matstream(input);
		if (matstream >> entry) break;
		cout << "Invalid entry, please try again." << endl;
	}
	return entry;
}

int Settings :: get_int() {
	string input;
	int entry;
	while (true) {
		cout << "Please input an int value: ";
		getline(cin, input);
		stringstream matstream(input);
		if (matstream >> entry) break;
		cout << "Invalid entry, please try again." << endl;
	}
	return entry;
}

void Settings :: add_mat(Setting& set) {
	//int dim = get_dim();
	float entry;
	for (int i = 0; i < dim*dim; i++) {
		entry = get_float();
		set.add(Setting::TypeFloat) = entry;
	}
}

void Settings :: add_vec(Setting& set) {
	//int dim = get_dim();
	float entry;
	for (int i = 0; i < dim; i++) {
		entry = get_float();
		set.add(Setting::TypeFloat) = entry;
	}
}

void Settings :: set_vint(vector<int>& vec, const int vsize) {
	vector<int> entries;
	vector<int> copies;
	
 	int unique;
	int entry;
	int copy;

	cout << "There is/are " << vsize << " maximum entries for this variable. \n";
	cout << "Please enter the unique values of the parameters. e.g. 100 200 300.\n"
			"If there are multiple copies of the same value you will be asked\n"
			"to enter these in the same order as previously. e.g. 5 2 1.\n"
			"The sum of the number of copies should equal " << vsize << ".\n\n";

	while (true) {
		if (vsize > 1) { 
			cout << "First, the enter number of unique values.\n";
			unique = get_int();
		}
		else { unique = 1; }
		cout << endl;

		cout << "Now input the values when requested.\n";
		for (int i = 0; i < unique; i++) {
			entry = get_int();
			entries.push_back(entry);
		}
		cout << endl;

		if (vsize > 1 && unique != 1 || vsize > 1 && unique != vsize) {
			cout << "Now the number of copies of each of the unique values:\n";
			for (int i = 0; i < unique; i++) {
				copy = get_int();
				copies.push_back(copy);
			}
		}
		else { copies.push_back((int)1); }
		cout << endl;
	
		int sumcopies = accumulate(copies.begin(), copies.end(), 0);	
		if (sumcopies == vsize) {break;}
		else {
		cout << "Number of copies is: " << copies.size() << endl << 
				"Number of entries is: " << entries.size() << endl <<
				"These should be equal." << endl;
		cout << "There seems to be a mistake, please try again.\n";
		copies.clear();
		entries.clear();
		}
	}

	for (int i = 0; i < copies.size(); i++) {
		int numcopy = copies[i];
		int nument = entries[i];
		for (int j = 0; j < numcopy; j++) {
			vec.push_back(nument);
		}
	}

	for (int i = 0; i < vec.size(); i++) {
		cout << "vec entries are " << vec[i] << endl;
	}
}
 

void Settings :: set_vfloat(vector<float>& vec, const int vsize) {
	vector<float> entries;
	vector<int> copies;
	
 	int unique;
	float entry;
	int copy;

	cout << "There is/are " << vsize << " maximum entries for this variable. \n";
	cout << "Please enter the unique values of the parameters. e.g. 0.75 0.85 0.95.\n"
			"If there are multiple copies of the same value you will be asked\n"
			"to enter these in the same order as previously. e.g. 5 2 1.\n"
			"The sum of the number of copies should equal " << vsize << "." << endl;

	while (true) {
		if (vsize > 1) { 
			cout << "First, the enter number of unique values.\n";
			unique = get_int();
		}
		else { unique = 1; }
		cout << endl;

		cout << "Now input the values when requested.\n";
		for (int i = 0; i < unique; i++) {
			entry = get_float();
			entries.push_back(entry);
		}
		cout << endl;

		if (vsize > 1 && unique != 1 || vsize > 1 && unique != vsize) {
			cout << "Now the number of copies of each of the unique values:\n";
			for (int i = 0; i < unique; i++) {
				copy = get_int();
				copies.push_back(copy);
			}
		}
		else { copies.push_back((int)1); }
		cout << endl;
		
		int sumcopies = accumulate(copies.begin(), copies.end(), 0);	
		if (sumcopies == vsize) {break;}
		else {
		cout << "Number of copies is: " << copies.size() << endl << 
				"Number of entries is: " << entries.size() << endl <<
				"These should be equal." << endl;
		cout << "There seems to be a mistake, please try again.\n";
		copies.clear();
		entries.clear();
		}
	}

	for (int i = 0; i < copies.size(); i++) {
		int numcopy = copies[i];
		float nument = entries[i];
		for (int j = 0; j < numcopy; j++) {
			vec.push_back(nument);
		}
	}

	for (int i = 0; i < vec.size(); i++) {
		cout << "vec entries are " << vec[i] << endl;
	}
}

void Settings :: make_file(Config& cfg) {
		
	Setting &base = cfg.getRoot();

	// Add some settings to cfg file
	Setting &model_info = base.add("model_info", Setting::TypeString) = 
	"Welcome to Baysian Dynamic Regression Trees.\n";

	Setting &file_path = base.add("file_path", Setting::TypeGroup);
	file_path.add("outpath", Setting::TypeString) = outpath;
	file_path.add("outfolder", Setting::TypeString) = outfolder;
	file_path.add("outname", Setting::TypeString) = outname;
	file_path.add("infolder", Setting::TypeString) = infolder;
	file_path.add("infilesx", Setting::TypeList);
	file_path.add("infilesy", Setting::TypeList);
	file_path.add("infilesz", Setting::TypeList);
	file_path.add("infilest", Setting::TypeList);
	file_path.add("outfiles", Setting::TypeList);

	base.add("basic_setup", Setting::TypeList);
	Setting &basic_setup = base["basic_setup"];
	bsetup(nmat);

	Setting &setinfo = base.add("setting_info", Setting::TypeGroup);
	setinfo.add("dim", Setting::TypeInt) = dim;
	setinfo.add("nloops", Setting::TypeInt) = nloops;
	setinfo.add("nforests", Setting::TypeInt) = nforest;
	setinfo.add("nrest", Setting::TypeInt) = nrest;
	setinfo.add("nmat", Setting::TypeInt) = nmat;

	bset_fill(basic_setup);
	cout << endl;

	Setting &matrix_setup = base.add("matrix_setup", Setting::TypeList);
	cout << "We now set up the matrices. Please enter the values when prompted.\n";
	mat_setup(matrix_setup, nmat);

	Setting &setname = base.add("setting_name", Setting::TypeList);
//	setname.add(Setting::TypeString) = "file_path";
	setname.add(Setting::TypeString) = "basic_setup";
	setname.add(Setting::TypeString) = "matrix_setup";

	const char * out = outname.c_str();
	// Write out the new configuration.
	try {
		cfg.writeFile(out);
		cerr << "New config file written to: " << outname << "\nin " <<
		outpath << endl;
	}
	catch (const FileIOException &fioex) {
		cerr << "I/O error while writing file: " << outname << endl;
	}
} 

void Settings:: bset_fill(Setting& basiclist) {
	for (int k = 0; k < nrest; k++) {
		for (int i = 0; i < nloops; i++) {
			for (int j = 0;  j < nforest; j++) {
			Setting &mylist = basiclist.add(Setting::TypeList);
			mylist.add(Setting::TypeInt) = dim;
			mylist.add(Setting::TypeInt) = loops[i];
			mylist.add(Setting::TypeInt) = forest[j];
			mylist.add(Setting::TypeInt) = dl[k];
			mylist.add(Setting::TypeInt) = ns[k];
			mylist.add(Setting::TypeFloat) = alpha[k];
			mylist.add(Setting::TypeFloat) = beta[k];
			}
		}
	}
}
	
void Settings:: mat_setup(Setting& matlist, int nmats) {
	for (int i = 0; i < nmats; i++) {
		Setting &nmatlist = matlist.add(Setting::TypeList);
		cout << "The number of matrices to fill are: " << nmats << endl;
		cout << "There are 7 entries for each matrix.\n";
		
		Setting &matarr = nmatlist.add(Setting::TypeArray);
		cout << "W0:\n";
		add_mat(matarr);
		cout << endl;
				
		Setting &matarr1 = nmatlist.add(Setting::TypeArray);
		cout << "WT:\n";
		add_mat(matarr1);
		cout << endl;
				
		Setting &matarr2 = nmatlist.add(Setting::TypeArray);
		cout << "VT:\n";
		add_mat(matarr2);
		cout << endl;

		Setting &matarr3 = nmatlist.add(Setting::TypeArray);
		cout << "HT:\n";
		add_mat(matarr3);
		cout << endl;
				
		Setting &matarr4 = nmatlist.add(Setting::TypeArray);
		cout << "FT:\n";
		add_mat(matarr4);
		cout << endl;

		Setting &matarr5 = nmatlist.add(Setting::TypeArray);
		cout << "SIG:\n";
		add_mat(matarr5);
		cout << endl;

		Setting &matarr6 = nmatlist.add(Setting::TypeArray);
		cout << "MU:\n";
		add_vec(matarr6);
		cout << endl;

		if (i < (nmats-1)) {
			cout << "Now for the next matrix ... " << endl << endl;
		}
	}
}

void Settings :: set_file() {
	cout << "Enter the name of the file to store the configuration.\n"
"For example 'test1'. The suffix will be added:\n ";
	while (cin.peek() != '\n') {
		cin >> outname;
	}
	outname = outname+".cfg";
}

void Settings :: bsetup(int nmats) {
	cout << "Here we need to modify all the settings. Press Enter to continue:";
	cin.ignore();

	cout << "Enter the dimension of the response variable:\n";
	dim = get_int();
	
	cout << "Enter the total number of LOOP experiments:\n";
	nloops = get_int();
	cout << "Now set the number of iterations of the unique LOOP\n" 
				"experiments and the number of copies of each.\n\n";
	set_vint(loops, nloops);
	cout << endl;

	cout << "Enter the total number of FOREST experiments:\n";
	nforest = get_int();
	cout << "Now set the ensemble size of the unique FOREST experiments\n" 
				"and the number of copies of each.\n\n";
	set_vint(forest, nforest);
	cout << endl;

	if (nmats > 1) { nrest = 1; }
	else {
		cout << "Enter the total number of examples for DATA_LENGTH, NSPLITS,\n"
					"ALPHA and BETA (this is a single integer):\n\n";
		nrest = get_int();
	}
	cout << "Now set the length of the unique predictor vectors, DATA_LENGTH,\n"
			"and the number of copies of each.\n\n";
	set_vint(dl, nrest);
	cout << endl;

	cout << "The total number of experiments for NSPLITS is: " << 
			nrest << ".\nYou can now enter the unique values to be used in\n"
			"the experiments and the number of copies of each.\n\n";
	set_vint(ns, nrest);
	cout << endl;

	cout << "The total number of experiments for ALPHA is: " << 
			nrest << ".\nYou can now enter the unique values to be used in\n"
			"the experiments and the number of copies of each.\n\n";
	set_vfloat(alpha, nrest);
	cout << endl;

	cout << "The total number of experiments for BETA is: " << 
			nrest << ".\nYou can now enter the unique values to be used in\n"
			"the experiments and the number of copies of each.\n\n";
	set_vfloat(beta, nrest);
	cout << endl;
}

void Settings :: getintArray(int marr[], const Setting& mat, 
		int start, int end) {
	for (int i = start; i < end; i++) {
		//cout << "INT " << i << "is " << (int)mat[i] << endl;
		marr[i] = mat[i];
	}
}

void Settings :: getdArray(double marr[], const Setting& mat, 
		int start, int end) {
	for (int i = start; i < end; i++) {
		//cout << "DBL " << i << " is " << (double)mat[i] << endl;
		marr[i] = mat[i];
		//cout << "marr " << i << " is " << marr[i] << endl;
		
	}
}
/*
int main() {
	Settings s;
	Config cfg;
	s.set_file();
	s.make_file(cfg);
//	s.example();
}*/

