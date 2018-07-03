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

void Settings :: set_vint(vector<int>& vec, const int vsize) {
	cout << "There are " << vsize << " to enter.\n";
	vec.clear();
	for (int i = 0; i < vsize; i++) {
		cout << "Please enter an int: \n";
		int input;
		cin >> input;
		if (cin.peek() == '\b') cin.clear();
		vec.push_back(input);
	}
}
 
void Settings :: set_vfloat(vector<float>& vec, const int vsize) {
	cout << "There are " << vsize << " to enter.\n";
	vec.clear();
	for (int i = 0; i < vsize; i++) {
		cout << "Please enter a float: \n";
		float input;
		cin >> input;
		if (cin.peek() == '\b') cin.clear();
		vec.push_back(input);
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
		Setting &matarr = nmatlist.add(Setting::TypeArray);
		double entry;
		cin.ignore();
		cout << "Firstly, W0:\n";
		while (cin.peek() != '\n') {
			cin >> entry;
			matarr.add(Setting::TypeFloat) = entry;
		}
		cout << endl;
		Setting &matarr1 = nmatlist.add(Setting::TypeArray);
		cout << "WT:\n";
		cin.ignore();
		while (cin.peek() != '\n') {
			cin >> entry;
			matarr1.add(Setting::TypeFloat) = entry;
		}
		cout << endl;
		Setting &matarr2 = nmatlist.add(Setting::TypeArray);
		cout << "VT:\n";
		cin.ignore();
		while (cin.peek() != '\n') {
			cin >> entry;
			matarr2.add(Setting::TypeFloat) = entry;
		}
		cout << endl;
		Setting &matarr3 = nmatlist.add(Setting::TypeArray);
		cout << "HT:\n";
		cin.ignore();
		while (cin.peek() != '\n') {
			cin >> entry;
			matarr3.add(Setting::TypeFloat) = entry;
		}
		cout << endl;
		Setting &matarr4 = nmatlist.add(Setting::TypeArray);
		cout << "FT:\n";
		cin.ignore();
		while (cin.peek() != '\n') {
			cin >> entry;
			matarr4.add(Setting::TypeFloat) = entry;
		}
		cout << endl;
		Setting &matarr5 = nmatlist.add(Setting::TypeArray);
		cout << "SIG:\n";
		cin.ignore();
		while (cin.peek() != '\n') {
			cin >> entry;
			matarr5.add(Setting::TypeFloat) = entry;
		}
		cout << endl;
		Setting &matarr6 = nmatlist.add(Setting::TypeArray);
		cout << "and finally MU:\n";
		cin.ignore();
		while (cin.peek() != '\n') {
			cin >> entry;
			matarr6.add(Setting::TypeFloat) = entry;
		}
	}
	cout << endl;
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
	cout << "Enter the number of Loop examples:\n";
	cin >> nloops;
	cout << "Enter the number of Forest examples:\n";
	cin >> nforest;
	if (nmats > 1) {
		nrest = 1;
	}
	else {
		cout << "Enter the number of examples for Data_Length, Nsplits, "
			"Alpha and Beta (this is a single integer):\n";
		cin >> nrest;
	}
	cout << "Finally, enter the dimension of the response variable:\n";
	cin >> dim;

	cout << "Set Loops.\n";
	set_vint(loops, nloops);
	cout << "Set Forest\n";
	set_vint(forest, nforest);
	cout << "Set Data_length\n";
	set_vint(dl, nrest);
	cout << "Set Nsplits\n";
	set_vint(ns, nrest);
	cout << "Set Alpha.\n";
	set_vfloat(alpha, nrest);
	cout << "Set Beta.\n";
	set_vfloat(beta, nrest);
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

