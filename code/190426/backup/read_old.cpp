// A function to read in data

# pragma once
# include <iostream>
# include <cmath>
# include <fstream>
# include <vector>
# include <string>
# include <map>
# include <Eigen/Dense>
# include <random>
# include <dirent.h>
# include <unistd.h>
# include "dist.h"

using namespace std;
using namespace Eigen;

// A function that recieves a string of input 
// and converts it to a double vector type

VectorXd read_line(int start, int end, ifstream& file) {
	int ncols = end - start;
	VectorXd outvec(ncols);
	double val;
	for (int i = start; i < end; i++) {
		file >> val;
		outvec(i%ncols) = val;
	}
	return outvec;
}

// A function that makes varaible names Xi based on number of predictors
vector<string> make_names(int bufsize, int p) {
	vector<string> predname;
	for (int i = 0; i < p; i++) {
    	char buffer[bufsize];
        sprintf(buffer, "X%d", i+1);
		string s = buffer;
		predname.push_back(s);
	}	
	return predname;
}

// A function that generates pairs of values given values and names
map<string, double> make_pairs(vector<string> keys, vector<double> values) {
	map<string, double> myNameVec;
	pair<string, double> myPair;
	for (int i = 0; i < values.size(); i++) {
		myPair = make_pair(keys[i], values[i]);
		myNameVec.insert(myPair);
	}
	return myNameVec;
}

// Get the index of a string in a vector of strings
int get_pred_index(string pred, vector<string> names) {
	int index = 0;
    for (int i = 0; i < names.size(); i++) {
		if (pred == names[i]) {
			index = i;
		}
	}
	return index;
}

// Returns a random key-value pair based on a vector of names and a 
// vector of values
pair<string, double> rand_pair(const MatrixXd& inmat, 
		vector<string> names, int DATA_LENGTH,
		int NSPLITS, int iter) {
	pair<string, double> outpair;
	if (iter == 0) {
		string mpred = names[rand()%DATA_LENGTH];
		int ind;
		ind = get_pred_index(mpred, names);
		double msplit;
		msplit = inmat(0, ind);
		outpair = make_pair(mpred, msplit);
		return outpair;
	}
	else if (iter > 0 && iter < NSPLITS) {
		string mpred = names[rand()%DATA_LENGTH];
		int ind;
		ind = get_pred_index(mpred, names);
		VectorXd lcol = inmat.col(ind);
		double msplit;
		msplit = lcol(rand()%iter); 
		outpair = make_pair(mpred, msplit);
		return outpair;
	}
	else {
		string mpred = names[rand()%DATA_LENGTH];
		int ind;
		ind = get_pred_index(mpred, names);
		VectorXd lcol = inmat.col(ind);
		double msplit;
		msplit = lcol(rand()%NSPLITS); 
		outpair = make_pair(mpred, msplit);
		return outpair;
	}	
}

// Simulate some y values based on predictors, observation variance,
// a matrix that assigns values to componets of output vector and,
// the system linear transformation matrix.
 
VectorXd sim_y(const VectorXd& vec, const MatrixXd& V, 
				const MatrixXd& amat, const MatrixXd& H) {
	VectorXd mu;
	mu = H*(amat*vec);
	
	// get normal random variate
	Dist d;
	VectorXd dev;
	dev =  d.mvn(1, mu, V);
	return dev;
}

// Simulate some predictors given a length of data 
// and write them to a text file

void sim_x(const int P, string filename, int stop) {
	Dist d;
	ofstream outfile;
	outfile.open(filename);
	int rand_int;
	int n = 0;
	while (n < stop) {
		for (int i = 1; i < (P+1); i++) {
			rand_int = d.runif(1, i*3);
			outfile << rand_int << " ";
		}
	outfile << "\n";
	n++;	
	}
}

void sim_normx(const int P, string filename, int stop, double var) {
	Dist d;
	ofstream outfile;
	outfile.open(filename);
	double rand_doub;
	int n = 0;
	while (n < stop) {
		for (int i = 1; i < (P+1); i++) {
			rand_doub = d.rnorm(i, var);
			outfile << rand_doub << " ";
		}
	outfile << "\n";
	n++;	
	}
}

// A function to count the number of words in a string
int wordcount(string s) {
	int count = 0;
	stringstream ss(s);
	string word;
	while (ss >> word) { count++;};
	return count;
}

// A function to get the current working directory
string get_wd() {
	string path;
	char nbuf[1024];
	path = getcwd(nbuf,1024);
	return path;
}

// A function to test if a file exists (exact matches only)
inline bool fexists(const char * name) {
	ifstream infile(name);
	return infile.good();
}

// A function to check if a file has a particular extension
bool has_ext(const string& name, const string& ext ) {
    return (name.size() >= ext.size()) && equal(ext.rbegin(), 
				ext.rend(), name.rbegin());    
}

// A function to search a directory for a file extension
void get_exts(const string& path, const string& ext) {
	DIR * dir = opendir(path.c_str());
	if (!dir) {
		cout << "Directory not found." << endl;
	}
	dirent * entry;
	while ((entry = readdir(dir)) != NULL) {
		if (has_ext(entry->d_name, ext)) {
			cout << entry->d_name << endl;
		}
	}
	cout << "The list is complete or the file does not exist." << endl;
	closedir(dir);
}

// A function to search a directory and delete all the files within.
void del_files(const string& path) {
	DIR * dir = opendir(path.c_str());
	char filepath[1000];
	if (!dir) {
		cout << "Directory not found." << endl;
	}
	dirent * entry;
	while ((entry = readdir(dir)) != NULL) {
		sprintf(filepath, "%s/%s", path.c_str(), entry->d_name);
		remove(filepath); 
	}
	cout << "The folder is empty." << endl;
	closedir(dir);
}

/*
// Simulate a response based on the Friedman function
void sim_f(string filename, int stop, int d, 
			const VectorXd& mu, const MatrixXd& sig) {
	ofstream outfile;
	outfile.open(filename);
	double x


int main() {
	string fileloc = "/home/michael/cpp/181120/output/testin/xsim_150001fixed.txt";
	int DL = 1;
	VectorXd data(DL);
	ifstream xfs(fileloc, ios::in);
	for (int i = 0; i < 10; i++) {
		data = read_line(i*DL, (i+1)*DL, xfs);
		cout << "data: " << data << endl;
	}
}


	//string mpath = get_wd();
	//get_exts(mpath, ".cfg");

	//bool test = fexists(".cpp");
	//cout << test << endl;
	//string fileloc = "/home/michael/Desktop/180416_cpp/output/xsim.txt";
	//sim_x(10, fileloc, 10);
	//ifstream ifs(fileloc, ios::in);	
	//VectorXd dataLine = read_line(0, 10, ifs);
	Vector2d mu_;
	mu_ << 0,0;
	Matrix2d sig_;
	sig_<< 1,0,0,1;//A_*A_.transpose();
	cout << "sig is: " << endl << sig_ << endl;
	int n_ = 1;
	MatrixXd mvn_;
	Dist d;
	int i = 0;
	while (i<10) {
	int randn = d.runif(0, 10);
	cout << randn << endl;
	mvn_ = d.mvn(n_, mu_, sig_);
	cout << "MVN test is: " << endl << mvn_ << endl; 
	i++;
	}
	int randn2 = d.runif(0, 10);
	cout << randn2 << endl;


//	VectorXd dataLine1 = read_line(10, 20, ifs);
//	cout << dataLine1 << endl;
//	VectorXd dataLine2 = read_line(20, 30, ifs);
//	cout << dataLine2 << endl;
//	VectorXd dataLine3 = read_line(30, 40, ifs);
//	cout << dataLine3 << endl;

	//for (int i = 0; i < dataLine.size(); i++) {
	//	cout << dataLine(i)/2 << " ";
	//}
	cout << endl;
	vector<string> myNames = make_names(100, 10);
	for (int i = 0; i < myNames.size(); i++) {
		cout << myNames[i] << " ";
	}
	cout << endl;


	//map<string, double> myNameVec = make_pairs(myNames, dataLine);
	//map<string, double>::iterator it;
	//pair<string, double> myPair;
	//for (int i = 0; i < dataLine.size(); i++) {
	//	myPair = make_pair(myNames[i], dataLine[i]);
	//	myNameVec.insert(myPair);
	//}
	//for (it = myNameVec.begin(); it != myNameVec.end(); ++it) {
	//	cout << "(" << it->first << ", " << it->second << ") ";
	//}
	//cout << endl;
	pair<string, double> nrandp = rand_pair(dataLine, myNames);
	cout << "(" << nrandp.first << ", " << nrandp.second << ")" << endl;
	
	string npred = nrandp.first;
	int myind = get_pred_index(npred, myNames);
	cout << "Index of pred is: " << myind << endl;
	cout << "Value of split is: " << dataLine[myind] << endl;



	MatrixXd A_mat = MatrixXd::Zero(2, 10);
	for (int k=0; k < A_mat.size(); k++) {
		A_mat(k) = rand()%2;
	}
	cout << "A_mat is: " << endl << A_mat << endl;

	VectorXd y_test;
	y_test = sim_y(dataLine, sig_, A_mat);
	cout << "Y_test is:" << endl << y_test << endl;
}*/
