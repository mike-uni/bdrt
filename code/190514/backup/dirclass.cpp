// A CLASS TO PERFORM FILE AND DIRECTORY FUNCTIONS

# pragma once
# include <iostream>
# include <fstream>
# include <vector>
# include <string>
# include <dirent.h>
# include <unistd.h>

using namespace std;

class dirf {
	int wordcount(string s);
	string get_wd();
	inline bool fexists(const char * name);
	bool has_ext(const string& name, const string& ext );
	void get_exts(const string& path, const string& ext);
	void del_files(const string& path);
}

// A function to count the number of words in a string
dirf::int wordcount(string s) {
	int count = 0;
	stringstream ss(s);
	string word;
	while (ss >> word) { count++;};
	return count;
}

// A function to get the current working directory
dirf::string get_wd() {
	string path;
	char nbuf[1024];
	path = getcwd(nbuf,1024);
	return path;
}

// A function to test if a file exists (exact matches only)
dirf::inline bool fexists(const char * name) {
	ifstream infile(name);
	return infile.good();
}

// A function to check if a file has a particular extension
dirf::bool has_ext(const string& name, const string& ext ) {
    return (name.size() >= ext.size()) && equal(ext.rbegin(), 
				ext.rend(), name.rbegin());    
}

// A function to search a directory for a file extension
dirf::void get_exts(const string& path, const string& ext) {
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
dirf::void del_files(const string& path) {
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
