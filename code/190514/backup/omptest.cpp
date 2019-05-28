//OMP TEST

#include<iostream>
#include<map>

using namespace std;

int main() {
	map<int, char*> tmap;
	char myarr[10] = {'a','b','c','d','e','f','g','h','i','j'};
#pragma omp parallel shared(tmap) 
	{	
#pragma omp for
	for (int i = 1; i <= 10; ++i) {
		tmap[i] = new char;
	}
#pragma omp for ordered
	for (int i  = 1; i <= 10; ++i) {
		*tmap[i] = myarr[i - 1];
		char j;
#pragma omp ordered
		j = *tmap[i];
		cout << "first: " << i << " second: " << tmap[i] << endl;
		cout << "char j: " << j << endl;
	}
	}
}
	

