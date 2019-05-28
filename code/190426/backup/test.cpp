// TEST ALGORITHMS

# pragma once
# include <iostream>
# include <cmath>
# include <algorithm>
# include <vector>
# include <set>
# include <map>
//# include "move.h"
//# include "tree.h"
//# include "node.h"

using namespace std;

int main() {
	map<int, int> tree;
	for (int i = 0; i < 5; i++) {
		tree[i] = i*10;
	}

	for(pair<int, int> item : tree) {
		cout << item.second << " ";
	}
	cout << endl;

	double num = 1.3456;
	int num1 = 3;
	cout << "num " << num%3 << "num1 " << num1%3 << endl;
}
/*
	vector<int> leaf_nodes = {2,3};
	int const num_nodes = leaf_nodes.size();

	vector<int>::iterator it;

	// Get leaves that are a pair.
	set<int> dual_leaf_set;
	for (int i = 0; i < num_nodes; ++i) {
		int j = leaf_nodes[i];
		it = find(leaf_nodes.begin(), leaf_nodes.end(), j+1);
		if (it != leaf_nodes.end()) {
			dual_leaf_set.insert(j);
			dual_leaf_set.insert(j+1);
		}
	}

	vector<int> dual_leaves;
	dual_leaves.assign(dual_leaf_set.begin(), dual_leaf_set.end());
	int num_dual = dual_leaves.size();

	int const rand_node = dual_leaves[rand()%num_dual];
	for (int i = 0; i < num_dual; ++i) {
		cout << dual_leaves[i] << endl;
	}
	
	map<int, Tree*> forest;
	for (int i = 0; i < 4; ++i) {
		forest[i]= new Tree(1, 1);
	}

	for (int i = 0; i < 4; ++i) {
		cout << forest[i]->tree[1]->nnode << endl;
	}
}*/
