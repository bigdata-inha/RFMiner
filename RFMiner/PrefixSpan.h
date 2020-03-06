#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#include <algorithm>
#include <bitset>
#include <chrono>

#include "Database.h"
#include "StructTypes.h"

using std::vector;
using std::string;
using std::ifstream;
using std::list;
using std::pair;
using std::unordered_map;
using std::ostream;
using std::ofstream;
using std::unordered_set;
using std::set;
using std::bitset;
using std::priority_queue;
using std::cout;

struct PseudoProjectedDatabaseEntry {
	int id, l;
	PseudoProjectedDatabaseEntry(int a_l, int a_id)
		:l(a_l), id(a_id) {
	}
};

class PrefixSpan{

public:
	PrefixSpan();
	~PrefixSpan();

	void ClearTraining();

	void LoadTrainDatabase(vector<Sequence> database);
	void LoadTestDatabase(vector<Sequence> database);
	vector<Pattern> GetFrequentPatterns();
	vector<unordered_set<int>> GetPatternCovers();

	void Run(double frequency_minsup);
	void Span(vector<int> pattern, vector<PseudoProjectedDatabaseEntry> &pre_pseudoprojected_database);

	void WriteFile(const string filename);
	void WriteQueryFile(const string filename);
	
	int cnt;
	int naive;
private:
	// Frequency version
	double frequency_minsup_;

	// sequence database
	vector<Sequence> training_database_;
	int training_database_size_;
	double double_training_db_sz;

	vector<Sequence> test_database_;
	int test_database_size_;

	// frequent patterns
	vector<Pattern> frequent_patterns_;
	vector<unordered_set<int>> pattern_covers_;

	// class support functions
	vector<PseudoProjectedDatabaseEntry> BuilProjectedDatabase(const int item, const vector<PseudoProjectedDatabaseEntry> &cur_pseudoprojected_database);
};