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

#include <boost/functional/hash.hpp>
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
	int offset;
	int transaction_id;
	PseudoProjectedDatabaseEntry(int a_offset, int a_transaction_id)
		:offset(a_offset), transaction_id(a_transaction_id) {
	}
};

struct DatabaseSequenceIterator {
	int transaction_id;
	vector<int> offsets;
	double current_score;
	bool offset_change;

	DatabaseSequenceIterator()
		:transaction_id(-1), current_score(0.0), offset_change(false) {}

	void PushBack(int offset) {
		offsets.push_back(offset);
		offset_change = true;
	}

	double ScoreCount(double len){
		if (offset_change) {
			current_score = 0.0;

			int sz = static_cast<int>(offsets.size());
			for (int i = 0; i < sz - 1; ++i) {
				double tmp = static_cast<double>(offsets[i + 1] - offsets[i]);
				current_score += 1.0 / tmp;
			}
			offset_change = false;
		}
		return (current_score + 1.0 / len)/static_cast<double>(offsets.size());
	}
};

class PrefixSpan{

public:

	// constructors, destructors, initilizeres -------------------
	PrefixSpan() {}
	~PrefixSpan() {}

	void clear_training() {
		frequent_patterns_.clear();
	}

	void load_train_database(vector<Sequence> database) {
		training_transaction_database_ = database;
		traing_transaction_database_size_ = static_cast<int>(training_transaction_database_.size());
	}

	void LoadTestDatabase(vector<Sequence> database) {
		test_database_ = database;
		test_database_size_ = static_cast<int>(test_database_.size());
	}

	// -- Support version
	void RunFrequencyVersion(double frequency_minsup);
	void SpanFrequency(const int &pre_item, const double &pre_support,  list<PseudoProjectedDatabaseEntry> &pre_pseudoprojected_database, const vector<Sequence>& transaction_database, vector<int>&pattern, const int &depth);

	// getters, setters --------------------------------------------
	vector<Pattern> GetFrequentPatterns() {
		return frequent_patterns_;
	}

	// -------------------------------------------------------------

	// IO functions ------------------------------------------------
	void WriteFile(const string &filename);
	void WriteQueryFile(const string &filename);
	// -------------------------------------------------------------

private:
	// Frequency version
	double frequency_minsup_;

	// transaction database
	vector<Sequence> training_transaction_database_;
	int traing_transaction_database_size_;

	int test_database_size_;
	vector<Sequence> test_database_;

	// frequent patterns
	vector<Pattern> frequent_patterns_;

	// class support functions
	unordered_map<int, unordered_set<int>> ScanCountSingleItemsInit();

	vector<Sequence> CreateOptimizedTransactionDB(unordered_map<int, unordered_set<int>> &item2seq_hash_table);
	list<PseudoProjectedDatabaseEntry> BuilProjectedDatabaseFrequency(const int &item, const list<PseudoProjectedDatabaseEntry> &cur_pseudoprojected_database, const unordered_set<int>& sidset, const vector<Sequence>& transaction_database);

	unordered_map<int, unordered_set<int>> FrequencyUtil(const list<PseudoProjectedDatabaseEntry> &cur_pseudoprojected_database, const vector<Sequence> &transaction_database);

	
	double PotentialSupport(const int &pattern_len, const int &num_frequency, const int &max_len);
};