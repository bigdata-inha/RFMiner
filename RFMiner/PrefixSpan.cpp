#include "Database.h"
#include "PrefixSpan.h"

// support version -------------------------------------------------------------------------------------------------------------------------------
void PrefixSpan::RunFrequencyVersion(double frequency_minsup) {
	frequency_minsup_ = frequency_minsup;
	cout << "running on minsup " << frequency_minsup_ << "\n";

	auto inverted_list = ScanCountSingleItemsInit();
	auto optimized_transaction_db = training_transaction_database_;

	// starting projected_database
	list<PseudoProjectedDatabaseEntry> init_pseudo_projected_database;

	int sz = static_cast<int>(optimized_transaction_db.size());
	for (int i = 0; i < sz; ++i) {
		init_pseudo_projected_database.push_back(PseudoProjectedDatabaseEntry(0, i));
	}

	printf("start");
	vector<int> pattern;

	int progress = 0;

	for (const auto& entry : inverted_list) {
		++progress;
		const int &item = entry.first;
		const auto &transaction_id_list = entry.second;
		const double &support = static_cast<double>(transaction_id_list.size()) / static_cast<double>(traing_transaction_database_size_);

		if (support >= frequency_minsup_) {
			pattern.push_back(item);
			//frequent_patterns_.push_back(Pattern(vector<int>(pattern), support));
			auto pseudoprojected_database = BuilProjectedDatabaseFrequency(item, init_pseudo_projected_database, transaction_id_list, optimized_transaction_db);
			SpanFrequency(item, support, pseudoprojected_database, optimized_transaction_db, pattern, 1);
			pattern.pop_back();
			printf("\r%.2lf complete", 100.0*((double)(progress)) / static_cast<double>(inverted_list.size()));
		}
	}

	cout << "\rNumber of frequent patterns: " << static_cast<int>(frequent_patterns_.size()) << "\n";
}

void PrefixSpan::SpanFrequency(const int &pre_item, const double &pre_support, list<PseudoProjectedDatabaseEntry> &pre_pseudoprojected_database, const vector<Sequence>& transaction_database, vector<int>&pattern, const int &depth) {


	auto inverted_list = FrequencyUtil(pre_pseudoprojected_database, transaction_database);

	cnt++;
	// for all counted items
	for (const auto&entry : inverted_list) {
		const int &item = entry.first;
		const auto &transaction_id_list = entry.second;
		const double &support = static_cast<double>(transaction_id_list.size()) / static_cast<double>(traing_transaction_database_size_);
	
		if (support >= frequency_minsup_) {
			pattern.push_back(item);
			frequent_patterns_.push_back(Pattern(3, vector<int>(pattern), support, support));
			auto pseudoprojected_database = BuilProjectedDatabaseFrequency(item, pre_pseudoprojected_database, transaction_id_list, transaction_database);
			SpanFrequency(item, support, pseudoprojected_database, transaction_database, pattern, depth + 1);
			pattern.pop_back();
		}
	}
}


unordered_map<int, unordered_set<int>> PrefixSpan::ScanCountSingleItemsInit() {
	unordered_map<int, unordered_set<int>> ret;
	
	for (int i = 0; i < traing_transaction_database_size_; ++i) {
		for (const auto& item : training_transaction_database_[i].sequence) {
			ret[item].insert(i);
		}
	}
	return ret;
}

unordered_map<int, unordered_set<int>> PrefixSpan::FrequencyUtil(const list<PseudoProjectedDatabaseEntry> &cur_pseudoprojected_database, const vector<Sequence>& transaction_database) {
	unordered_map<int, unordered_set<int>> ret;
	for (auto it = cur_pseudoprojected_database.begin(); it != cur_pseudoprojected_database.end(); ++it) {
		int offset = it->offset;
		int sid = it->transaction_id;
		for (int i = offset; i < transaction_database[sid].sequence.size(); ++i) {
			ret[transaction_database[sid].sequence[i]].insert(sid); // this stops counting repeated items cuz set DT.
		}
	}
	return ret;
}

vector<Sequence> PrefixSpan::CreateOptimizedTransactionDB(unordered_map<int, unordered_set<int>> &inverted_list) {
	vector<Sequence> ret;
	for (const auto& entry : training_transaction_database_) {
		int sequence_id = entry.sequence_id;
		vector<int> optimized_sequence;
		for (const auto& item : entry.sequence) {
			if (static_cast<int>(inverted_list[item].size()) >= frequency_minsup_) {
				optimized_sequence.push_back(item);
			}
		}
		if (optimized_sequence.empty()) continue;
		ret.push_back(Sequence(sequence_id, optimized_sequence));
	}
	return ret;
}


list<PseudoProjectedDatabaseEntry> PrefixSpan::BuilProjectedDatabaseFrequency(const int &item, const list<PseudoProjectedDatabaseEntry> &cur_pseudoprojected_database, const unordered_set<int>& sidset, const vector<Sequence>& transaction_database) {
	list<PseudoProjectedDatabaseEntry> ret;

	for (auto it = cur_pseudoprojected_database.begin(); it != cur_pseudoprojected_database.end(); ++it) {
		int offset = it->offset;
		int transaction_id = it->transaction_id;

		if (sidset.find(transaction_id) == sidset.end()) continue;
		
		auto &transaction = transaction_database[transaction_id].sequence;
		int transaction_length = static_cast<int>(transaction.size());
		for (int i = offset; i < transaction_length; ++i) {
			if (transaction[i] == item) {
				ret.push_back(PseudoProjectedDatabaseEntry(i + 1, transaction_id));
				break;
			}
		}
	}
	return ret;
}

void PrefixSpan::WriteFile(const string &filename) {
	ofstream outfile(filename);
	if (frequent_patterns_.empty()) return;
	sort(frequent_patterns_.begin(), frequent_patterns_.end());
	if (frequent_patterns_[0].type == 1) {
		int sz = static_cast<int>(frequent_patterns_.size());
		for (int i = 0; i < sz; ++i) {
			const auto &entry = frequent_patterns_[i];
			int pattern_length = static_cast<int>(entry.pattern.size());
			string line;
			for (int j = 0; j < pattern_length; ++j) {
				line += std::to_string(entry.pattern[j]) + " ";
			}
			line += "-1 " + std::to_string(entry.interestingness) + " " + std::to_string(entry.frequency) + "\n";
			for (int j = 0; j < pattern_length - 1; ++j) {
				line += std::to_string(entry.gap_sequence[j]);
				if (j + 1 != pattern_length - 1) line += " ";
			}
			if (i + 1 != sz) line += "\n";
			outfile << line;
		}
	}
	else if (frequent_patterns_[0].type == 2) {
		int sz = static_cast<int>(frequent_patterns_.size());
		for (int i = 0; i < sz; ++i) {
			const auto &entry = frequent_patterns_[i];
			int pattern_length = static_cast<int>(entry.pattern.size());
			string line;
			for (int j = 0; j < pattern_length; ++j) {
				line += std::to_string(entry.pattern[j]) + " ";
			}
			line += "-1 " + std::to_string(entry.interestingness) + " " + std::to_string(entry.frequency);
			if (i + 1 != sz) line += "\n";
			outfile << line;
		}
	}
	else if (frequent_patterns_[0].type == 3) {
		int sz = static_cast<int>(frequent_patterns_.size());
		for (int i = 0; i < sz; ++i) {
			const auto &entry = frequent_patterns_[i];
			int pattern_length = static_cast<int>(entry.pattern.size());
			string line;
			for (int j = 0; j < pattern_length; ++j) {
				line += std::to_string(entry.pattern[j]) + " ";
			}
			line += "-1 " + std::to_string(entry.frequency);
			if (i + 1 != sz) line += "\n";
			outfile << line;
		}
	}
	else {
		printf("Error in type: WriteFile()\n");
		getchar();
		exit(-1);
	}
}

void PrefixSpan::WriteQueryFile(const string &filename) {
	unordered_map<int, unordered_set<int>> inverted_list;
	for (int i = 0; i < test_database_size_; ++i) {
		for (const auto& item : test_database_[i].sequence) {
			inverted_list[item].insert(i);
		}
	}
	vector<Sequence> opt_test_database;
	for (const auto& entry : test_database_) {
		vector<int> optimized_sequence;
		int sequence_id = entry.sequence_id;
		for (const auto& item : entry.sequence) {
			if (static_cast<int>(inverted_list[item].size()) >= frequency_minsup_) {
				optimized_sequence.push_back(item);
			}
		}
		if (optimized_sequence.empty()) continue;
		opt_test_database.push_back(Sequence(sequence_id, optimized_sequence));
	}

	ofstream outfile(filename);
	if (opt_test_database.empty()) return;
	int sz = static_cast<int>(opt_test_database.size());
	for (int i = 0; i < sz; ++i) {
		const auto& entry = opt_test_database[i];
		outfile << entry.sequence_id << " ";

		for (const auto& item : entry.sequence) {
			outfile << item << " ";
		}
		outfile << "-1";
		if (i + 1 != sz) outfile << "\n";
	}
}