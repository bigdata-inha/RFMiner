#include "Database.h"
#include "PrefixSpan.h"

PrefixSpan::PrefixSpan() { 
	node_cnt = 1ULL;
	naive = 0;
}
PrefixSpan::~PrefixSpan(){}

void PrefixSpan::ClearTraining() {
	frequent_patterns_.clear();
}

int PrefixSpan::GetNodeCnt() {
	return node_cnt;
}

void PrefixSpan::LoadTrainDatabase(vector<Sequence> database) {
	training_database_ = database;
	training_database_size_ = static_cast<int>(training_database_.size());
	double_training_db_sz = static_cast<double>(training_database_size_);
}

void PrefixSpan::LoadTestDatabase(vector<Sequence> database) {
	test_database_ = database;
	test_database_size_ = static_cast<int>(test_database_.size());
}

vector<Pattern> PrefixSpan::GetFrequentPatterns() {
	return frequent_patterns_;
}

vector<unordered_set<int>> PrefixSpan::GetPatternCovers() {
	return pattern_covers_;
}


void PrefixSpan::Run(double frequency_minsup) {
	frequency_minsup_ = frequency_minsup;
	node_cnt++;
	cout << "running on minsup " << frequency_minsup_ << "\n";

	// Counting
	unordered_map<int, unordered_set<int>> inverted_list;
	for (int id = 0; id < training_database_size_; ++id) {
		for (const auto& item : training_database_[id].sequence) {
			inverted_list[item].insert(id);
		}
	}

	// starting projected_database
	vector<PseudoProjectedDatabaseEntry> init_pseudo_projected_database;
	for (int id = 0; id < training_database_size_; ++id) {
		init_pseudo_projected_database.push_back(PseudoProjectedDatabaseEntry(0, id));
	}

	printf("start");
	vector<int> pattern;

	int progress = 0;

	unordered_set<int> candidate_set;
	for (const auto& entry : inverted_list) {
		const int &item = entry.first;
		const auto &transaction_id_list = entry.second;
		const double &support = static_cast<double>(transaction_id_list.size()) / double_training_db_sz;
		if (support >= frequency_minsup_) {
			candidate_set.insert(item);
		}
	}

	for (const auto& entry : candidate_set) {
		++progress;
		const int &item = entry;
		const auto &transaction_id_list = inverted_list[item];
		const double &support = static_cast<double>(transaction_id_list.size()) / double_training_db_sz;

		if (support >= frequency_minsup_) {
			pattern.push_back(item);
			//frequent_patterns_.push_back(Pattern(3, vector<int>(pattern), support, support));
			auto pseudoprojected_database = BuilProjectedDatabase(item, init_pseudo_projected_database);
			Span(pattern, pseudoprojected_database, support);
			pattern.pop_back();
			printf("\r%.2lf complete", 100.0*((double)(progress)) / static_cast<double>(inverted_list.size()));
		}
	}

	cout << "\rNumber of frequent patterns: " << static_cast<int>(frequent_patterns_.size()) << "\n";
}

void PrefixSpan::Span(vector<int>pattern, vector<PseudoProjectedDatabaseEntry> &pre_pseudoprojected_database, double pre_support) {

	// Counting from projected database
	unordered_map<int, unordered_set<int>> inverted_list;
	unordered_set<int> candidate_set;

	for (auto it = pre_pseudoprojected_database.begin(); it != pre_pseudoprojected_database.end(); ++it) {
		int offset = it->l;
		int sid = it->id;
		vector<int> &S = training_database_[sid].sequence;
		int sz = static_cast<int>(S.size());
		for (int i = offset; i < sz; ++i) {
			inverted_list[S[i]].insert(sid); // this stops counting repeated items cuz set DT.
			double a = static_cast<double>(inverted_list[S[i]].size());
			double support = a / double_training_db_sz;
			if (support >= frequency_minsup_) candidate_set.insert(S[i]);
		}
	}

	node_cnt++;
	// for all counted items
	for (const auto&entry : candidate_set) {
		const int &item = entry;
		const auto &id_list = inverted_list[item];
		const double support = static_cast<double>(inverted_list[entry].size())/double_training_db_sz;
		const double confidence = support / pre_support;

		pattern.push_back(item);
		frequent_patterns_.push_back(Pattern(3, vector<int>(pattern), support, support, confidence));
		if(naive) pattern_covers_.push_back(inverted_list[entry]);
		auto pseudoprojected_database = BuilProjectedDatabase(item, pre_pseudoprojected_database);
		Span(pattern, pseudoprojected_database, support);
		pattern.pop_back();
		
	}
}


vector<PseudoProjectedDatabaseEntry> PrefixSpan::BuilProjectedDatabase(const int item, const vector<PseudoProjectedDatabaseEntry> &cur_pseudoprojected_database) {
	vector<PseudoProjectedDatabaseEntry> ret;

	for (auto it = cur_pseudoprojected_database.begin(); it != cur_pseudoprojected_database.end(); ++it) {
		int offset = it->l;
		int transaction_id = it->id;
		
		auto &transaction = training_database_[transaction_id].sequence;
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

void PrefixSpan::WriteFile(const string filename) {
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
			line += "-1 " + std::to_string(entry.frequency) + " " + std::to_string(entry.confidence);
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

void PrefixSpan::WriteFile2(const string filename) {
	ofstream outfile(filename);
	if (frequent_patterns_.empty()) return;
	sort(frequent_patterns_.begin(), frequent_patterns_.end());
	int sz = static_cast<int>(frequent_patterns_.size());
	for (int i = 0; i < sz; ++i) {
		const auto &entry = frequent_patterns_[i];
		int pattern_length = static_cast<int>(entry.pattern.size());
		string line;
		for (int j = 0; j < pattern_length; ++j) {
			line += std::to_string(entry.pattern[j]) + " ";
		}
		line += "-1 " + std::to_string(entry.frequency) + " " + std::to_string(entry.confidence);
		if (i + 1 != sz) line += "\n";
		outfile << line;
	}
	
}



void PrefixSpan::WriteQueryFile(const string filename) {
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