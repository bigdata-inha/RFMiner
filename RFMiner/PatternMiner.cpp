#include "PatternMiner.h"

PatternMiner::PatternMiner() {
	efficient_upperbound = 1;
	node_cnt_ = 0ULL;
}

PatternMiner::~PatternMiner() {}

int PatternMiner::GetNodeCnt() {
	return node_cnt_;
}

void PatternMiner::ClearPatterns() {
	sequential_patterns_.clear();
	debug_ = 0;
}

void PatternMiner::LoadTrainingDatabase(vector<Sequence> database) {
	training_sequence_database_ = database;
	training_sequencde_database_sz_ = static_cast<int>(training_sequence_database_.size());
	double_training_db_sz = static_cast<double>(training_sequencde_database_sz_);
}

void PatternMiner::LoadTestDatabase(vector<Sequence> database) {
	test_sequence_database_ = database;
	test_sequence_database_sz_ = static_cast<int>(test_sequence_database_.size());
}

void PatternMiner::Run(double init_threshold, double threshold, int option) {
	init_threshold_ = init_threshold;
	threshold_ = threshold;
	printf("threshold: %lf\n", threshold_);
	option_ = option;
	node_cnt_++;

	unordered_map<int, unordered_set<int>> inverted_list;
	unordered_map<int, vector<PatternInstance>> item2pattern_instances;

	// Build candidate set and instance db
	for (int id = 0; id < training_sequencde_database_sz_; ++id) {
		vector<int> &S = training_sequence_database_[id].sequence;
		int sz = static_cast<int>(S.size());
		for (int i = 0; i < sz; ++i) {
			const int e = S[i];
			inverted_list[e].insert(id);
			const double single_event_support = static_cast<double>(inverted_list[e].size()) / double_training_db_sz;
			if (single_event_support >= init_threshold_) candidate_events_.insert(e);
			item2pattern_instances[e].push_back(PatternInstance(option_, id, i, i, 0, vector<int>({i})));
		}
	}

	vector<int> pattern;
	double progress = 0;
	printf("start");

	for (const int e : candidate_events_) {
		vector<PatternInstance> pi_set = item2pattern_instances[e];
		pattern.push_back(e);
		const double single_event_support = static_cast<double>(inverted_list[e].size()) / double_training_db_sz;
		RFGrowth(pattern, pi_set, single_event_support);
		pattern.pop_back();
		progress++;
		printf("\r%.2lf complete", 100.0*(progress) / static_cast<double>(candidate_events_.size()));
	}
	if (option_ == 1) cout << "\rNumber of recency patterns: " << static_cast<int>(sequential_patterns_.size()) << "\n";
	else if(option_ == 2) cout << "\rNumber of compactness patterns: " << static_cast<int>(sequential_patterns_.size()) << "\n";
}

vector<Pattern> PatternMiner::GetSequentialPatterns() {
	return sequential_patterns_;
}

//void PatternMiner::WritePatternFile(string filename) {
//	ofstream outfile(filename);
//	if (sequential_patterns_.empty()) return;
//	sort(sequential_patterns_.begin(), sequential_patterns_.end());
//	if (sequential_patterns_[0].type == 1) {
//		int sz = static_cast<int>(sequential_patterns_.size());
//		for (int i = 0; i < sz; ++i) {
//			const auto &entry = sequential_patterns_[i];
//			int pattern_length = static_cast<int>(entry.pattern.size());
//			string line;
//			for (int j = 0; j < pattern_length; ++j) {
//				line += std::to_string(entry.pattern[j]) + " ";
//			}
//			line += "-1 " + std::to_string(entry.interestingness) + " " + std::to_string(entry.frequency) + "\n";
//			for (int j = 0; j < pattern_length - 1; ++j) {
//				line += std::to_string(entry.gap_sequence[j]);
//				if (j + 1 != pattern_length - 1) line += " ";
//			}
//			if (i + 1 != sz) line += "\n";
//			outfile << line;
//		}
//	}
//	else if (sequential_patterns_[0].type == 2) {
//		int sz = static_cast<int>(sequential_patterns_.size());
//		for (int i = 0; i < sz; ++i) {
//			const auto &entry = sequential_patterns_[i];
//			int pattern_length = static_cast<int>(entry.pattern.size());
//			string line;
//			for (int j = 0; j < pattern_length; ++j) {
//				line += std::to_string(entry.pattern[j]) + " ";
//			}
//			line += "-1 " + std::to_string(entry.interestingness) + " " + std::to_string(entry.frequency) + "\n";
//			line += std::to_string(entry.gap_sequence.front()) + "\n";
//			outfile << line;
//		}
//	}
//	else if (sequential_patterns_[0].type == 3) {
//		int sz = static_cast<int>(sequential_patterns_.size());
//		for (int i = 0; i < sz; ++i) {
//			const auto &entry = sequential_patterns_[i];
//			int pattern_length = static_cast<int>(entry.pattern.size());
//			string line;
//			for (int j = 0; j < pattern_length; ++j) {
//				line += std::to_string(entry.pattern[j]) + " ";
//			}
//			line += "-1 " + std::to_string(entry.frequency);
//			if (i + 1 != sz) line += "\n";
//			outfile << line;
//		}
//	}
//	else {
//		printf("Error in type: WriteFile()\n");
//		getchar();
//		exit(-1);
//	}
//}

void PatternMiner::WritePatternFile(string filename) {
	ofstream outfile(filename);
	outfile.precision(20);
	if (sequential_patterns_.empty()) return;
	sort(sequential_patterns_.begin(), sequential_patterns_.end());
	if (sequential_patterns_[0].type == 1) {
		int sz = static_cast<int>(sequential_patterns_.size());
		for (int i = 0; i < sz; ++i) {
			const auto &entry = sequential_patterns_[i];
			int pattern_length = static_cast<int>(entry.pattern.size());
			string line;
			for (int j = 0; j < pattern_length; ++j) {
				outfile << entry.pattern[j] << " ";
				//line += std::to_string(entry.pattern[j]) + " ";
			}
			outfile << "-1 " << std::setprecision(20) << entry.interestingness << " " << std::setprecision(20) << entry.frequency << " " << std::setprecision(20)<<entry.confidence << "\n";
			//line += "-1 " + std::to_string(entry.interestingness) + " " + std::to_string(entry.frequency) + " " + std::to_string(entry.confidence) + "\n";
			for (int j = 0; j < pattern_length - 1; ++j) {
				outfile << entry.gap_sequence[j];
				//line += std::to_string(entry.gap_sequence[j]);
				if (j + 1 != pattern_length - 1) {
					//line += " ";
					outfile << " ";
				}
			}
			if (i + 1 != sz) {
				outfile << "\n";
				//line += "\n";
			}
			//outfile << line;
		}
	}
	else if (sequential_patterns_[0].type == 2) {
		int sz = static_cast<int>(sequential_patterns_.size());
		for (int i = 0; i < sz; ++i) {
			const auto &entry = sequential_patterns_[i];
			int pattern_length = static_cast<int>(entry.pattern.size());
			for (int j = 0; j < pattern_length; ++j) {
				outfile << entry.pattern[j] << " ";
			}
			outfile << "-1 " << std::setprecision(20) << entry.interestingness << " " << std::setprecision(20) << entry.frequency << " " << std::setprecision(20) << entry.confidence << "\n";
			outfile << entry.gap_sequence.front() << "\n";
		}
	}
}

void PatternMiner::WriteQueryFile(string filename) {
	ofstream outfile(filename);
	if (test_sequence_database_.empty()) return;
	for (int i = 0; i < test_sequence_database_sz_; ++i) {
		const auto& entry = test_sequence_database_[i];
		outfile << entry.sequence_id << " ";

		for (const auto& item : entry.sequence) {
			outfile << item << " ";
		}
		outfile << "-1";
		if (i + 1 != test_sequence_database_sz_) outfile << "\n";
	}
}


void PatternMiner::RFGrowth(vector<int> pattern, vector<PatternInstance> pi_set, double pre_interesting) {
	node_cnt_++;

	const int pattern_length = static_cast<int>(pattern.size());
	const int plus_pattern_len = pattern_length + 1;
	const int gap_length = pattern_length;

	if (debug_) {
		printf("Before Pattern: <");
		for (int i = 0; i < pattern.size(); ++i) printf("%d ", pattern[i]);
		printf(">\n");
		/*for (const auto &entry : pi_set) {
			printf("[%d %d]\n", entry.l, entry.r);
		}*/
	}

	unordered_map<int, unordered_map<int, int>> max_length_map;
	unordered_map<int, double> min_instance;
	unordered_map<int, unordered_set<int>> inverted_list;

	unordered_map<int, vector<PatternInstance>> item2pattern_instances = GrowNew(pattern, pi_set, pattern.back(), max_length_map, min_instance, inverted_list);

	for (const auto &entry : item2pattern_instances) {
		const int e = entry.first;
		pattern.push_back(e); // P+ <- appending event e to pattern P
		vector<PatternInstance> cur_pi_set = entry.second; // Grow()
		
		if (debug_) {
			printf("Current Pattern: <");
			for (int i = 0; i < pattern.size(); ++i) printf("%d ", pattern[i]);
			printf(">\n");
		}

		vector<double> pattern_gaps;
		pattern_gaps.resize(pattern_length);
		double interestingness = 0.0;
		double len = 0.0;

		if (option_ == 1) {
			unordered_map<int, pair<vector<double>, double>> compact_gaps;
			unordered_map<int, pair<double, double>> compact_weights;
			
			for (const auto& instance : cur_pi_set) {
				const int id = instance.sid;
				vector<double> gaps = instance.GapSequence();			
				double ppl = static_cast<double>(plus_pattern_len);
				compact_weights[id].first += instance.value;
				compact_weights[id].second++;
				if (compact_gaps[id].second == 0) {
					compact_gaps[id].first.resize(gap_length);
				}
				int sz = static_cast<int>(gaps.size());
				for (int i = 0; i < sz; ++i) compact_gaps[id].first[i] += gaps[i];
				compact_gaps[id].second++;
			}

			// Gap Sequence Calculation
			for (auto &entry : compact_gaps) {
				int sz = static_cast<int>(entry.second.first.size());
				for (int i = 0; i < sz; ++i) {
					pattern_gaps[i] += (entry.second.first[i] / entry.second.second);
				}
			}
			for (auto &entry : pattern_gaps) entry /= static_cast<double>(compact_gaps.size());

			// Interestingness
			for (auto &entry : compact_weights) interestingness += entry.second.first / static_cast<double>(entry.second.second);
			interestingness /= static_cast<double>(training_sequencde_database_sz_);
		}
		else if(option_ == 2){
			unordered_map<int, pair<double, double>> tracker;
			unordered_map<int, pair<double, double>> length_tracker;

			for (const auto & entry : cur_pi_set) {
				int id = entry.sid;
				tracker[id].first += entry.value;
				tracker[id].second++;
				length_tracker[id].first += entry.denom;
				length_tracker[id].second++;
			}

			// Compactness [l: r]. average(r - l + 1)
			for (const auto & entry : length_tracker) len += (entry.second.first / entry.second.second);
			len /= static_cast<double>(length_tracker.size());
			
			// Interestingness
			for (const auto & entry : tracker) interestingness += (entry.second.first / entry.second.second);
			interestingness /= static_cast<double>(training_sequencde_database_sz_);
		}
	
		unordered_set<int> vis;
		for (const auto & entry : cur_pi_set) vis.insert(entry.sid);
		
		if(option_ == 3) interestingness = static_cast<double>(vis.size()) / static_cast<double>(training_sequencde_database_sz_);
		double frequency_support = static_cast<double>(vis.size()) / static_cast<double>(training_sequencde_database_sz_);

		double confidence = interestingness / pre_interesting;
		double max1_confidence = confidence < 1.0 ? confidence : 1.0;

		if (interestingness >= threshold_) {
			if (option_ == 1) sequential_patterns_.push_back(Pattern(option_, pattern, interestingness, frequency_support, pattern_gaps, confidence));
			else if (option_ == 2) {
				vector<double> avglen;
				avglen.push_back(len);
				sequential_patterns_.push_back(Pattern(option_, pattern, interestingness, frequency_support, avglen, confidence));
			}
			RFGrowth(pattern, cur_pi_set, interestingness);
		}
		else if(option_ != 3){
			double upper_bound;
			if (!efficient_upperbound) upper_bound = frequency_support;
			else upper_bound = UpperBound(static_cast<double>(plus_pattern_len), max_length_map[e], min_instance[e]);
			if(upper_bound >= threshold_)RFGrowth(pattern, cur_pi_set, interestingness);
		}
		pattern.pop_back();
	}
}

unordered_map<int, vector<PatternInstance>> PatternMiner::GrowNew(vector<int> pattern, vector<PatternInstance> pi_set, int last_event, unordered_map<int, unordered_map<int, int >> &max_length_map, unordered_map<int, double> &min_instance, unordered_map<int, unordered_set<int>> &inverted_list) {
	unordered_map<int, vector<PatternInstance>> ret;

	int pi_set_sz = static_cast<int>(pi_set.size());

	for (int i = 0; i < pi_set_sz; ++i) {
		const PatternInstance &pi = pi_set[i];
		const int id = pi.sid;
		vector<int> landmark = pi.landmark;

		const vector<int> &S = training_sequence_database_[id].sequence;
		int sz = static_cast<int>(S.size());

		unordered_set<int> leftmost_vis;

		// check if ajacent pattern instances have the same id
		if (i + 1 < pi_set_sz && id == pi_set[i + 1].sid) {
			sz = std::min(sz, pi_set[i + 1].r + 1);
		}

		bool first = false;
		for (int j = pi.r + 1; j < sz; ++j) {
			const int e = S[j];

			// if not candidate events, then discard
			if (candidate_events_.find(e) == candidate_events_.end()) continue;

			// [a x x x b x x c x x c]
			if (leftmost_vis.find(e) == leftmost_vis.end()) {
				leftmost_vis.insert(e);
				int ext_len = static_cast<int>(S.size()) - j;
				vector<int> landmark;
				landmark.push_back(j);
				int pos = pattern.size() - 1;
				for (int k = j - 1; k > pi.l; --k) {
					if (pos == 0) break;
					if (S[k] == pattern[pos]) {
						--pos;
						landmark.push_back(k);
					}
				}
				landmark.push_back(pi.l);
				std::reverse(landmark.begin(), landmark.end());
				ret[e].push_back(PatternInstance(option_, id, pi.l, j, ext_len, landmark));
				if(ext_len != 0) max_length_map[e][id] = std::max(max_length_map[e][id], ext_len);
				if (min_instance.find(e) == min_instance.end()) {
					min_instance[e] = 100000000.0;
				}
				min_instance[e] = std::min(min_instance[e], ret[e].back().denom);
			}
		}
	}
	return ret;
}

double PatternMiner::UpperBound(double plus_pattern_len, unordered_map<int, int > max_length_map, double min_instance) {
	double max_upper = 0.0;
	map<int, double> upper_bound_map;

	// calculate N_k
	for (auto it = max_length_map.begin(); it != max_length_map.end(); ++it) {
		const int ext_length = it->second;
		upper_bound_map[-ext_length]++;
	}

	for (auto it = upper_bound_map.begin(); it != upper_bound_map.end(); ++it) {
		int cnt = it->second;
		++it;
		if (it == upper_bound_map.end()) break;
		it->second += cnt;
		--it;
	}
	for (auto it = upper_bound_map.begin(); it != upper_bound_map.end(); ++it) {
		auto it2 = it;
		it++;
		if (it == upper_bound_map.end()) break;
		assert(it->second >= it2->second);
	}

	if (option_ == RECENCY) {
		double P = (double)(plus_pattern_len);
		for (const auto &entry : upper_bound_map) {
			double K = (double)(-entry.first);
			double N = entry.second;
			double V = min_instance;
			double a = (1.0 + P + K)*(P + K) / 2.0;
			double b = K*P;
			double c = (1.0 + K)*K / 2.0;
			double single_best_val = a / (V + b + c);
			max_upper = std::max(max_upper, single_best_val*N);
		}
		max_upper /= static_cast<double>(training_sequencde_database_sz_);
	}
	else if (option_ == COMPACTNESS) {
		double P = (double)(plus_pattern_len);
		for (const auto &entry : upper_bound_map) {
			double K = (double)(-entry.first);
			double N = entry.second;
			double V = min_instance;
			double a = P + K;
			double b = V + K;
			double single_best_val = a / b;
			max_upper = std::max(max_upper, single_best_val*N);
		}
		max_upper /= static_cast<double>(training_sequencde_database_sz_);
	}
	return max_upper;
}

void PatternMiner::SetDebug() {
	debug_ = 1;
}