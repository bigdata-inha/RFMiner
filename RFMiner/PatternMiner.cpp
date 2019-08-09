#include "PatternMiner.h"

PatternMiner::PatternMiner() {}

PatternMiner::~PatternMiner() {}

void PatternMiner::ClearPatterns() {
	sequential_patterns_.clear();
	debug_ = 0;
}

void PatternMiner::LoadTrainingDatabase(vector<Sequence> database) {
	training_sequence_database_ = database;
	training_sequencde_database_sz_ = static_cast<int>(training_sequence_database_.size());
}

void PatternMiner::LoadTestDatabase(vector<Sequence> database) {
	test_sequence_database_ = database;
	test_sequence_database_sz_ = static_cast<int>(test_sequence_database_.size());
}

void PatternMiner::Run(double init_threshold, double threshold, int option) {
	init_threshold_ = init_threshold;
	threshold_ = threshold;
	option_ = option;

	auto inverted_list = ScanCountSingleItemsInit();
	vector<int> pattern;
	int progress = 0;
	printf("start");
	for (const auto &entry : inverted_list) {

		const double support = static_cast<double>(entry.second.size()) / static_cast<double>(training_sequencde_database_sz_);
		
		if (support >= init_threshold_) {
			pattern.push_back(entry.first);

			// construct initial pattern instance set
			unordered_set<int> sid_set;
			unordered_set<PI::PatternInstance, boost::hash<PI::PatternInstance>> pi_set;
			for (const auto &sid : entry.second) {
				/*
					Frequency measure does not need multiple pattern instances in a single sequence.
				*/
				if (option_ == 3) {
					if (sid_set.find(sid) != sid_set.end()) continue;
					sid_set.insert(sid);
				}
				
				vector<int> &S = training_sequence_database_[sid].sequence;
				int sz = static_cast<int>(S.size());
				for (int i = 0; i < sz; ++i) {
					if (S[i] == entry.first) pi_set.insert(PI::PatternInstance(sid, i, i));
				}
			}

			RFGrowth(pattern, pi_set);
			pattern.pop_back();
		}
		progress++;
		printf("\r%.2lf complete", 100.0*((double)(progress)) / static_cast<double>(inverted_list.size()));
	}
}

vector<Pattern> PatternMiner::GetSequentialPatterns() {
	return sequential_patterns_;
}

void PatternMiner::WritePatternFile(string filename) {
	ofstream outfile(filename);
	if (sequential_patterns_.empty()) return;
	sort(sequential_patterns_.begin(), sequential_patterns_.end());
	if (sequential_patterns_[0].type == 1) {
		int sz = static_cast<int>(sequential_patterns_.size());
		for (int i = 0; i < sz; ++i) {
			const auto &entry = sequential_patterns_[i];
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
	else if (sequential_patterns_[0].type == 2) {
		int sz = static_cast<int>(sequential_patterns_.size());
		for (int i = 0; i < sz; ++i) {
			const auto &entry = sequential_patterns_[i];
			int pattern_length = static_cast<int>(entry.pattern.size());
			string line;
			for (int j = 0; j < pattern_length; ++j) {
				line += std::to_string(entry.pattern[j]) + " ";
			}
			line += "-1 " + std::to_string(entry.interestingness) + " " + std::to_string(entry.frequency) + "\n";
			line += std::to_string(entry.gap_sequence.front()) + "\n";
			outfile << line;
		}
	}
	else if (sequential_patterns_[0].type == 3) {
		int sz = static_cast<int>(sequential_patterns_.size());
		for (int i = 0; i < sz; ++i) {
			const auto &entry = sequential_patterns_[i];
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

unordered_map<int, unordered_set<int>> PatternMiner::ScanCountSingleItemsInit() {
	unordered_map<int, unordered_set<int>> ret;

	for (int i = 0; i < training_sequencde_database_sz_; ++i) {
		for (const auto& item : training_sequence_database_[i].sequence) {
			ret[item].insert(i);
		}
	}
	return ret;
}

void PatternMiner::RFGrowth(vector<int> pattern, unordered_set<PI::PatternInstance, boost::hash<PI::PatternInstance>> pi_set) {
	/*
		Find all possible candidates events for each sequence.
		S = [A x x x B o o o o o]
		'o' s are possible condidates
	*/
	unordered_set<int> candidate_events;
	for (const auto &entry : pi_set) {
		int id = entry.sid;
		vector<int> &S = training_sequence_database_[id].sequence;
		int sz = static_cast<int>(S.size());
		for (int i = entry.r + 1; i < sz; ++i) {
			candidate_events.insert(S[i]);
		}
	}

	vector<double> pgaps;
	int pattern_length = static_cast<int>(pattern.size());

	if (debug_) {
		printf("Pattern: ");
		for (int i = 0; i < pattern.size(); ++i) printf("%d ", pattern[i]);
		printf("\n");
		/*for (const auto &entry : pi_set) {
			printf("[%d %d]\n", entry.l, entry.r);
		}*/
	}

	for (const auto &e : candidate_events) {
		pattern.push_back(e); // P+ <- appending event e to pattern P
		unordered_set<PI::PatternInstance, boost::hash<PI::PatternInstance>> cur_pi_set = Grow(e, pi_set); // Grow()
		
		vector<double> pattern_gaps;
		pattern_gaps.resize(pattern_length);
		double max_upper = 0.0;
		double support = 0.0;
		pair<double, double> tmp;

		if (option_ == 1) {
			unordered_map<int, pair<vector<double>, double>> compact_gaps;
			unordered_map<int, pair<double, double>> compact_weights;
			map<int, double> upper_bound_map;

			unordered_map<int, int> max_length_map;
			double best_slice_sum = 100000000.0;

			// for every "pattern instance" in the "current pattern instance set"
			for (const auto& pattern_instance : cur_pi_set) {
				const int seq_id = pattern_instance.sid;
				const int ext_length = training_sequence_database_[seq_id].size() - pattern_instance.r;

				// tracking max_length
				if (max_length_map.find(seq_id) == max_length_map.end()) max_length_map[seq_id] = ext_length;
				else max_length_map[seq_id] = std::max(max_length_map[seq_id], ext_length);

				double weight = 0.0;

				// copy to pattern instance vector "pi"
				vector<int> pi;
				for (int i = pattern_instance.l; i <= pattern_instance.r; ++i) pi.push_back(training_sequence_database_[seq_id][i]);
				int pisz = static_cast<int>(pi.size());

				// receny measure calculation
				vector<double> gaps;
				int pre = pisz - 1;
				int cur = pattern_length - 1;
				for (int i = pisz - 2; i >= 0; --i) {
					if (pi[i] == pattern[cur]) {
						if (cur == 0 && i != 0) continue;
						gaps.push_back(static_cast<double>(pre - i));
						pre = i;
						--cur;
					}
				}
				std::reverse(gaps.begin(), gaps.end());

				int slice_sum = 1;
				int cur_val = 1;
				for (int i = gaps.size() - 1; i >= 0; --i) {
					cur_val += gaps[i];
					slice_sum += cur_val;
				}

				best_slice_sum = std::min(best_slice_sum, (double)slice_sum);

				compact_weights[pattern_instance.sid].first += (((double)pattern_length + 2.0)*((double)pattern_length + 1.0) / 2.0)/(double)slice_sum;
				compact_weights[pattern_instance.sid].second++;
				
				if (compact_gaps[pattern_instance.sid].second == 0) {
					compact_gaps[pattern_instance.sid].first.resize(pattern_length);
				}
				int sz = static_cast<int>(gaps.size());
				for (int i = 0; i < sz; ++i) compact_gaps[pattern_instance.sid].first[i] += gaps[i];
				compact_gaps[pattern_instance.sid].second++;

				double pg = 0.0;
				for (int i = 0; i < sz; ++i) pg += gaps[i];
				pgaps.push_back(pg);
			} // end of iterating all pattern instances


			for (auto &entry : compact_gaps) {
				int sz = entry.second.first.size();
				for (int i = 0; i < sz; ++i) {
					pattern_gaps[i] += (entry.second.first[i] / entry.second.second);
				}
			}

			for (auto &entry : pattern_gaps) entry /= static_cast<double>(compact_gaps.size());

			for (auto &entry : compact_weights) {
				entry.second.first /= static_cast<double>(entry.second.second);
				support += entry.second.first;
			}
			
			support /= static_cast<double>(training_sequencde_database_sz_);

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
			double P = (double)(pattern_length + 1);
			for (const auto &entry : upper_bound_map) {
				double K = (double)(-entry.first);
				double N = entry.second;
				double V = best_slice_sum;
				double a = (1.0 + P + K)*(P + K) / 2.0;
				double b = K*P;
				double c = (1.0 + K)*K / 2.0;
				double single_best_val = a / (V + b + c);
				max_upper = std::max(max_upper, single_best_val*N);
			}
			max_upper /= static_cast<double>(training_sequencde_database_sz_);
		}
		else {
			tmp = CalculateWeight(cur_pi_set, pattern_length + 1, option_);
			support = tmp.first / static_cast<double>(training_sequencde_database_sz_);
		}

		double frequency_support = CalculateWeight(cur_pi_set, pattern_length + 1, 3).first / static_cast<double>(training_sequencde_database_sz_);
		
		if (support >= threshold_) {
			if (option_ == 1) sequential_patterns_.push_back(Pattern(option_, pattern, support, frequency_support, pattern_gaps));
			else if (option_ == 2) {
				vector<double> avglen;
				avglen.push_back(tmp.second);
				sequential_patterns_.push_back(Pattern(option_, pattern, support, frequency_support, avglen));
			}
			else if (option_ == 3) sequential_patterns_.push_back(Pattern(option_, pattern, support, frequency_support));
			RFGrowth(pattern, cur_pi_set);
		}
		else{
			if(option_ == 1 && max_upper >= threshold_)RFGrowth(pattern, cur_pi_set);
			else if (option_ == 2 && UpperBound(cur_pi_set, pgaps, pattern_length + 1) >= threshold_) RFGrowth(pattern, cur_pi_set); // UpperBound function check
		}
		pattern.pop_back();
	}
}


unordered_set<PI::PatternInstance, boost::hash<PI::PatternInstance>> PatternMiner::Grow(int e, unordered_set<PI::PatternInstance, boost::hash<PI::PatternInstance>> pi_set) {
	unordered_set<int> sid_set;
	unordered_map<int, pair<int, int>> prev_tracker;
	for (const auto &PI : pi_set) sid_set.insert(PI.sid);
	for (const auto &sid : sid_set) {
		prev_tracker[sid].first = -1;
		prev_tracker[sid].second = -1;
	}
	set<PI::PatternInstance> sorted_pi_set;
	unordered_map<int, vector<pair<int, int>>> tracker;
	for (const auto &entry : pi_set) sorted_pi_set.insert(entry);
	for (const auto &entry : sorted_pi_set) {
		int id = entry.sid;
		int next_e = -1;

		const vector<int> &S = training_sequence_database_[id].sequence;
		int sz = static_cast<int>(S.size());
		for (int i = entry.r + 1; i < sz; ++i) {
			if (S[i] == e) {
				next_e = i;
				break;
			}
		}

		if (next_e != -1) {
			int curs = entry.l;
			int &sprev = prev_tracker[id].first;
			int &eprev = prev_tracker[id].second;
			if (sprev <= curs && next_e <= eprev) {
				tracker[id].pop_back();
			}
			tracker[id].push_back({ curs, next_e });
			sprev = curs;
			eprev = next_e;
		}
	}
	unordered_set<PI::PatternInstance, boost::hash<PI::PatternInstance>> ret;
	for (const auto &entry : tracker) {
		for (const auto &range_instance : entry.second) {
			ret.insert(PI::PatternInstance(entry.first, range_instance.first, range_instance.second));
		}
	}
	return ret;
}

pair<double, double> PatternMiner::CalculateWeight(unordered_set<PI::PatternInstance, boost::hash<PI::PatternInstance>> pi_set, int pattern_length, int option) {
	if (option == 2) {
		unordered_map<int, pair<double, double>> tracker;
		unordered_map<int, pair<double, double>> length_tracker;

		for (const auto & entry : pi_set) {
			int id = entry.sid;
			double w = static_cast<double>(pattern_length) / static_cast<double>(entry.r - entry.l + 1);
			tracker[id].first += w;
			tracker[id].second++;
			length_tracker[id].first += static_cast<double>(entry.r - entry.l + 1);
			length_tracker[id].second++;
		}
		double ret = 0.0;
		double len = 0.0;
		for (const auto & entry : tracker) {
			ret += (entry.second.first / entry.second.second);
		}
		for (const auto & entry : length_tracker) {
			len += (entry.second.first / entry.second.second);
		}
		return { ret, len/static_cast<double>(length_tracker.size())};
	}
	else if (option == 3) {
		unordered_set<int> vis;
		for (const auto & entry : pi_set) {
			int id = entry.sid;
			if (vis.find(id) != vis.end()) continue;
			vis.insert(id);
		}
		return { static_cast<double>(vis.size()), static_cast<double>(vis.size()) };
	}
	else {
		printf("Error in option CalculateCompactSupport()\n");
		getchar();
		exit(-1);
	}
}


double PatternMiner::UpperBound(unordered_set<PI::PatternInstance, boost::hash<PI::PatternInstance>> pi_set, vector<double> pgaps, int pattern_length) {
	unordered_map<int, pair<PI::PatternInstance, int>> mp;

	int cnt = 0;
	for (auto it = pi_set.begin(); it != pi_set.end(); ++it, ++cnt) {
		int sid = it->sid;
		int l = it->l;
		int r = it->r;
		if (mp.find(sid) == mp.end()) {
			mp[sid] = { *it, cnt };
		}
		else if (l < mp[sid].first.l) {
			mp[sid] = { *it, cnt };
		}
	}

	double F = 0.0;
	
	if (option_ == 1) {
		int max_len = 0;

		for (auto it = mp.begin(); it != mp.end(); ++it) {
			int sid = it->second.first.sid;
			int l = it->second.first.l;
			int r = it->second.first.r;
			int idx = it->second.second;
			int sz = static_cast<int>(training_sequence_database_[sid].size());
			int tl = sz - r;
			max_len = std::max(max_len, tl);
			F += (pgaps[idx] + static_cast<double>(tl)) / static_cast<double>(pattern_length - 1 + tl);
		}
		F /= static_cast<double>(training_sequencde_database_sz_);

		int frequency = static_cast<int>(CalculateWeight(pi_set, pattern_length, 3).first);
		return Potential(pattern_length, frequency, max_len);
	}
	else if(option_ == 2){
		for (auto it = mp.begin(); it != mp.end(); ++it) {
			int sid = it->second.first.sid;
			int l = it->second.first.l;
			int r = it->second.first.r;
			int sz = static_cast<int>(training_sequence_database_[sid].size());
			int tl = sz - r;
			F += static_cast<double>(pattern_length - 1 + tl) / static_cast<double>(r - l + tl);
		}
		F /= static_cast<double>(training_sequencde_database_sz_);
	}
	else {
		printf("Error in option UpperBound()\n");
		getchar();
		exit(-1);
	}
	return F;
}

double PatternMiner::Potential(int pattern_length, int frequency, int max_len) {
	double l = static_cast<double>(pattern_length);
	double k = static_cast<double>(max_len);
	double n = static_cast<double>(frequency);
	return ((l + k)*threshold_ - k * n / static_cast<double>(training_sequencde_database_sz_)) / l;
}

void PatternMiner::SetDebug() {
	debug_ = 1;
}