#include "Predictor.h"

Predictor::Predictor(vector<Pattern> interesting_patterns)
	: interesting_patterns_(interesting_patterns) {
	interesting_patterns_sz = static_cast<int>(interesting_patterns_.size());
	top_patterns_ = interesting_patterns_;
	top_patterns_sz_ = interesting_patterns_sz;
	debug_ = 0;
}

Predictor::~Predictor() {}

void Predictor::SetTopPatternNumber(int pattern_num_lim) {
	top_patterns_sz_ = pattern_num_lim;
	int sz = std::min(interesting_patterns_sz, top_patterns_sz_);
	sort(interesting_patterns_.begin(), interesting_patterns_.end());
	vector<Pattern> vec;
	for (int i = 0; i < sz; ++i) {
		vec.push_back(interesting_patterns_[i]);
	}
	top_patterns_ = vec;
	top_patterns_sz_ = sz;
}

void Predictor::GenerateRules(const vector<int> &query_sequence, int option, const int &repeat) {
	for (int i = 0; i < top_patterns_sz_; ++i) {
		Rule tmp = CreateRule(top_patterns_[i], query_sequence, option, repeat);
		
		/* A valid rule_type is either 1, 2, or 3 */
		if (tmp.rule_type == 1 || tmp.rule_type == 2 || tmp.rule_type == 3) {
			rule_list_.push_back(tmp);
		}
	}
	if (static_cast<int>(rule_list_.size()) > top_patterns_sz_) {
		printf("rule list size error!");
		getchar();
		exit(1);
	}
	sort(rule_list_.begin(), rule_list_.end());
}

vector<int> Predictor::ReturnRulePrediction(const int &top_k, const int &repeat) {
	unordered_set<int> vis;
	vector<int> ret;
	int rsz = static_cast<int>(rule_list_.size());
	for (int i = 0; i < rsz; ++i) {
		for (const auto &item : rule_list_[i].consequent) {
			if (repeat == 0 && vis.find(item) != vis.end()) continue;
			if (repeat == 0) vis.insert(item);
			ret.push_back(item);
			if (static_cast<int>(ret.size()) == top_k) break;
		}
		if (static_cast<int>(ret.size()) == top_k) break;
	}
	return ret;
}

void Predictor::ClearRuleList() {
	rule_list_.clear();
};

void Predictor::PrintRuleList(int topK) {
	int sz = static_cast<int>(rule_list_.size());
	for (int i = 0; i < sz; ++i) {
		if (i == topK) break;
		rule_list_[i].Print();
	}
}

Rule Predictor::CreateRule(const Pattern &pattern, const vector<int> &query_sequence, int option, const int &repeat) {

	const vector<int>& pattern_sequence = pattern.pattern;
	const int pattern_sz = static_cast<int>(pattern_sequence.size());
	const int query_sz = static_cast<int>(query_sequence.size());
	const vector<double>& pattern_gaps = pattern.gap_sequence;

	// Find minimal landmarks of a pattern P's prefix [1: |P| - 1] in query sequence
	vector<PatternInstance> pi_set;
	vector<PatternInstance> final_pi_set;

	/* Find all indices of query sequence, where the event at that index is the first event in the pattern sequence */
	for (int i = 0; i < query_sz; ++i) {
		if (query_sequence[i] == pattern_sequence[0]) pi_set.push_back(PatternInstance(0, i, i));
	}
	
	// If there is no match
	if (pi_set.empty()) return Rule(0);

	for (int i = 1; i < pattern_sz; ++i) {
		pi_set = CompactInstances(pattern_sequence[i], pi_set, query_sequence, final_pi_set, i, pattern_sz);
	}

	// If there is no match
	if (final_pi_set.empty()) return Rule(0);



	double max_transition_score = -1.0;
	double max_ratio_score = -1.0;
	double max_frequency_score = -1.0;
	int x = -1;
	vector<double> max_score_gaps;
	
	for (const auto &pi : final_pi_set) {
		const int &last_matched_pattern_offset = pattern_sz - 2;

		// P = <a, b, c, d>
		// Q = <a, x, x, b, c, x, x>

		if (debug_) printf("[%d %d] pattern_offset: %d\n", pi.l, pi.r, last_matched_pattern_offset);

		if (option == 1) {
			double transition_difference = 0.0;
			double transition_score = -1.0;
			vector<double> gaps;
			double weight = -1.0;	
		
			if(debug_) printf("unique\n");

			int pre = query_sz;
			int cur = last_matched_pattern_offset;
			for (int i = query_sz - 1; i >= 0 && cur >= 0; --i) {
				if (query_sequence[i] == pattern_sequence[cur]) {
					gaps.push_back(static_cast<double>(pre - i));
					pre = i;
					--cur;
				}
			}
			std::reverse(gaps.begin(), gaps.end());

			vector<double> gap_matches;
			int gaps_sz = static_cast<int>(gaps.size());
			for (int i = 0; i < gaps_sz; ++i) {
				gap_matches.push_back(abs(gaps[i] - pattern_gaps[i]));
			}
			int gap_matches_sz = static_cast<int>(gap_matches.size());
			for (int i = 0; i < gap_matches_sz; ++i) {
				transition_difference += gap_matches[i];
			}
			transition_score = 1.0 / (transition_difference + 1.0) * pattern.interestingness;

			if (max_transition_score < transition_score) {
				max_transition_score = transition_score;
				x = last_matched_pattern_offset;
				max_score_gaps = gap_matches;
			}
		}
		else if (option == 2) {
			double ratio_difference = abs(pattern.gap_sequence.front() - static_cast<double>((query_sz - pi.l)));
			double ratio_score = 1.0 / (ratio_difference + 1.0) * pattern.interestingness;
			if (max_ratio_score < ratio_score) {
				max_ratio_score = ratio_score;
				x = last_matched_pattern_offset;
			}
		}
		else if (option == 3) {
			double frequency_score = pattern.frequency;
			if (max_frequency_score < frequency_score) {
				max_frequency_score = frequency_score;
				x = last_matched_pattern_offset;
			}
		}
		else {
			printf("Error in option CreateRule()");
			getchar();
			exit(-1);
		}
	}

	vector<int> precedent;
	vector<double> precedent_gap;
	vector<int> consequent;
	vector<double> consequent_gap;
	vector<double> query_gap;

	for (int i = 0; i < pattern_sz; ++i) {
		if (i <= x) precedent.push_back(pattern_sequence[i]);
		else consequent.push_back(pattern_sequence[i]);
	}

	return Rule(precedent, consequent, max_score_gaps, max_transition_score, max_ratio_score, max_frequency_score, pattern.frequency, option);
}

vector<PatternInstance> Predictor::CompactInstances(int e, vector<PatternInstance> pi_set, const vector<int> &query_sequence, vector<PatternInstance> &final_pi_set, int pattern_offset, int pattern_sz) {
	
	pair<int, int> prev_tracker(-1, -1);
	set<PatternInstance> sorted_pi_set;
	vector<pair<pair<int, int>, int>> tracker;
	unordered_set<int> keep;

	// sort landmarks in ascending order of the first index
	for (const auto &entry : pi_set) sorted_pi_set.insert(entry);

	int num = 0;
	for (const auto &entry : sorted_pi_set) {
		int next_e = -1;

		int sz = static_cast<int>(query_sequence.size());
		for (int i = entry.r + 1; i < sz; ++i) {
			if (query_sequence[i] == e) {
				next_e = i;
				break;
			}
		}

		if (next_e != -1) {
			int curs = entry.l;
			int &sprev = prev_tracker.first;
			int &eprev = prev_tracker.second;
			if (sprev <= curs && next_e <= eprev) {
				keep.insert(tracker.back().second);
				tracker.pop_back();
			}
			tracker.push_back({{ curs, next_e }, num});
			sprev = curs;
			eprev = next_e;
		}
		else if(pattern_offset == pattern_sz - 1){
			final_pi_set.push_back(entry);
		}
		++num;
	}
	vector<PatternInstance> ret;

	for (const auto &range_instance : tracker) {
		ret.push_back(PatternInstance(0, range_instance.first.first, range_instance.first.second));
	}
	
	return ret;
}

void Predictor::SetDebug() {
	debug_ = 1;
}
