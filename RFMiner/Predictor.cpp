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
		Rule r = CreateRule(top_patterns_[i], query_sequence, option, repeat);
		if (r.rule_type == EMPTY) continue;
		rule_list_.push_back(r);
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
	vector<double> max_score_gaps;
	
	//assert(final_pi_set.size() == 1);
	for (const auto &pi : final_pi_set) {

		if (option == RECENCY) {
			double transition_difference = 0.0;
			double transition_score = -1.0;
			vector<double> gaps;
			double weight = -1.0;	

			int j = pattern_sz - 2;
			for (int i = pi.r; i >= pi.l; --i) {
				if (query_sequence[i] == pattern_sequence[j]) {
					if (j == 0) assert(i == pi.l);
					gaps.push_back(static_cast<double>(query_sz - i + 1));
					--j;
				}
			}
			std::reverse(gaps.begin(), gaps.end());

			vector<double> gap_matches;
			int gaps_sz = static_cast<int>(gaps.size());
			const vector<double>& pattern_gaps = pattern.gap_sequence;
			for (int i = 0; i < gaps_sz; ++i) {
				gap_matches.push_back(abs(gaps[i] - pattern_gaps[i]));
			}
			int gap_matches_sz = static_cast<int>(gap_matches.size());
			//transition_difference = 1000000.0;
			for (int i = 0; i < gap_matches_sz; ++i) {
				transition_difference += gap_matches[i];
				//transition_difference = std::min(transition_difference, gap_matches[i]);
			}
			// this is newly added feature
			transition_difference /= static_cast<double>(gap_matches.size());
			transition_score = 1.0 / (transition_difference + 1.0) * pattern.confidence;

			if (max_transition_score < transition_score) {
				max_transition_score = transition_score;
				max_score_gaps = gap_matches;
			}
		}
		else if (option == COMPACTNESS) {
			double ratio_difference = abs(pattern.gap_sequence.front() - static_cast<double>((query_sz - pi.l + 1)));
			double ratio_score = 1.0 / (ratio_difference + 1.0) * pattern.confidence;
			if (max_ratio_score < ratio_score) {
				max_ratio_score = ratio_score;
			}
		}
		else if (option == PRESENCE) {
			double frequency_score = (pattern.pattern.size() -1) * pattern.confidence;
			if (max_frequency_score < frequency_score) {
				max_frequency_score = frequency_score;
			}
		}
	}

	vector<int> precedent;
	vector<double> precedent_gap;
	vector<int> consequent;
	vector<double> consequent_gap;
	vector<double> query_gap;

	// Push [1:|P| - 1] to precedent
	for (int i = 0; i < pattern_sz - 1; ++i) precedent.push_back(pattern_sequence[i]);
	// Push the last element of P to consequent
	consequent.push_back(pattern_sequence.back());

	return Rule(precedent, consequent, max_score_gaps, max_transition_score, max_ratio_score, max_frequency_score, pattern.frequency, option);
}

vector<PatternInstance> Predictor::CompactInstances(int e, vector<PatternInstance> pi_set, const vector<int> &query_sequence, vector<PatternInstance> &final_pi_set, int pattern_offset, int pattern_sz) {
	vector<PatternInstance> ret;
	int pi_set_sz = static_cast<int>(pi_set.size());
	for (int i = 0; i < pi_set_sz; ++i) {
		PatternInstance &pi = pi_set[i];
		const int id = pi.sid;

		int sz = static_cast<int>(query_sequence.size());
		if (i + 1 < pi_set_sz) {
			sz = std::min(sz, pi_set[i + 1].r + 1);
		}

		bool flag = false;
		for (int j = pi.r + 1; j < sz; ++j) {
			if (query_sequence[j] == e) {
				flag = true;
				ret.push_back(PatternInstance(id, pi.l, j));
				break;
			}
		}

		// if we have looked for the last element in the pattern, but failed to find it,
		// then the query string contains the prefix[1:|P| - 1] of the pattern, and not the whole pattern [1: |P|]
		if (pattern_offset == pattern_sz - 1 && !flag) {
			final_pi_set.push_back(pi);
		}
	}
	return ret;
}

void Predictor::SetDebug() {
	debug_ = 1;
}
