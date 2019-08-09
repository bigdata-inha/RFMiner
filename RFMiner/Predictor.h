#ifndef PREDICTOR_H_
#define PREDICTOR_H_

#include "StructTypes.h"

struct Rule {
	vector<int> precedent, consequent;
	vector<double> gaps;
	int rule_type;
	double transition_score, ratio_score, frequency_score, frequency;

	Rule(int a_rule_type)
		: transition_score(0.0), ratio_score(0.0), frequency(0.0), rule_type(a_rule_type) {}

	Rule(vector<int> a_precedent, vector<int> a_consequent, vector<double> a_gaps, double a_transition_score, double a_ratio_score, double a_frequency_score, double a_frequency, int a_type)
		: precedent(a_precedent), consequent(a_consequent), gaps(a_gaps), transition_score(a_transition_score), ratio_score(a_ratio_score), frequency_score(a_frequency_score), frequency(a_frequency), rule_type(a_type) {}

	void Print() {
		if (rule_type == 1) printf("T");
		else if (rule_type == 2) printf("R");
		else if (rule_type == 3) printf("F");
		printf("{ ");
		for (const auto &item : precedent) printf("%d ", item);
		printf("} => { ");
		for (const auto &item : consequent) printf("%d ", item);
		printf("}(score: ");
		if (rule_type == 1) printf("[%lf], frequency: [%lf])\n", transition_score, frequency);
		else if (rule_type == 2) printf("[%lf], frequency: [%lf])\n", ratio_score, frequency);
		else if (rule_type == 3)printf("[%lf], frequency: [%lf])\n", frequency_score, frequency);
		if (rule_type == 1) {
			printf("[ ");
			for (const auto &w : gaps) printf("%.2lf ", w);
			printf("]\n");
		}
	}

	bool operator<(const Rule &rhs) const {
		if (rule_type == 1)
			return transition_score > rhs.transition_score;
		else if (rule_type == 2)
			return ratio_score > rhs.ratio_score;
		else if (rule_type == 3)
			return frequency_score > rhs.frequency_score;
	}
};

class Predictor {
public:

	/* Constructor, Destructor */
	Predictor(vector<Pattern> interesting_patterns);
	~Predictor();

	/* Sets number of patterns used for recommendation */
	void SetTopPatternNumber(int top_pattern_number);

	/* Calculates prediction quality of each pattern, given a query sequence */
	void GenerateRules(const vector<int> &query_sequence, int option, const int &repeat);
	
	/* Return topK events for prediction */
	vector<int> ReturnRulePrediction(const int &top_k, const int &repeat);

	/* Print topK rules for given query sequence, should be called after GenerateRules() */
	void PrintRuleList(int topK);

	/* Resets recommendation score calculation done by GenerateRules() */
	void ClearRuleList();

	void SetDebug();
private:
	vector<Rule> rule_list_;

	vector<Pattern> interesting_patterns_;
	int interesting_patterns_sz;
	vector<Pattern> top_patterns_;
	int top_patterns_sz_;

	int debug_;

	Rule CreateRule(const Pattern &pattern, const vector<int> &query_sequence, int option, const int &repeat);
	vector<PatternInstance> CompactInstances(int e, vector<PatternInstance> pi_set, const vector<int> &query_sequence, vector<PatternInstance> &final_pi_set, int pattern_offset, int pattern_sz);
};

#endif // !PREDICTOR_H_