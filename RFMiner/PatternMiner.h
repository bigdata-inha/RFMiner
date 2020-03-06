#ifndef PATTERNMINER_H_
#define PATTERNMINER_H_

#include "StructTypes.h"

// Iterate over the database sequences and extracts interesting sequential patterns.
//		1 - Recency
//		2 - Compactness
//		3 - Presense(Frequency)

class PatternMiner {
public:
	PatternMiner();
	~PatternMiner();
	void ClearPatterns();
	void LoadTrainingDatabase(vector<Sequence> database);
	void LoadTestDatabase(vector<Sequence> database);
	void Run(double init_threshold, double threshold, int option);
	vector<Pattern> GetSequentialPatterns();
	void WritePatternFile(string filename);
	void WriteQueryFile(string filename);

	int cnt;
	void SetDebug();
	int efficient_upperbound;
private:
	double init_threshold_;;
	double threshold_;
	unordered_set<int> candidate_events_;

	int training_sequencde_database_sz_;
	double double_training_db_sz;
	vector<Sequence> training_sequence_database_;

	int test_sequence_database_sz_;
	vector<Sequence> test_sequence_database_;

	vector<Pattern> sequential_patterns_;

	int option_;
	int debug_;

	void RFGrowth(vector<int> pattern, vector<PatternInstance> pi_set);
	unordered_map<int, vector<PatternInstance>> GrowNew(vector<PatternInstance> pi_set, int last_event, unordered_map<int, unordered_map<int, int>> &max_length_map, unordered_map<int, double> &min_instance, unordered_map<int, unordered_set<int>> &inverted_list);
	double UpperBound(double plus_pattern_len, unordered_map<int, int > max_length_map, double min_instance);

};

#endif // !PATTERNMINER_H_