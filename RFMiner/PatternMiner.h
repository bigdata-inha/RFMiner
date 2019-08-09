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

	void SetDebug();
private:
	double init_threshold_;;
	double threshold_;

	int training_sequencde_database_sz_;
	vector<Sequence> training_sequence_database_;

	int test_sequence_database_sz_;
	vector<Sequence> test_sequence_database_;

	vector<Pattern> sequential_patterns_;

	int option_;

	int debug_;

	unordered_map<int, unordered_set<int>> ScanCountSingleItemsInit();
	void RFGrowth(vector<int> pattern, vector<PatternInstance> pi_set);
	vector<PatternInstance> Grow(int e, vector<PatternInstance> pi_set);
	pair<double, double> CalculateWeight(vector<PatternInstance> pi_set, int pattern_length, int option);
	double UpperBound(vector<PatternInstance> pi_set, vector<double> pgaps, int pattern_length);
	double Potential(int pattern_length, int frequency, int max_len);
};

#endif // !PATTERNMINER_H_