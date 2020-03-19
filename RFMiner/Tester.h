#include <vector>
#include <algorithm>

#include "Database.h"
#include "PrefixSpan.h"

#pragma once
using std::vector;
using std::max;
using std::reverse;
using std::to_string;
using std::queue;

struct Measure {
	double accuracy;
	double precision;
	double recall;
	double ndcg;
	double r_precision;
	double tf_idf;
	double relevance_score;
	double coverage;
	int total_queries;
	int matched_queries;
	char buffer[1000];

	Measure()
		:accuracy(0.0), precision(0.0), recall(0.0), ndcg(0.0), r_precision(0.0), tf_idf(0.0), relevance_score(0.0), total_queries(0), matched_queries(0) {}

	void Calculate() {
		double total = static_cast<double>(total_queries);
		double match = static_cast<double>(matched_queries);
		accuracy = 100.0*accuracy / total;
		precision = 100.0*precision / total;
		r_precision = 100.0*r_precision / total;
		ndcg = 100.0*ndcg / total;
		coverage = 100.0*match / total;
		recall = 100.0*recall / total;
		tf_idf = tf_idf / total;
	}

	void PrintStat() {
		sprintf(buffer, "Accuracy: [%.2lf] Precision: [%.2lf], Rprecision: [%.2lf], NDCG: [%.2lf], Coverage: [%.1lf], Recall: [%.2lf], TFIDF: [%.6lf]\n", 
			accuracy, precision, r_precision, ndcg, coverage, recall, tf_idf);
		cout << buffer;
	}
};

struct PatternSetAnalysis {
	vector<Pattern> patterns;
	vector<Pattern> active_patterns;
	unordered_set<int> unique_items;
	double avg_pattern_length;
	int number_unique_items;
	int pattern_num_lim;

	PatternSetAnalysis(vector<Pattern> a_patterns)
		:patterns(a_patterns), pattern_num_lim(1000000) {
		// patterns are sorted in descending order according to their score
		sort(patterns.begin(), patterns.end());
		avg_pattern_length = 0.0;
		number_unique_items = 0;
		active_patterns = patterns;
	}

	void Count() {
		avg_pattern_length = 0.0;
		number_unique_items = 0;
		for (const auto& p : active_patterns) {
			double len = static_cast<double>(p.pattern.size());
			avg_pattern_length += len;
			for (const auto &item : p.pattern) {
				unique_items.insert(item);
			}
		}
		avg_pattern_length /= static_cast<double>(active_patterns.size());
		number_unique_items = static_cast<int>(unique_items.size());
	}

	void set_pattern_num_lim(int a_pattern_num_lim) {
		if (a_pattern_num_lim == -1) pattern_num_lim = static_cast<int>(patterns.size());
		else pattern_num_lim = a_pattern_num_lim;
		int sz = std::min(static_cast<int>(patterns.size()), pattern_num_lim);
		vector<Pattern> vec;
		for (int i = 0; i < sz; ++i) {
			vec.push_back(patterns[i]);
		}
		active_patterns = vec;
	}

	void PrintTopPatterns(int topK) {
		for (int i = 0; i < topK; ++i) {
			patterns[i].PrintPattern();
		}
	}

	double GetLastThreshold() {
		return active_patterns.back().interestingness;
	}

	double get_avg_pattern_len() { return avg_pattern_length; }
	int get_number_unique_items() { return number_unique_items; }

};


class Tester {
public:
	Tester() {};
	~Tester() {};

	void test_loadable(int fold_size, double split_ratio, int top_k, int patter_num_lim, string path, int repeated_events, int total_experiment, int offset);
	void fold_divider(int fold_size, double transition_ratio_init_threshold, double transition_threshold, double ratio_threshold, double frequency_threshold, const string& path, Database &database);

	double Accuracy(const vector<int> &rec, const vector<int> &ans);
	double Rprecision(const vector<int> &rec, const vector<int> &ans);
	double NDCG(const int sequence_id, const vector<int> &rec, const vector<int> & ans);
	double Precision(const vector<int> &rec, const vector<int> &ans);
	double Recall(const vector<int> &rec, const vector<int> &ans);
	double TFIDF(const int sequence_id, const vector<int> &rec, const vector<int> &ans);

	void set_result_filename(const string &filename) { result_filename = filename; }
	vector<Pattern> LoadFrequentPatterns(const string &filename, const int option);

	void LoadDatabase(Database &database) {
		database_ = database;
	}

	void WriteTopPatternInfo(int fold_size, string path, vector<int> pattern_numbers);

	void ExperimentTopK(int fold_size, string path, int top_k, int repeated_events, vector<int> pattern_numbers, vector<double> split_ratio);

	void ExecutionTimeTest(int fold_size, string path, string dataset_path, vector<int> pattern_numbers);

	vector<double> KthThreshold(int fold_size, string path, int top_patterns);
	vector<double> KthRunTime(string path, int top_patterns, vector<double> thresholds, Database &database);

	double NaiveCal(string path, int top_patterns, double thresholds, Database &database);

	void ExperimentTimeNaive(int fold_size, string path, string dataset_path, vector<int> pattern_numbers);

	void CaseStudy(string path);


private:
	vector<int> Q, S;
	Database database_;
	
	const double truth_split = 0.1;

	string result_filename;

	vector<Measure> measures_[3];

	vector<Sequence> LoadTestSequences(const string &file);

};