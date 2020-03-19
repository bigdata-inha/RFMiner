#include "Tester.h"
#include "Predictor.h"
#include "PatternMiner.h"

template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
	os << "[";
	for (int i = 0; i < v.size(); ++i) {
		os << v[i];
		if (i != v.size() - 1)
			os << ", ";
	}
	os << "]";
	return os;
}

void Tester::fold_divider(int fold_size, double transition_ratio_init_threshold, double transition_threshold, double ratio_threshold, double frequency_threshold, const string& path, Database &database) {
	
	printf("Dividing dataset %d folds\n", fold_size);
	
	PatternMiner transition_version, ratio_version;
	PrefixSpan frequency_version;

	database.set_fold(fold_size);
	database.print_stats();

	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto mining_time_transition = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	auto mining_time_ratio = mining_time_transition;
	auto mining_time_frequency = mining_time_transition;

	int avg_transition_pattern_set_sz = 0;
	int avg_ratio_pattern_set_sz = 0;
	int avg_frequency_pattern_set_sz = 0;

	char buffer[1000];

	string config = "[" + to_string(fold_size) + ", " + to_string(transition_threshold) + ", " + to_string(ratio_threshold) + ", " + to_string(frequency_threshold) + "]";
	ofstream outfile(path + "/Experiment/Mining.txt");

	for (int i = 0; i < fold_size; ++i) {

		transition_version.ClearPatterns();
		ratio_version.ClearPatterns();
		frequency_version.ClearTraining();

		transition_version.LoadTrainingDatabase(database.get_train_db(i));
		ratio_version.LoadTrainingDatabase(database.get_train_db(i));
		frequency_version.LoadTrainDatabase(database.get_train_db(i));
		
		start = std::chrono::high_resolution_clock::now();
		transition_version.Run(transition_ratio_init_threshold, transition_threshold, 1);
		stop = std::chrono::high_resolution_clock::now();
		mining_time_transition += std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

		start = std::chrono::high_resolution_clock::now();
		ratio_version.Run(transition_ratio_init_threshold, ratio_threshold, 2);
		stop = std::chrono::high_resolution_clock::now();
		mining_time_ratio += std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

		start = std::chrono::high_resolution_clock::now();
		frequency_version.Run(frequency_threshold);
		stop = std::chrono::high_resolution_clock::now();
		mining_time_frequency += std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

		avg_transition_pattern_set_sz += static_cast<int>(transition_version.GetSequentialPatterns().size());
		avg_ratio_pattern_set_sz += static_cast<int>(ratio_version.GetSequentialPatterns().size());
		avg_frequency_pattern_set_sz += static_cast<int>(frequency_version.GetFrequentPatterns().size());

		auto test_db = database.get_test_db(i);
		string filenameA = path + "foldA" + to_string(i + 1);
		string filenameB = path + "foldB" + to_string(i + 1);
		string filenameC = path + "foldC" + to_string(i + 1);
		string testfilename = path + "test" + to_string(i + 1);

		transition_version.WritePatternFile(filenameA);
		ratio_version.WritePatternFile(filenameB);
		frequency_version.WriteFile(filenameC);

		frequency_version.LoadTestDatabase(test_db);
		frequency_version.WriteQueryFile(testfilename);
	}

	double fs = static_cast<double>(fold_size);

	mining_time_transition /= fs;
	mining_time_ratio /= fs;
	mining_time_frequency /= fs;

	avg_transition_pattern_set_sz = avg_transition_pattern_set_sz/ fold_size;
	avg_ratio_pattern_set_sz = avg_ratio_pattern_set_sz/ fold_size;
	avg_frequency_pattern_set_sz = avg_frequency_pattern_set_sz/ fold_size;

	sprintf(buffer, "Threshold: %lf\tk-fold: %d\t\n", transition_threshold, fold_size);
	outfile << buffer;
	sprintf(buffer, "\tNumber of patterns\tExecution Time(ms)\tDP function called\n");
	outfile << buffer;
	sprintf(buffer, "Transition\t%d\t%lf\n", avg_transition_pattern_set_sz, static_cast<double>(mining_time_transition.count()));
	outfile << buffer;
	sprintf(buffer, "Ratio\t%d\t%lf\n", avg_ratio_pattern_set_sz, static_cast<double>(mining_time_ratio.count()));
	outfile << buffer;
	sprintf(buffer, "Frequency\t%d\t%lf\n", avg_frequency_pattern_set_sz, static_cast<double>(mining_time_frequency.count()));
	outfile << buffer;

	sprintf(buffer, "Transition: [#%d] [%lfms]\n", avg_transition_pattern_set_sz, static_cast<double>(mining_time_transition.count()));
	cout << buffer;
	sprintf(buffer, "Ratio: [#%d] [%lfms]\n", avg_ratio_pattern_set_sz, static_cast<double>(mining_time_ratio.count()));
	cout << buffer;
	sprintf(buffer, "Frequency: [#%d] [%lfms]\n", avg_frequency_pattern_set_sz, static_cast<double>(mining_time_frequency.count()));
	cout << buffer;
}

vector<Pattern> Tester::LoadFrequentPatterns(const string &filename, const int option) {
	ifstream infile(filename);
	vector<Pattern> ret;
	int x;
	double y;
	double frequency;
	double confidence;
	int pattern_id = 0;
	vector<int> vec;
	while (infile >> x) {
		if (x == -1) {
			infile >> y;
			if(option != 3) infile >> frequency;
			infile >> confidence;
			vector<double> tmp_weights;
			double weight;
			if (option == 1) {
				for (int i = 0; i < vec.size() - 1; ++i) {
					infile >> weight;
					tmp_weights.push_back(weight);
				}
			}
			else if (option == 2) {
				infile >> weight;
				tmp_weights.push_back(weight);
			}
			if(option == 1)ret.push_back(Pattern(option, vec, y, frequency, tmp_weights, confidence));
			else if (option == 2) ret.push_back(Pattern(option, vec, y, frequency, tmp_weights, confidence));
			else if (option == 3) ret.push_back(Pattern(option, vec, y, y, confidence));
			ret.back().pattern_id = pattern_id++;
			vec.clear();
			continue;
		}
		vec.push_back(x);
	}
	return ret;
}

vector<Sequence> Tester::LoadTestSequences(const string &file) {
	ifstream infile(file);
	vector<Sequence> ret;
	int sequence_id;
	while (infile >> sequence_id) {
		int x;
		vector<int> sequence;
		while (infile >> x) {
			if (x == -1) {
				ret.push_back(Sequence(sequence_id, sequence));
				break;
			}
			sequence.push_back(x);
		}
	}
	return ret;
}

void Tester::test_loadable(int fold_size, double split_ratio, int top_k, int pattern_num_lim, string path, int repeated_events, int total_experiment, int offset) {
	int mismatch_because_insufficient_K = 0;
	int dump_query = 0;
	int no_rep = 0;
	double success = 0;
	double win = 0.0;
	double lose = 0.0;
	double avg_query_length = 0.0;
	double avg_query_set_sz = 0.0;
	double transition_avg_pattern_length = 0.0;
	double ratio_avg_pattern_length = 0.0;
	double frequency_avg_pattern_length = 0.0;
	double transition_avg_unique_item = 0.0;
	double ratio_avg_unique_item = 0.0;
	double frequency_avg_unique_item = 0.0;

	Measure transition_measure;
	Measure ratio_measure;
	Measure frequency_measure;

	char buffer[1000];
	vector<int> query_sequence, ground_truth_sequence , transition_retrieved, ratio_retrieved, frequency_retrieved;

	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto query_time_transition = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	auto query_time_ratio = query_time_transition;
	auto query_time_frequency = query_time_transition;

	string output_filename = result_filename;
	ofstream outfile(output_filename);

	// Each fold is different
	for (int i = 0; i < fold_size; ++i) {
		sprintf(buffer, "\n%d(th) testing...\n", i + 1);
		cout << buffer;
		outfile << buffer;

		string fileA = path + "foldA" + to_string(i + 1);
		string fileB = path + "foldB" + to_string(i + 1);
		string fileC = path + "foldC" + to_string(i + 1);

		string testfile = path + "test" + to_string(i + 1);

		cout << fileA << " " << fileB << " " << fileC << " " << testfile << "\n";

		auto transition_pattern_set = LoadFrequentPatterns(fileA, 1);
		auto ratio_pattern_set = LoadFrequentPatterns(fileB, 2);
		auto frequency_pattern_set = LoadFrequentPatterns(fileC, 3);

		auto query_database = LoadTestSequences(testfile);

		PatternSetAnalysis TM(transition_pattern_set);
		PatternSetAnalysis RM(ratio_pattern_set);
		PatternSetAnalysis FM(frequency_pattern_set);
		if (pattern_num_lim != -1) {
			TM.set_pattern_num_lim(pattern_num_lim);
			RM.set_pattern_num_lim(pattern_num_lim);
			FM.set_pattern_num_lim(pattern_num_lim);
		}

		TM.Count();
		RM.Count();
		FM.Count();

		avg_query_set_sz += static_cast<double>(query_database.size());
		transition_avg_pattern_length += TM.avg_pattern_length;
		ratio_avg_pattern_length += RM.avg_pattern_length;
		frequency_avg_pattern_length += FM.avg_pattern_length;
		transition_avg_unique_item += static_cast<double>(TM.number_unique_items);
		ratio_avg_unique_item += static_cast<double>(RM.number_unique_items);
		frequency_avg_unique_item += static_cast<double>(FM.number_unique_items);

		sprintf(buffer, "transition pattern set: %d, ratio pattern set: %d, frequency pattern set: %d, query set: %d\n", transition_pattern_set.size(), ratio_pattern_set.size(), frequency_pattern_set.size(), query_database.size());
		cout << buffer;
		outfile << buffer;
		sprintf(buffer, "Transition: avg pattern length: %lf unique items: %d\n", TM.avg_pattern_length, TM.number_unique_items);
		cout << buffer;
		outfile << buffer;
		sprintf(buffer, "Ratio: avg pattern length: %lf unique items: %d\n", RM.avg_pattern_length, RM.number_unique_items);
		cout << buffer;
		outfile << buffer;
		sprintf(buffer, "Frequency: avg pattern length: %lf unique items: %d\n", FM.avg_pattern_length, FM.number_unique_items);
		cout << buffer;
		outfile << buffer;

		Predictor transition_predictor(transition_pattern_set);
		Predictor ratio_predictor(ratio_pattern_set);
		Predictor frequency_predictor(frequency_pattern_set);
		if (pattern_num_lim != -1) {
			transition_predictor.SetTopPatternNumber(pattern_num_lim);
			ratio_predictor.SetTopPatternNumber(pattern_num_lim);
			frequency_predictor.SetTopPatternNumber(pattern_num_lim);
		}

		int query_database_sz = static_cast<int>(query_database.size());


		double valid_queries = 0.0;
		// queries that aren't suppose to be answered 
		double dump_queries = 0.0;

		for (int t = 0; t < query_database_sz; ++t) {
			int sequence_id = query_database[t].sequence_id;
			vector<int> &tmp = query_database[t].sequence;
			if (tmp.size() == 1) {
				dump_query++;
				dump_queries++;
				continue;
			}

			query_sequence.clear();
			ground_truth_sequence.clear();
			transition_predictor.ClearRuleList();
			ratio_predictor.ClearRuleList();
			frequency_predictor.ClearRuleList();
			
			transition_retrieved.clear();
			ratio_retrieved.clear();
			frequency_retrieved.clear();

			int sz = (int)((double)tmp.size()*(1.0 - truth_split));
			vector<int> tmp_seqeuence;

			for (int k = 0; k < tmp.size(); ++k) {
				if (k < sz) tmp_seqeuence.push_back(tmp[k]);
				else ground_truth_sequence.push_back(tmp[k]);
			}

			sz = (int)((double)tmp_seqeuence.size()*split_ratio);
			for (int k = tmp_seqeuence.size() - sz; k < tmp_seqeuence.size(); ++k) {
				query_sequence.push_back(tmp[k]);
			}

			if (query_sequence.size() == 0 || ground_truth_sequence.size() == 0) {
				dump_query++;
				dump_queries++;
				continue;
			}
			
			valid_queries++;
			transition_measure.valid_queries++;
			ratio_measure.valid_queries++;
			frequency_measure.valid_queries++;
			
			start = std::chrono::high_resolution_clock::now();
			stop = std::chrono::high_resolution_clock::now();
			auto query_time_transition_tmp = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
			auto query_time_ratio_tmp = query_time_transition_tmp;
			auto query_time_frequency_tmp = query_time_transition_tmp;

			auto start1 = std::chrono::high_resolution_clock::now();
			transition_predictor.GenerateRules(query_sequence, 1, repeated_events);
			transition_retrieved = transition_predictor.ReturnRulePrediction(top_k, repeated_events);
			auto stop1 = std::chrono::high_resolution_clock::now();
			query_time_transition_tmp += std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1);

			auto start2 = std::chrono::high_resolution_clock::now();
			ratio_predictor.GenerateRules(query_sequence, 2, repeated_events);
			ratio_retrieved = ratio_predictor.ReturnRulePrediction(top_k, repeated_events);
			auto stop2 = std::chrono::high_resolution_clock::now();
			query_time_ratio_tmp += std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - start2);

			auto start3 = std::chrono::high_resolution_clock::now();
			frequency_predictor.GenerateRules(query_sequence, 3, repeated_events);
			frequency_retrieved = frequency_predictor.ReturnRulePrediction(top_k, repeated_events);
			auto stop3 = std::chrono::high_resolution_clock::now();
			query_time_frequency_tmp += std::chrono::duration_cast<std::chrono::milliseconds>(stop3 - start3);

			printf("\r%.2lf%, [%d][%d]", 100.0*(double)t / query_database_sz, transition_measure.matched_queries, dump_query + (transition_measure.valid_queries - transition_measure.matched_queries));
			
		/*	if (transition_retrieved.size() < top_k || ratio_retrieved.size() < top_k || frequency_retrieved.size() < top_k) {
				mismatch_because_insufficient_K++;
				cout << "query sequence: " << query_sequence << "\n";
				cout << "ground truth sequence: " << ground_truth_sequence << "\n";
				cout << "transition: " << transition_retrieved << "\n";
				cout << "ratio: " << ratio_retrieved << "\n";
				cout << "frequency: " << frequency_retrieved << "\n";
				getchar();
				continue;
			}*/

			query_time_transition += query_time_transition_tmp;
			query_time_ratio += query_time_ratio_tmp;
			query_time_frequency += query_time_frequency_tmp;
			avg_query_length += static_cast<double>(query_sequence.size());
			success += 1.0;
		
			if (!transition_retrieved.empty()) {
				transition_measure.correct_count += Accuracy(transition_retrieved, ground_truth_sequence);
				transition_measure.precision += Precision(transition_retrieved, ground_truth_sequence);
				transition_measure.recall += Recall(transition_retrieved, ground_truth_sequence);
				transition_measure.r_precision += Rprecision(transition_retrieved, ground_truth_sequence);
				transition_measure.weighted_hit += WeightedHit(sequence_id, transition_retrieved, ground_truth_sequence);
				transition_measure.tf_idf += TFIDF(sequence_id, transition_retrieved, ground_truth_sequence);
				transition_measure.matched_queries++;
				transition_measure.predictions += transition_retrieved.size();
			}
			if (!ratio_retrieved.empty()) {
				ratio_measure.correct_count += Accuracy(ratio_retrieved, ground_truth_sequence);
				ratio_measure.precision += Precision(ratio_retrieved, ground_truth_sequence);
				ratio_measure.recall += Recall(ratio_retrieved, ground_truth_sequence);
				ratio_measure.r_precision += Rprecision(ratio_retrieved, ground_truth_sequence);
				ratio_measure.weighted_hit += WeightedHit(sequence_id, ratio_retrieved, ground_truth_sequence);
				ratio_measure.tf_idf += TFIDF(sequence_id, ratio_retrieved, ground_truth_sequence);
				ratio_measure.matched_queries++;
				ratio_measure.predictions += ratio_retrieved.size();
			}
			if (!frequency_retrieved.empty()) {
				frequency_measure.correct_count += Accuracy(frequency_retrieved, ground_truth_sequence);
				frequency_measure.precision += Precision(frequency_retrieved, ground_truth_sequence);
				frequency_measure.recall += Recall(frequency_retrieved, ground_truth_sequence);
				frequency_measure.r_precision += Rprecision(frequency_retrieved, ground_truth_sequence);
				frequency_measure.weighted_hit += WeightedHit(sequence_id, frequency_retrieved, ground_truth_sequence);
				frequency_measure.tf_idf += TFIDF(sequence_id, frequency_retrieved, ground_truth_sequence);
				frequency_measure.matched_queries++;
				frequency_measure.predictions += frequency_retrieved.size();
			}

		}

		if (transition_measure.matched_queries != 0) {
			//printf("\n");

			//transition_measure.PrintStat();
			//ratio_measure.PrintStat();
			//frequency_measure.PrintStat();

			//printf("transition average complete fit: %.2lf\n", static_cast<double>(transition_predictor.cnt_)/ static_cast<double>(transition_measure.matched_queries));
			//printf("ratio average complete fit: %.2lf\n", static_cast<double>(ratio_predictor.cnt_)/static_cast<double>(ratio_measure.matched_queries));
			//printf("frequency average complete fit: %.2lf\n", static_cast<double>(frequency_predictor.cnt_)/ static_cast<double>(frequency_measure.matched_queries));

			sprintf(buffer, "AVG Query Time: [%lfms], [%lfms], [%lfms]\n", static_cast<double>(query_time_transition.count()), static_cast<double>(query_time_ratio.count()), static_cast<double>(query_time_frequency.count()));
			cout << "AVG Query Time: [" << query_time_transition.count() << "ms] [" << query_time_ratio.count() << "ms] [" << query_time_frequency.count() << "ms]\n";
		}
		else
			printf("None\n");
	}
	transition_measure.Calculate();
	ratio_measure.Calculate();
	frequency_measure.Calculate();
	measures_[0][offset] = transition_measure;
	measures_[1][offset] = ratio_measure;
	measures_[2][offset] = frequency_measure;
	double fsz = static_cast<double>(fold_size);
	avg_query_set_sz /= fsz;
	transition_avg_pattern_length /= fsz;
	ratio_avg_pattern_length /= fsz;
	frequency_avg_pattern_length /= fsz;
	transition_avg_unique_item /= fsz;
	ratio_avg_unique_item /= fsz;
	frequency_avg_unique_item /= fsz;

	printf("split ratio: %.2lf Top patterns: %d\n", split_ratio, pattern_num_lim);
	printf("avg query length: %.2lf avg query set size: %.2lf\n", avg_query_length / success, avg_query_set_sz);
	printf("Transition: avg pattern length: %.2lf unique items: %.2lf\n", transition_avg_pattern_length, transition_avg_unique_item);
	printf("Ratio: avg pattern length: %.2lf unique items: %.2lf\n", ratio_avg_pattern_length, ratio_avg_unique_item);
	printf("Frequency: avg pattern length: %.2lf unique items: %.2lf\n", frequency_avg_pattern_length, frequency_avg_unique_item);
}

double Tester::Accuracy(const vector<int> &rec, const vector<int> &ans) {
	unordered_set<int> st;
	for (const auto& e : ans) st.insert(e);
	for (const auto&e : rec) {
		if (st.find(e) != st.end()) return 1.0;
	}
	return 0.0;
}

double Tester::Rprecision(const vector<int> &rec, const vector<int> &ans) {
	if (rec.size() <= ans.size()) return Precision(rec, ans);
	unordered_map<int, int> mp;
	for (const auto& e : ans) mp[e]++;
	int cnt = 0;
	for (int i = 0; i < ans.size(); ++i) {
		if (mp[rec[i]] > 0) {
			mp[rec[i]]--;
			cnt++;
		}
	}
	return (double)cnt / (double)ans.size();
}

double Tester::WeightedHit(const int sequence_id, const vector<int> &rec, const vector<int> & ans) {
	unordered_set<int> st;
	for (const auto& e : ans) st.insert(e);
	for (int i = 0; i< rec.size(); ++i) {
		const int &e = rec[i];
		if (st.find(e) != st.end()) {
			return (double)(ans.size() - i - 1) / (double)(ans.size());
		}
	}
	return 0.0;
}

//double Tester::NDCG(const int sequence_id, const vector<int> &rec, const vector<int> & ans) {
//	unordered_map<int, int> mp;
//	for (const auto& e : ans) mp[e]++;
//
//	double dcg = 0.0;
//	vector<double> vec;
//	for (int i = 0; i < rec.size(); ++i) {
//		const int &e = rec[i];
//		if (mp[e] > 0) {
//			mp[e]--;
//			dcg += static_cast<double>(database_.get_tfidf_score(e, sequence_id)) / log2(static_cast<double>(i) + 2.0);
//			vec.push_back(static_cast<double>(database_.get_tfidf_score(e, sequence_id)));
//		}
//	}
//	sort(vec.begin(), vec.end());
//	reverse(vec.begin(), vec.end());
//	double idcg = 0.0;
//	for (int i = 0; i < vec.size(); ++i) idcg += vec[i] / log2(static_cast<double>(i) + 2.0);
//	if (idcg == 0.0) return 0.0;
//	return dcg / idcg;
//}

double Tester::Precision(const vector<int> &rec, const vector<int> &ans) {
	unordered_map<int, int> mp;
	for (const auto& e : ans) mp[e]++;
	int cnt = 0;
	for (const auto&e : rec) {
		if (mp[e] > 0) {
			mp[e]--;
			cnt++;
		}
	}
	return (double)cnt / (double)rec.size();
}

double Tester::TFIDF(const int sequence_id, const vector<int> &rec, const vector<int> &ans) {
	unordered_map<int, int> mp;
	for (const auto& e : ans) mp[e]++;
	double val = 0.0;
	for (const auto&e : rec) {
		if (mp[e] > 0) {
			mp[e]--;
			val += database_.get_tfidf_score(e, sequence_id);
		}
	}
	return val / static_cast<double>(rec.size());
}


double Tester::Recall(const vector<int> &rec, const vector<int> &ans) {
	unordered_map<int, int> mp;
	for (const auto& e : ans) mp[e]++;
	int cnt = 0;
	for (const auto&e : rec) {
		if (mp[e] > 0) {
			mp[e]--;
			cnt++;
		}
	}
	return (double)cnt / (double)ans.size();
}

void Tester::WriteTopPatternInfo(int fold_size, string path, vector<int> pattern_numbers) {
	vector<double> transition[2];
	vector<double> ratio[2];
	vector<double> frequency[2];

	int sz = static_cast<int>(pattern_numbers.size());
	for (int i = 0; i < 2; ++i) {
		transition[i].resize(sz);
		ratio[i].resize(sz);
		frequency[i].resize(sz);
		for (int j = 0; j < sz; ++j) {
			transition[i][j] = 0.0;
			ratio[i][j] = 0.0;
			frequency[i][j] = 0.0;
		}
	}

	double total_pattern_num[3] = { 0.0, 0.0, 0.0 };


	for (int i = 0; i < fold_size; ++i) {
		string fileA = path + "foldA" + to_string(i + 1);
		string fileB = path + "foldB" + to_string(i + 1);
		string fileC = path + "foldC" + to_string(i + 1);

		auto transition_pattern_set = LoadFrequentPatterns(fileA, 1);
		auto ratio_pattern_set = LoadFrequentPatterns(fileB, 2);
		auto frequency_pattern_set = LoadFrequentPatterns(fileC, 3);


		PatternSetAnalysis TM(transition_pattern_set);
		PatternSetAnalysis RM(ratio_pattern_set);
		PatternSetAnalysis FM(frequency_pattern_set);

		total_pattern_num[0] += static_cast<int>(transition_pattern_set.size());
		total_pattern_num[1] += static_cast<int>(ratio_pattern_set.size());
		total_pattern_num[2] += static_cast<int>(frequency_pattern_set.size());

		for (int j = 0; j < sz; ++j) {
			int pattern_num_lim = pattern_numbers[j];
			TM.set_pattern_num_lim(pattern_num_lim);
			RM.set_pattern_num_lim(pattern_num_lim);
			FM.set_pattern_num_lim(pattern_num_lim);
			TM.Count();
			RM.Count();
			FM.Count();
			transition[0][j] += TM.avg_pattern_length;
			transition[1][j] += TM.number_unique_items;
			ratio[0][j] += RM.avg_pattern_length;
			ratio[1][j] += RM.number_unique_items;
			frequency[0][j] += FM.avg_pattern_length;
			frequency[1][j] += FM.number_unique_items;
		}
	}

	for(int i = 0; i <=2; ++i) total_pattern_num[i] /= static_cast<double>(fold_size);

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < sz; ++j) {
			transition[i][j] /= static_cast<double>(fold_size);
			ratio[i][j] /= static_cast<double>(fold_size);
			frequency[i][j] /= static_cast<double>(fold_size);
		}
	}

	char buffer[1000];
	ofstream outfile(path + "/Experiment/PatternAnalysis2.txt");
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < sz; ++j) {
			if (pattern_numbers[j] == -1) outfile << total_pattern_num[i] << "\t";
			else outfile << pattern_numbers[j] << "\t";
		}
	}
	outfile << "\n";
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < sz; ++j) {
			if (i == 0) outfile << transition[0][j] << "\t";
			else if(i == 1) outfile << ratio[0][j] << "\t";
			else if(i == 2) outfile << frequency[0][j] << "\t";
		}
	}
	outfile << "\n";
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < sz; ++j) {
			if (i == 0) outfile << transition[1][j] << "\t";
			else if (i == 1) outfile << ratio[1][j] << "\t";
			else if (i == 2) outfile << frequency[1][j] << "\t";
		}
	}
}

void Tester::ExperimentTopK(int fold_size, string path, int top_k, int repeated_events, vector<int> pattern_numbers, vector<double> split_ratio) {
	int split_sz = static_cast<int>(split_ratio.size());
	int top_patterns_sz = static_cast<int>(pattern_numbers.size());
	for (int i = 0; i < split_sz; ++i) {

		for (int j = 0; j < 3; ++j) {
			measures_[j].clear();
			measures_[j].resize(top_patterns_sz);
		}

		ofstream outfile(path + "/Experiment/SplitRatio" + to_string(split_ratio[i]) + ".txt");
		for (int j = 0; j < top_patterns_sz; ++j) {
			int number_of_active_patterns = pattern_numbers[j];

			// must check if test_loadable is renewed
			test_loadable(fold_size, split_ratio[i], top_k, number_of_active_patterns, path, repeated_events, 1, j);
		}

		// measures_ must be filled and foreach Measure struct, Calculate must be called.
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < top_patterns_sz; ++k) {
				outfile << pattern_numbers[k];
				if (!(j == 2 && k == top_patterns_sz - 1)) outfile << "\t";
			}
		}
		outfile << "\n";

		// Coverage
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < top_patterns_sz; ++k) {
				outfile << measures_[j][k].coverage;
				if (!(j == 2 && k == top_patterns_sz - 1)) outfile << "\t";
			}
		}
		outfile << "\n";

		// Accuracy
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < top_patterns_sz; ++k) {
				outfile << measures_[j][k].accuracy;
				if (!(j == 2 && k == top_patterns_sz - 1)) outfile << "\t";
			}
		}
		outfile << "\n";

		// Precision
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < top_patterns_sz; ++k) {
				outfile << measures_[j][k].precision;
				if (!(j == 2 && k == top_patterns_sz - 1)) outfile << "\t";
			}
		}
		outfile << "\n";

		// Recall
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < top_patterns_sz; ++k) {
				outfile << measures_[j][k].recall;
				if (!(j == 2 && k == top_patterns_sz - 1)) outfile << "\t";
			}
		}
		outfile << "\n";

		// R-Precision
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < top_patterns_sz; ++k) {
				outfile << measures_[j][k].r_precision;
				if (!(j == 2 && k == top_patterns_sz - 1)) outfile << "\t";
			}
		}
		outfile << "\n";

		// weightedF1
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < top_patterns_sz; ++k) {
				outfile << measures_[j][k].weightedF1;
				if (!(j == 2 && k == top_patterns_sz - 1)) outfile << "\t";
			}
		}
		outfile << "\n";

		// TF-IDF
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < top_patterns_sz; ++k) {
				outfile << measures_[j][k].tf_idf;
				if (!(j == 2 && k == top_patterns_sz - 1)) outfile << "\t";
			}
		}
		outfile << "\n";
	}

}

void Tester::ExecutionTimeTest(int fold_size, string path, string dataset_path, vector<int> pattern_numbers) {

	const int pn_sz = static_cast<int>(pattern_numbers.size());
	vector<vector<double>> threshold_matrix(3);
	for (int i = 0; i < 3; ++i) threshold_matrix[i].resize(pn_sz);

	cout << "Number of K : " << pn_sz << "\n";
	for (int i = 0; i < pn_sz; ++i) {
		vector<double> ret = KthThreshold(fold_size, path, pattern_numbers[i]);
		for (int j = 0; j < 3; ++j) {
			threshold_matrix[j][i] = ret[j];
		}
		cout << i + 1 << "/" << pn_sz << " finished.\n";
	}

	// Save Threshold Matrix
	ofstream outfile1(path + "/Experiment/ThresholdMatrix.txt");
	for (int i = 0; i < pn_sz; ++i) {
		outfile1 << pattern_numbers[i];
		if (i + 1 < pn_sz) outfile1 << "\t";
		else outfile1 << "\n";
	}
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < pn_sz; ++j) {
			outfile1 << threshold_matrix[i][j];
			if (j + 1 < pn_sz) outfile1 << "\t";
		}
		if (i + 1 < 3) outfile1 << "\n";
	}

	Database database;
	database.ReadDataSpmfFormat(dataset_path);
	vector<vector<double>> execution_time_matrix(3);
	for (int i = 0; i < 3; ++i) execution_time_matrix[i].resize(pn_sz);

	for (int i = 0; i < pn_sz; ++i) {
		vector<double> thresholds(3);
		for (int j = 0; j < 3; ++j) thresholds[j] = threshold_matrix[j][i];

		vector<double> ret = KthRunTime(path, pattern_numbers[i], thresholds, database);
		for (int j = 0; j < 3; ++j) {
			execution_time_matrix[j][i] = ret[j];
		}
		cout << i + 1 << "/" << pn_sz << "finished.\n";
	}

	// Save Execution Time Matrix
	ofstream outfile2(path + "/Experiment/ExecutionMatrix.txt");
	for (int i = 0; i < pn_sz; ++i) {
		outfile2 << pattern_numbers[i];
		if (i + 1 < pn_sz) outfile2 << "\t";
		else outfile2 << "\n";
	}
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < pn_sz; ++j) {
			outfile2 << execution_time_matrix[i][j];
			if (j + 1 < pn_sz) outfile2 << "\t";
		}
		if (i + 1 < 3) outfile2 << "\n";
	}

	cout << "finished.\n";
}

vector<double> Tester::KthThreshold(int fold_size, string path, int top_patterns) {

	// parameter, string = dataset, K = 1000 := number of patterns.
	
	vector<double> k_th_threshold(3);

	// fold 1 ~5, Kth threshold for each measure
	for (int i = 0; i < fold_size; ++i) {

		string fileA = path + "foldA" + to_string(i + 1);
		string fileB = path + "foldB" + to_string(i + 1);
		string fileC = path + "foldC" + to_string(i + 1);

		auto transition_pattern_set = LoadFrequentPatterns(fileA, 1);
		auto ratio_pattern_set = LoadFrequentPatterns(fileB, 2);
		auto frequency_pattern_set = LoadFrequentPatterns(fileC, 3);

		PatternSetAnalysis TM(transition_pattern_set);
		PatternSetAnalysis RM(ratio_pattern_set);
		PatternSetAnalysis FM(frequency_pattern_set);
		if (top_patterns != -1) {
			TM.set_pattern_num_lim(top_patterns);
			RM.set_pattern_num_lim(top_patterns);
			FM.set_pattern_num_lim(top_patterns);
		}

		k_th_threshold[0] += TM.GetLastThreshold();
		k_th_threshold[1] += RM.GetLastThreshold();
		k_th_threshold[2] += FM.GetLastThreshold();
	}

	for (int i = 0; i < 3; ++i) k_th_threshold[i] /= static_cast<double>(fold_size);
	return k_th_threshold;
}

vector<double> Tester::KthRunTime(string path, int top_patterns, vector<double> thresholds, Database &database) {
	printf("Execution Time Test for K: %d\n", top_patterns);

	PatternMiner recency_miner, compactness_miner;
	PrefixSpan presence_miner;

	// upperbound set to frequency support
	recency_miner.efficient_upperbound = 1;
	compactness_miner.efficient_upperbound = 1;

	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto mining_time_rec = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	auto mining_time_comp = mining_time_rec;
	auto mining_time_pre = mining_time_rec;

	recency_miner.LoadTrainingDatabase(database.get_full_db());
	compactness_miner.LoadTrainingDatabase(database.get_full_db());
	presence_miner.LoadTrainDatabase(database.get_full_db());

	start = std::chrono::high_resolution_clock::now();
	recency_miner.Run(thresholds[0], thresholds[0], 1);
	stop = std::chrono::high_resolution_clock::now();
	mining_time_rec = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

	start = std::chrono::high_resolution_clock::now();
	compactness_miner.Run(thresholds[1], thresholds[1], 2);
	stop = std::chrono::high_resolution_clock::now();
	mining_time_comp = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

	start = std::chrono::high_resolution_clock::now();
	presence_miner.Run(thresholds[2]);
	stop = std::chrono::high_resolution_clock::now();
	mining_time_pre = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

	vector<double> mining_time(3);

	mining_time[0] = static_cast<double>(mining_time_rec.count());
	mining_time[1] = static_cast<double>(mining_time_comp.count());
	mining_time[2] = static_cast<double>(mining_time_pre.count());

	return mining_time;
}


void Tester::ExperimentTimeNaive(int fold_size, string path, string dataset_path, vector<int> pattern_numbers) {
	const int pn_sz = static_cast<int>(pattern_numbers.size());
	vector<double> thresholds(pn_sz);

	cout << "Number of K : " << pn_sz << "\n";
	for (int i = 0; i < pn_sz; ++i) {
		vector<double> ret = KthThreshold(fold_size, path, pattern_numbers[i]);
		thresholds[i] = ret[0];
		cout << i + 1 << "/" << pn_sz << " finished.\n";
	}

	// Save Threshold Matrix
	ofstream outfile1(path + "/Experiment/Thresholds.txt");
	for (int i = 0; i < pn_sz; ++i) {
		outfile1 << pattern_numbers[i];
		if (i + 1 < pn_sz) outfile1 << "\t";
		else outfile1 << "\n";
	}
	{
		for (int j = 0; j < pn_sz; ++j) {
			outfile1 << thresholds[j];
			if (j + 1 < pn_sz) outfile1 << "\t";
		}

		Database database;
		database.ReadDataSpmfFormat(dataset_path);
		vector<double> execution_time(pn_sz);

		for (int i = 0; i < pn_sz; ++i) {
			double ret = NaiveCal(path, pattern_numbers[i], thresholds[i], database);
			execution_time[i] = ret;
			cout << i + 1 << "/" << pn_sz << "finished.\n";
		}

		// Save Execution Time Matrix
		ofstream outfile2(path + "/Experiment/ExecutionNaive.txt");
		for (int i = 0; i < pn_sz; ++i) {
			outfile2 << pattern_numbers[i];
			if (i + 1 < pn_sz) outfile2 << "\t";
			else outfile2 << "\n";
		}

		for (int j = 0; j < pn_sz; ++j) {
			outfile2 << execution_time[j];
			if (j + 1 < pn_sz) outfile2 << "\t";
		}

		cout << "finished.\n";
	}
}

double Tester::NaiveCal(string path, int top_patterns, double thresholds, Database &database) {
	printf("(Naive) Execution Time Test for K: %d\n", top_patterns);

	PrefixSpan presence_miner;
	vector<Sequence> db = database.get_full_db();
	presence_miner.LoadTrainDatabase(db);

	auto start = std::chrono::high_resolution_clock::now();
	presence_miner.naive = 1;
	presence_miner.Run(thresholds);
	auto patterns = presence_miner.GetFrequentPatterns();
	auto pattern_covers = presence_miner.GetPatternCovers();

	int patterns_sz = static_cast<int>(patterns.size());
	for (int k = 0; k < patterns_sz;++k) {
		auto pattern = patterns[k];
		double pattern_recency = 0.0;

		for (int sid : pattern_covers[k]) {
			vector<PatternInstance> pi_set;
			vector<PatternInstance> final_pi_set;


			const vector<int> &P = pattern.pattern;
			const int psz = P.size();

			vector<int> S = db[sid].sequence;
			const int id = sid;
			int sz = static_cast<int>(S.size());
			

			for (int i = 0; i < sz; ++i) {
				if (S[i] == P[0]) pi_set.push_back(PatternInstance(id, i, i));
			}

			for (int i = 1; i < psz; ++i) {
				vector<PatternInstance> ret;
				int pi_set_sz = static_cast<int>(pi_set.size());
				for (int i = 0; i < pi_set_sz; ++i) {
					PatternInstance &pi = pi_set[i];
					const int id = pi.sid;

					if (i + 1 < pi_set_sz) {
						sz = std::min(sz, pi_set[i + 1].r + 1);
					}

					bool flag = false;

					for (int j = pi.r + 1; j < sz; ++j) {
						if (S[j] == P[i]) {
							flag = true;
							ret.push_back(PatternInstance(id, pi.l, j));
							break;
						}
					}

					if (i == psz - 1 && !flag) {
						final_pi_set.push_back(pi);
					}
				}
				pi_set = ret;
			}

			if (final_pi_set.empty()) continue;

			double sum = 0.0;
			for (const auto &pi : final_pi_set) {
				sum += pi.value;
			}
			sum /= static_cast<double>(final_pi_set.size());
			pattern_recency += sum;
		}
		pattern_recency /= static_cast<double>(database.get_full_db().size());
	}


	auto stop = std::chrono::high_resolution_clock::now();
	auto mining_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

	return mining_time.count();
}

void Tester::CaseStudy(string path) {
	string fileA = path + "foldA1";
	string fileB = path + "foldB1";
	string fileC = path + "foldC1"; 

	vector<Pattern> pattern_set[3];
	cout << "Loading...";
	pattern_set[0] = LoadFrequentPatterns(fileA, 1);
	pattern_set[1] = LoadFrequentPatterns(fileB, 2);
	pattern_set[2] = LoadFrequentPatterns(fileC, 3);
	cout << "Run\n";

	vector<vector<int>> V = { {17, 46}, {17, 1}, {17, 48}, {155, 147}, {17, 46, 48}, {17, 45}, {46, 45, 15} };

	for (int i = 0; i < V.size(); ++i) {
		for (int k = 0; k < 3; ++k) {
			for (int j = 0; j < pattern_set[k].size(); ++j) {
				auto &P = pattern_set[k][j];
				if (P.pattern == V[i]) {
					cout << "[" << j + 1 << "] ";
					for (auto e : P.pattern) cout << e << " ";
					cout << "(";
					cout << P.interestingness << ") <";
					for (auto g : P.gap_sequence) cout << g << " ";
					cout << ">\n";
					break;
				}
			}
		}
	}


}