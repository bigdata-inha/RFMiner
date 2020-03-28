#include "Database.h"
#include "prefixSpan.h"
#include "Tester.h"
#include "Predictor.h"
#include "TestFunctions.h"
#include "PatternMiner.h"

void BibleTest() {
	Tester T;
	Database database;
	database.ReadDataSpmfFormat("Data/Bible/Bible.txt");
	database.print_stats();

	database.CalculateTFIDF();
	T.LoadDatabase(database);

	int fold_size = 5;
	double transition_ratio_init_threshold = 0.02;
	double transition_threshold = 0.02;
	double ratio_threshold = 0.02;
	double frequency_threshold = 0.02;
	int pattern_num_lim = 50;
	double split_ratio = 0.5;
	int topk = 3;
	string config = "[" + to_string(fold_size) + ", " + to_string(transition_threshold) + ", " + to_string(ratio_threshold) + ", " + to_string(frequency_threshold) + "]";

	string path = "Data/Bible/Test/";
	string result_filename = path + "result" + config + ".txt";

	T.set_result_filename(result_filename);

	//T.fold_divider(fold_size, transition_ratio_init_threshold, transition_threshold, ratio_threshold, frequency_threshold, path, database);

	vector<int> top_pattern_numbers = { 50, 100, 150, 200, 250 };
	//T.WriteTopPatternInfo(fold_size, path, top_pattern_numbers);

	vector<int> exe_top_pattern_numbers = { 1000, 2000, 3000, 4000, 5000 };
	//T.ExecutionTimeTest(fold_size, path, "Data/BMS1/BMSWebView(KDDCUP)59601.txt", exe_top_pattern_numbers);

	vector<int> naive_pattern_numbers = { 1000, 2000, 3000, 4000, 5000 };
	//T.ExperimentTimeNaive(fold_size, path, "Data/BMS1/BMSWebView(KDDCUP)59601.txt", naive_pattern_numbers);

	vector<double> split_ratio_vec = { 0.1, 0.3, 0.5, 0.7, 0.9 };
	T.ExperimentTopK(fold_size, path, topk, database.get_rep(), top_pattern_numbers, split_ratio_vec);
}

void BMSView1Test() {
	Tester T;
	Database database;
	database.ReadDataSpmfFormat("Data/BMS1/BMSWebView(KDDCUP)59601.txt");
	database.print_stats();

	database.CalculateTFIDF();
	T.LoadDatabase(database);

	int fold_size = 5;
	double transition_ratio_init_threshold = 0.0005;
	double transition_threshold = 0.0005;
	double ratio_threshold = 0.0005;
	double frequency_threshold = 0.0005;
	int pattern_num_lim = 200;
	double split_ratio = 0.5;
	int topk = 1;
	string config = "[" + to_string(fold_size) + ", " + to_string(transition_threshold) + ", " + to_string(ratio_threshold) + ", " + to_string(frequency_threshold) + "]";

	string path = "Data/BMS1/Test/";
	string result_filename = path + "result" + config + ".txt";

	T.set_result_filename(result_filename);

	//T.fold_divider(fold_size, transition_ratio_init_threshold, transition_threshold, ratio_threshold, frequency_threshold, path, database);
	
	vector<int> top_pattern_numbers = { 1000, 1500, 2000, 2500, 3000 };
	//T.WriteTopPatternInfo(fold_size, path, top_pattern_numbers);

	vector<int> exe_top_pattern_numbers = { 1000, 2000, 3000, 4000, 5000 };
	//T.ExecutionTimeTest(fold_size, path, "Data/BMS1/BMSWebView(KDDCUP)59601.txt", exe_top_pattern_numbers);

	vector<int> naive_pattern_numbers = { 1000, 2000, 3000, 4000, 5000 };
	//T.ExperimentTimeNaive(fold_size, path, "Data/BMS1/BMSWebView(KDDCUP)59601.txt", naive_pattern_numbers);
	
	vector<double> split_ratio_vec = { 0.1, 0.3, 0.5, 0.7, 0.9 };
	T.ExperimentTopK(fold_size, path, topk, database.get_rep(), top_pattern_numbers, split_ratio_vec);
}

void BMSView2Test() {
	Tester T;
	Database database;
	database.ReadDataSpmfFormat("Data/BMS2/BMSWebView2(KDDCUP)77512.txt");
	database.print_stats();

	database.CalculateTFIDF();
	T.LoadDatabase(database);

	int fold_size = 5;
	double transition_ratio_init_threshold = 0.0005;
	double transition_threshold = 0.0005;
	double ratio_threshold = 0.0005;
	double frequency_threshold = 0.0005;
	int pattern_num_lim = 2000;
	double split_ratio = 0.1;
	int topk = 1;
	string config = "[" + to_string(fold_size) + ", " + to_string(transition_threshold) + ", " + to_string(ratio_threshold) + ", " + to_string(frequency_threshold) + "]";

	string path = "Data/BMS2/Test/";
	string result_filename = path + "result" + config + ".txt";

	T.set_result_filename(result_filename);

	T.fold_divider(fold_size, transition_ratio_init_threshold, transition_threshold, ratio_threshold, frequency_threshold, path, database);
	
	vector<int> top_pattern_numbers = { 1000, 1500, 2000, 2500, 3000 };
	//T.WriteTopPatternInfo(fold_size, path, top_pattern_numbers);

	vector<int> exe_top_pattern_numbers = { 10000, 20000, 30000, 40000, 50000 };
	int efficient_upperbound = 0;
	//T.ExecutionTimeTest(fold_size, path, "Data/BMS2/BMSWebView2(KDDCUP)77512.txt", exe_top_pattern_numbers);

	vector<int> naive_pattern_numbers = { 1000, 2000, 3000, 4000, 5000 };
	//T.ExperimentTimeNaive(fold_size, path, "Data/BMS2/BMSWebView2(KDDCUP)77512.txt", naive_pattern_numbers);
	
	vector<double> split_ratio_vec = { 0.1, 0.3, 0.5, 0.7, 0.9 };
	//T.ExperimentTopK(fold_size, path, topk, database.get_rep(), top_pattern_numbers, split_ratio_vec);
}

void KosarakTest() {
	Tester T;
	Database database;
	database.ReadDataSpmfFormat("Data/Kosarak/Kosarak25000.txt");
	database.print_stats();

	database.CalculateTFIDF();
	T.LoadDatabase(database);

	int fold_size = 5;
	double transition_ratio_init_threshold = 0.002;
	double transition_threshold = 0.002;
	double ratio_threshold = 0.002;
	double frequency_threshold = 0.002;
	int pattern_num_lim = 20;
	double split_ratio = 0.5;
	int topk = 1;
	string config = "[" + to_string(fold_size) + ", " + to_string(transition_threshold) + ", " + to_string(ratio_threshold) + ", " + to_string(frequency_threshold) + "]";

	string path = "Data/Kosarak/Test/";
	string result_filename = path + "result" + config + ".txt";

	T.set_result_filename(result_filename);

	//T.fold_divider(fold_size, transition_ratio_init_threshold, transition_threshold, ratio_threshold, frequency_threshold, path, database);
	vector<int> top_pattern_numbers = { 100, 200, 300, 400, 500 };
	//T.WriteTopPatternInfo(fold_size, path, top_pattern_numbers);
	//vector<int> top_pattern_numbers = { 100, 200, 300, 400, 500, -1 };
	//T.WriteTopPatternInfo(fold_size, path, top_pattern_numbers);
	vector<double> split_ratio_vec = { 0.1, 0.3, 0.5, 0.7, 0.9 };
	T.ExperimentTopK(fold_size, path, topk, database.get_rep(), top_pattern_numbers, split_ratio_vec);
}

void MSNBCTest() {
	Tester T;
	Database database;
	database.ReadDataSpmfFormat("Data/MSNBC/MSNBC(31790).txt");
	database.print_stats();

	database.CalculateTFIDF();
	T.LoadDatabase(database);

	int fold_size = 5;
	double transition_ratio_init_threshold = 0.08;
	double transition_threshold = 0.08;
	double ratio_threshold = 0.08;
	double frequency_threshold = 0.08;
	int pattern_num_lim = 400;
	double split_ratio = 0.5;
	int topk = 1;
	string config = "[" + to_string(fold_size) + ", " + to_string(transition_threshold) + ", " + to_string(ratio_threshold) + ", " + to_string(frequency_threshold) + "]";

	string path = "Data/MSNBC/Test/";
	string result_filename = path + "result" + config + ".txt";

	T.set_result_filename(result_filename);

	T.fold_divider(fold_size, transition_ratio_init_threshold, transition_threshold, ratio_threshold, frequency_threshold, path, database);
	vector<int> top_pattern_numbers = { 100, 200, 300, 400};
	T.WriteTopPatternInfo(fold_size, path, top_pattern_numbers);
	vector<double> split_ratio_vec = { 0.1, 0.3, 0.5, 0.7, 0.9 };
	T.ExperimentTopK(fold_size, path, topk, database.get_rep(), top_pattern_numbers, split_ratio_vec);

}

void FIFATest() {
	Tester T;
	Database database;
	database.ReadDataSpmfFormat("Data/FIFA/FIFA.txt");
	database.print_stats();

	database.CalculateTFIDF();
	T.LoadDatabase(database);

	int fold_size = 5;
	double transition_ratio_init_threshold = 0.01;
	double transition_threshold = 0.01;
	double ratio_threshold = 0.01;
	double frequency_threshold = 0.01;
	int pattern_num_lim = 20;
	double split_ratio = 0.1;
	int topk = 1;
	string config = "[" + to_string(fold_size) + ", " + to_string(transition_threshold) + ", " + to_string(ratio_threshold) + ", " + to_string(frequency_threshold) + "]";

	string path = "Data/FIFA/Test/";
	string result_filename = path + "result" + config + ".txt";

	T.set_result_filename(result_filename);

	//T.fold_divider(fold_size, transition_ratio_init_threshold, transition_threshold, ratio_threshold, frequency_threshold, path, database);
	vector<int> top_pattern_numbers = { 1000, 1500, 2000, 2500, 3000 };
	//T.WriteTopPatternInfo(fold_size, path, top_pattern_numbers);
	
	vector<int> exe_top_pattern_numbers = { 5000, 8000, 11000, 14000, 17000 };
	int efficient_upperbound = 0;
	//T.ExecutionTimeTest(fold_size, path, "Data/FIFA/FIFA.txt", exe_top_pattern_numbers);

	vector<int> naive_pattern_numbers = { 5000, 8000, 11000, 14000, 17000 };
	//T.ExperimentTimeNaive(fold_size, path, "Data/FIFA/FIFA.txt", naive_pattern_numbers);

	vector<double> split_ratio_vec = { 0.1, 0.3, 0.5, 0.7, 0.9 };
	T.ExperimentTopK(fold_size, path, topk, database.get_rep(), top_pattern_numbers, split_ratio_vec);
}

void ExampleTest() {
	Database database;
	database.ReadDataSpmfFormat("b.txt");
	database.print_stats();

	double transition_ratio_init_threshold = 0.168;
	double transition_threshold = 0.168;
	double ratio_threshold = 0.168;
	double frequency_threshold = 0.5;
	
	PatternMiner transition, ratio, frequency;
	transition.LoadTrainingDatabase(database.get_full_db());
	ratio.LoadTrainingDatabase(database.get_full_db());
	frequency.LoadTrainingDatabase(database.get_full_db());

	transition.Run(transition_ratio_init_threshold, transition_ratio_init_threshold, 1);
	ratio.Run(transition_ratio_init_threshold, transition_ratio_init_threshold, 2);
	frequency.Run(frequency_threshold, frequency_threshold, 3);

	transition.WritePatternFile("TransitionResult.txt");
	ratio.WritePatternFile("RatioResult.txt");
	frequency.WritePatternFile("FrequencyResult.txt");
}

int main() {
	//BibleTest();
	//KosarakTest();
	BMSView2Test();
	//BMSView1Test();
	//MSNBCTest();
	//FIFATest();
	//FifaDataProcess("Data/FIFA/FIFA20450.txt");
	//DataGenerate();

	//ExperimentExecutionTime("Data/Kosarak/Kosarak25000.txt");
	//DebugMinerCorrectnessCheck();
	//TfIdfTableCheck();
	//ExampleTest();
	//TransitionWeightClassTest();
	//QueryTest();
	//ChoiceTest();
	
	//PrefixSpan ps;
	//PatternMiner rec;
	//const double threshold = 0.0023354;
	//Database db;
	//db.ReadDataSpmfFormat("Data/BMS2/BMSWebView2(KDDCUP)77512.txt");
	//db.print_stats();
	//
	//ps.LoadTrainDatabase(db.get_full_db());
	//rec.LoadTrainingDatabase(db.get_full_db());

	//auto start = std::chrono::high_resolution_clock::now();
	//auto stop = std::chrono::high_resolution_clock::now();
	//auto mining_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	//auto mining_time2 = mining_time;


	//start = std::chrono::high_resolution_clock::now();
	//ps.Run(threshold);
	//stop = std::chrono::high_resolution_clock::now();
	//mining_time2 += std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	//cout << "\n";
	//cout << mining_time2.count() << "\n";
	//cout << ps.cnt << "\n";

	//start = std::chrono::high_resolution_clock::now();
	//rec.Run(threshold, threshold, 1);
	//stop = std::chrono::high_resolution_clock::now();
	//mining_time += std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

	//cout << "\n";
	//cout << mining_time.count() << "\n";
	//cout << rec.cnt   << "\n";
	//Tester T;
	//T.CaseStudy("Data/FIFA/Test/");

	getchar();
}