#include "Database.h"
#include "prefixSpan.h"
#include "Tester.h"
#include "Predictor.h"
#include "TestFunctions.h"
#include "PatternMiner.h"

void BMSView1Test() {
	Tester testerAgent;
	Database database;
	string testResultOutputPath = "Data/BMS1/Test";

	database.readNegOneDelimiterFormat("Data/BMS1/BMSWebView(KDDCUP)59601.txt");
	database.printDatabaseInformation();

	testerAgent.LoadDatabase(database);

	// Fold Dividing and Prediction Test
	// Dividing Dabase into folds
	const int numberOfFolds = 5;
	const double foldDividingTestThreshold = 0.00045;
	const double recencyAndCompactnessInitialThreshold = 0.00045;
	const double recencyThreshold = 0.0005;
	const double compatnessThreshold = 0.00045;
	const double presenceThreshold = 0.0006;
	testerAgent.fold_divider(numberOfFolds, recencyAndCompactnessInitialThreshold, recencyThreshold, compatnessThreshold, presenceThreshold, testResultOutputPath, database);

	// Prediction Test
	int topk = 1;
	vector<double> split_ratio_vec = { 0.1, 0.3, 0.5, 0.7, 0.9 };
	vector<int> numberOfTopPatternsForPredictionTest = { 1000, 1500, 2000, 2500, 3000 };
	testerAgent.ExperimentTopK(numberOfFolds, testResultOutputPath, topk, database.isThereDuplicatedItermsInSingleEventSequence(), numberOfTopPatternsForPredictionTest, split_ratio_vec);
}

void BMSView2Test() {
	Tester testerAgent;
	Database database;
	string testResultOutputPath = "Data/BMS2/Test";

	database.readNegOneDelimiterFormat("Data/BMS2/BMSWebView2(KDDCUP)77512.txt");
	database.printDatabaseInformation();

	testerAgent.LoadDatabase(database);

	// Fold Dividing and Prediction Test
	// Dividing Dabase into folds
	const int numberOfFolds = 5;
	const double foldDividingTestThreshold = 0.0013;
	const double recencyAndCompactnessInitialThreshold = foldDividingTestThreshold;
	const double recencyThreshold = foldDividingTestThreshold;
	const double compatnessThreshold = foldDividingTestThreshold;
	const double presenceThreshold = foldDividingTestThreshold;
	testerAgent.fold_divider(numberOfFolds, recencyAndCompactnessInitialThreshold, recencyThreshold, compatnessThreshold, presenceThreshold, testResultOutputPath, database);

	// Prediction Test
	int topk = 1;
	vector<double> split_ratio_vec = { 0.1, 0.3, 0.5, 0.7, 0.9 };
	vector<int> numberOfTopPatternsForPredictionTest = { 1000, 1500, 2000, 2500, 3000 };
	testerAgent.ExperimentTopK(numberOfFolds, testResultOutputPath, topk, database.isThereDuplicatedItermsInSingleEventSequence(), numberOfTopPatternsForPredictionTest, split_ratio_vec);
}

void FIFATest() {
	Tester testerAgent;
	Database database;
	string testResultOutputPath = "Data/FIFA/Test";

	database.readNegOneDelimiterFormat("Data/FIFA/FIFA.txt");
	database.printDatabaseInformation();

	testerAgent.LoadDatabase(database);

	// Fold Dividing and Prediction Test
	// Dividing Dabase into folds
	const int numberOfFolds = 5;
	const double foldDividingTestThreshold = 0.008;
	const double recencyAndCompactnessInitialThreshold = foldDividingTestThreshold;
	const double recencyThreshold = foldDividingTestThreshold;
	const double compatnessThreshold = foldDividingTestThreshold;
	const double presenceThreshold = foldDividingTestThreshold;
	testerAgent.fold_divider(numberOfFolds, recencyAndCompactnessInitialThreshold, recencyThreshold, compatnessThreshold, presenceThreshold, testResultOutputPath, database);

	// Prediction Test
	int topk = 1;
	vector<double> split_ratio_vec = { 0.1, 0.3, 0.5, 0.7, 0.9 };
	vector<int> numberOfTopPatternsForPredictionTest = { 1000, 1500, 2000, 2500, 3000 };
	testerAgent.ExperimentTopK(numberOfFolds, testResultOutputPath, topk, database.isThereDuplicatedItermsInSingleEventSequence(), numberOfTopPatternsForPredictionTest, split_ratio_vec);
}

int main() {
	BMSView1Test();
	BMSView2Test();
	FIFATest();
	//KosarakTest();
	//BibleTest();
	//MSNBCTest();
	//SIGNTest();
	//retailTest();
	//foodmartTest();
	//chainstoreTest();
	getchar();
}

// archive

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
	BMSView1Test();
	BMSView2Test();
	FIFATest();
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