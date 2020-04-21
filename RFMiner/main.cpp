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
	database.readNegOneDelimiterFormat("b.txt");
	database.printDatabaseInformation();

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