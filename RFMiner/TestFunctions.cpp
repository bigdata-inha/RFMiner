#include "TestFunctions.h"

#define TEST_FILE "Test/test.txt"
#define TEST_RECENCY_FILE "Test/RecencyResult.txt"
#define TEST_COMPACTNESS_FILE "Test/CompactnessResult.txt"
#define TEST_PRESENCE_FILE "Test/PresenceResult.txt"
#define ANSWER_RECENCY_FILE "Test/RecencyAnswer.txt"
#define ANSWER_COMPACTNESS_FILE "Test/CompactnessAnswer.txt"
#define ANSWER_PRESENCE_FILE "Test/PresenceAnswer.txt"
#define THRESHOLD 0.13

void DebugMinerCorrectnessCheck() {
	MinerCorrectnessCheck();
	FileIdenticalCheck(TEST_RECENCY_FILE, ANSWER_RECENCY_FILE);
	FileIdenticalCheck(TEST_COMPACTNESS_FILE, ANSWER_COMPACTNESS_FILE);
	FileIdenticalCheck(TEST_PRESENCE_FILE, ANSWER_PRESENCE_FILE);
	cout << "Correct!\n";
}

void MinerCorrectnessCheck() {
	Database database;
	database.ReadDataSpmfFormat(TEST_FILE);
	database.print_stats();
	//database.printDatabase();

	PatternMiner recencyMiner;
	PatternMiner compactnessMiner;
	PrefixSpan presenceMiner;

	// Loading databases to pattern miners
	recencyMiner.LoadTrainingDatabase(database.get_full_db());
	compactnessMiner.LoadTrainingDatabase(database.get_full_db());
	presenceMiner.LoadTrainDatabase(database.get_full_db());

	// Run Pattern Miners
	recencyMiner.Run(THRESHOLD, THRESHOLD, RECENCY);
	compactnessMiner.Run(THRESHOLD, THRESHOLD, COMPACTNESS);
	presenceMiner.Run(THRESHOLD);

	// Wrtie test result files
	recencyMiner.WritePatternFile(TEST_RECENCY_FILE);
	compactnessMiner.WritePatternFile(TEST_COMPACTNESS_FILE);
	presenceMiner.WriteFile(TEST_PRESENCE_FILE);
}

void FileIdenticalCheck(const string &source, const string &target) {
	ifstream source_stream(source);
	ifstream target_stream(target);
	string source_string, target_string;
	while (std::getline(source_stream, source_string) && std::getline(target_stream, target_string)) {
		if (source_string == target_string) continue;
		else {
			cout << source << " is not identical with " << target << "\n";
			getchar();
			exit(-1);
		}
	}
}

void DataGenerate() {
	ofstream outfile("b.txt");
	int x;
	int noise = 100;
	while (1) {
		cout << "1. insert new sequence, 2. finish\n";
		std::cin >> x;
		if (x == 1) {
			vector<int> vec;
			while (1) {
				cout << "insert element, -1. insert noise, -2.finish\n";
				std::cin >> x;
				if (x == -2) {
					vec.push_back(-2);
					break;
				}
				else if (x == -1) {
					cout << "how many times?\n";
					std::cin >> x;
					for (int i = 0; i < x; ++i) {
						vec.push_back(noise++);
						vec.push_back(-1);
					}
				}
				else {
					vec.push_back(x);
					vec.push_back(-1);
				}
			}
			cout << "inserting sequence is: ";
			for (int i = 0; i < vec.size(); ++i) cout << vec[i] << " ";
			cout << "\ninsert how many times?\n";
			std::cin >> x;
			for (int j = 0; j < x; ++j) {
				for (int i = 0; i < vec.size(); ++i) {
					if (vec[i] >= 131) vec[i] = vec[i] + j * 5;
					outfile << vec[i] << " ";
				}
				outfile << "\n";
			}
		}
		else break;
	}
}

void TfIdfTableCheck() {
	Database database;
	database.ReadDataSpmfFormat("DebugTestSet/Example.txt");
	database.print_stats();
	database.CalculateTFIDF();
	database.WriteTfIdfTable("DebugTestSet/TFIDF.txt");
}

void QueryTest() {
	vector<int> query = { 1, 2, 3, 1, 1, 2, 2, 3, 3, 4, 1, 2, 3, 2, 3, 4 };
	vector<int> pattern = { 1, 2, 3, 4, 5};
	vector<double> gaps = {1.0, 1.0, 1.0, 1.0};
	
	vector<Pattern> sequential_patterns;
	sequential_patterns.push_back(Pattern(1, pattern, 1.0, 1.0, gaps));

	Predictor predictor(sequential_patterns);
	predictor.SetDebug();

	predictor.GenerateRules(query, 1, 0);
	
	predictor.PrintRuleList(3);
}

void ChoiceTest() {
	vector<int> query = { 1, 10, 2, 5, 10, 6, 10 };
	vector<int> pattern1 = { 1, 2, 3 };
	vector<int> pattern2 = { 5, 6, 7 };
	vector<double> gaps = { 2.0, 5.0 };

	vector<Pattern> sequential_patterns;
	sequential_patterns.push_back(Pattern(1, pattern1, 1.0, 1.0, gaps));
	sequential_patterns.push_back(Pattern(1, pattern2, 1.0, 1.0, gaps));

	Predictor predictor(sequential_patterns);
	predictor.SetDebug();

	predictor.GenerateRules(query, 1, 0);

	predictor.PrintRuleList(3);
}

void FifaDataProcess(const string &filename) {
	ofstream outfile("Data/FIFA/FIFA.txt");
	ifstream infile(filename);
	int x;
	vector<int> seq;
	while (!infile.eof()) {
		infile >> x;
		if (x == -1) continue;
		else if (x == -2) {
			int sz = (int)ceil((double)seq.size()*(1.0 / 4));
			for (int i = 0; i < sz; ++i) {
				outfile << seq[i] << " -1 ";
			}
			outfile << "-2\n";
			seq.clear();
			continue;
		}
		seq.push_back(x);
	}
}