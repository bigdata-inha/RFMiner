#include "Database.h"

void Database::readNegOneDelimiterFormat(string file) {
	unordered_set<int> items;
	unordered_set<int> rep_check;
	ifstream infile(file);
	int len = 0;
	int x = 0;
	double sum = 0.0;
	int i = 0;
	vector<int> seq;
	int sequence_id = 0;
	while (!infile.eof()) {
		infile >> x;
		if (x == -1) continue;
		else if (x == -2) {
			sequence_database_.push_back(Sequence(sequence_id++, seq));
			sum += (double)seq.size();
			if (max_len < seq.size()) max_len = seq.size();
			seq.clear();
			rep_check.clear();
			continue;
		}
		seq.push_back(x);
		items.insert(seq.back());
		if (rep_check.find(x) != rep_check.end()) rep = 1;
		rep_check.insert(seq.back());
		if (maximum_index < x) maximum_index = x;
	}
	N = sequence_database_.size();
	unique_items = items.size();
	avg_len = sum / (double)N;
}

void Database::set_fold(unsigned int k) {
	k_fold = k;
	vector<int> shuf_vec;
	int fold_size = (N + k_fold - 1)/ k_fold;
	std::srand(unsigned(std::time(0)));
	for (int i = 0; i < sequence_database_.size(); ++i) shuf_vec.push_back(i);
	std::random_shuffle(shuf_vec.begin(), shuf_vec.end());

	int idx = 0;
	for (int i = 0; i < k_fold; ++i) {
		for (int j = 0; j < fold_size; ++j) {
			fold_db[i].push_back(sequence_database_[shuf_vec[idx]]);
			idx++;
			if (idx == N) break;
		}
		if (idx == N) break;
	}


	printf("Dividing %d folds\n", k);
	for (int i = 0; i < k_fold; ++i) {
		printf("%d ", fold_db[i].size());
	}
	printf("\n");
}

vector<Sequence> Database::get_full_db() {
	return sequence_database_;
}

vector<Sequence> Database::get_train_db(int fold_num) {
	train_db.clear();
	for (int i = 0; i < k_fold; ++i) {
		if (i == fold_num) continue;
		for (int j = 0; j < fold_db[i].size(); ++j) {
			train_db.push_back(fold_db[i][j]);
		}
	}
	return train_db;
}

vector<Sequence> Database::get_test_db(int fold_num) {
	test_db.clear();
	test_db = fold_db[fold_num];
	return test_db;
}

unsigned int Database::size() {
	return N;
}

void Database::CalculateTFIDF() {

	unordered_map<int, double> inverse_document_frequency;
	unordered_map<int, unordered_set<int>> document_frequency;

	int transaction_database_sz = static_cast<int>(sequence_database_.size());
	for (int i = 0; i < transaction_database_sz; ++i) {
		int sequence_id = sequence_database_[i].sequence_id;
		vector<int> &sequence = sequence_database_[i].sequence;
		for (const auto&item : sequence) {
			document_frequency[item].insert(sequence_id);
			frequency_per_transaction_[{item, sequence_id}] += 1.0;
		}
		unordered_set<int> vis;
		for (const auto&item : sequence) {
			if (vis.find(item) != vis.end()) continue;
			vis.insert(item);
			frequency_per_transaction_[{item, sequence_id}] /= static_cast<double>(sequence.size());
		}
	}

	for (auto it = document_frequency.begin(); it != document_frequency.end(); ++it) {
		int denominator = static_cast<int>(it->second.size()) + 1;
		if (denominator > transaction_database_sz) denominator = transaction_database_sz;
		inverse_document_frequency[it->first] = log2(static_cast<double>(transaction_database_sz)/(static_cast<double>(denominator)));
	}
	
	for (int i = 0; i < transaction_database_sz; ++i) {
		int sequence_id = sequence_database_[i].sequence_id;
		vector<int> &transaction = sequence_database_[i].sequence;
		for (const auto&item : transaction) {
			if (frequency_per_transaction_[{item, sequence_id}] * inverse_document_frequency[item] < 0.0) {
				printf("negative value in tfidf!");
				getchar();
				exit(1);
			}
			tfidf_[{item, sequence_id}] = frequency_per_transaction_[{item, sequence_id}] * inverse_document_frequency[item];
		}
	}
}

void Database::WriteTfIdfTable(const string &filename) {
	std::ofstream outfile(filename);
	int sz = static_cast<int>(sequence_database_.size());
	for (int i = 0; i < sz; ++i) {
		int sequence_id = sequence_database_[i].sequence_id;
		outfile << sequence_id << ":";
		for (const int &item : sequence_database_[i].sequence) {
			double val = tfidf_[{item, sequence_id}];
			outfile << " [" << item << ":" << val << "] ";
		}
		outfile << "\n";
	}
}


void Database::printDatabaseInformation() {
	cout << "loaded filename: " << filename << "\n";
	cout << "number of transcation: " << N << "\n";
	cout << "number of distinct items: " << unique_items << "\n";
	cout << "average sequence length: " << avg_len << "\n";
	cout << "maximum item index " << maximum_index << "\n";
	cout << "maximum sequence length " << max_len << "\n";
	cout << "repeated items in a single sequence " << rep << "\n";
}