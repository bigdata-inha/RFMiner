#ifndef DATABASE_H_
#define DATABASE_H_
#include "StructTypes.h"

class Database {
public:
	Database() {};
	~Database() {};

	void readNegOneDelimiterFormat(string file);

	void set_fold(unsigned int k);
	vector<Sequence> get_full_db();
	vector<Sequence> get_train_db(int fold_num);
	vector<Sequence> get_test_db(int fold_num);

	unsigned int size();
	void printDatabaseInformation();

	unsigned int get_max_seq_len() { return max_len; }
	int isThereDuplicatedItermsInSingleEventSequence() { return rep; }

	void CalculateTFIDF();

	double get_tfidf_score(int item, int sequence_id) {
		return tfidf_[{item, sequence_id}];
	}

	bool isEmpty() { return tfidf_.empty(); }
	void WriteTfIdfTable(const string &filename);

	void printDatabase() {
		int num = 1;
		for (Sequence &s : sequence_database_) {
			printf("%d: ", num++);
			for (int i = 0; i < s.size(); ++i) {
				if(s[i] >= 100) printf("x ");
				else printf("%c ", 'A' + s[i] - 1);
			}
			printf("\n");
		}
	}

private:
	struct HASH {
		size_t operator()(const pair<int, int>&x)const {
			return std::hash<long long>()(((long long)x.first) ^ (((long long)x.second) << 32));
		}
	};

	string filename;
	unsigned int N;
	unsigned int unique_items;
	unsigned int k_fold;
	unsigned fold_size;
	unsigned int max_len;
	unsigned int maximum_index;
	double avg_len;
	int rep;

	vector<Sequence> sequence_database_;
	vector<Sequence> fold_db[15];
	vector<Sequence> train_db;
	vector<Sequence> test_db;

	unordered_map<pair<int, int>, double, HASH> frequency_per_transaction_; // length adjusted
	unordered_map<pair<int, int>, double, HASH> tfidf_;

};

#endif // !DATABASE_H_