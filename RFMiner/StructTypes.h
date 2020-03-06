#ifndef STRUCTTYPES_H_
#define STRUCTTYPES_H_

#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <map>
#include <algorithm>

using std::vector;
using std::pair;
using std::string;
using std::set;
using std::unordered_set;
using std::unordered_map;
using std::cout;
using std::ofstream;
using std::ifstream;
using std::map;

struct Sequence {
	int sequence_id;
	vector<int> sequence;

	Sequence(int a_sequence_id, vector<int> a_sequence)
		: sequence_id(a_sequence_id), sequence(a_sequence) {}

	int &operator[] (int x) {
		return sequence[x];
	}

	void push_back(int x) {
		sequence.push_back(x);
	}

	int size() {
		return static_cast<int>(sequence.size());
	}
};

struct Pattern {
	int pattern_id;
	vector<int> pattern;
	vector<double> gap_sequence;
	double interestingness;
	double frequency;
	int type;

	Pattern(int a_type, vector<int> a_pattern, double a_interestingness, double a_frequency)
		: type(a_type), pattern(a_pattern), interestingness(a_interestingness), frequency(a_frequency) {
	}

	Pattern(int a_type, vector<int> a_pattern, double a_interestingness, double a_frequency, vector<double> a_gap_sequence)
		: type(a_type), pattern(a_pattern), interestingness(a_interestingness), frequency(a_frequency), gap_sequence(a_gap_sequence) {
	}

	bool operator <(const Pattern &rhs) const {
		if (type == 1 || type ==2 || type == 3) {
			return interestingness > rhs.interestingness;
		}
		else {
			printf("Error in pattern type.");
			getchar();
			exit(-1);
		}
	}

	void PrintPattern() {
		int psz = static_cast<int>(pattern.size());
		int wsz = static_cast<int>(gap_sequence.size());
		printf("[");
		for (int i = 0; i < psz; ++i) printf("%d ", pattern[i]);
		if (type == 3) {
			printf("] [%lf]\n", frequency);
		}
		else if (type == 2) {
			printf("] [%lf]\n", interestingness);
		}
		else if (type == 1) {
			printf("] [%lf]\n", interestingness);
			for (int i = 0; i < wsz; ++i) printf("%lf ", gap_sequence[i]);
			printf("\n");
		}
		else {
			printf("Error in pattern type.");
			getchar();
			exit(-1);
		}
	}
};


struct PatternInstance {
	int sid, l, r, ext_len, type;
	vector<int> landmark;
	double value, denom;

	PatternInstance()
		:sid(-1), l(100000), r(100000) {}

	PatternInstance(int a_sid, int a_l, int a_r)
		:sid(a_sid), l(a_l), r(a_r) {};

	PatternInstance(int a_type, int a_sid, int a_l, int a_r, int a_ext_len, vector<int> a_landmark)
		:type(a_type), sid(a_sid), l(a_l), r(a_r), ext_len(a_ext_len), landmark(a_landmark){
		int sz = static_cast<int>(landmark.size());
		double nominator = 0.0;
		double denominator = 0.0;

		if (type == 1) {
			for (int i = 0; i < sz; ++i) {
				nominator += static_cast<double>((i + 1));
				denominator += static_cast<double>((landmark.back() - landmark[i] + 1));
			}
		}
		else if (type == 2) {
			nominator = static_cast<double>(landmark.size());
			denominator = static_cast<double>((landmark.back() - landmark.front() + 1));
		}
		denom = denominator;
		value = nominator / denominator;
	};

	bool operator <(const PatternInstance &rhs) const {
		if (l == rhs.l) return sid < rhs.sid;
		return l < rhs.l;
	}

	vector<double> GapSequence() const{
		vector<double> gap_sequence;
		int lsz = static_cast<int>(landmark.size());
		for (int i = 0; i < lsz - 1; ++i) {
			gap_sequence.push_back(static_cast<double>(landmark[i + 1] - landmark[i]));
		}
		return gap_sequence;
	}
};

#endif // !STRUCTTYPES_H_