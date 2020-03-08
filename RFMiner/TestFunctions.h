#include "Database.h"
#include "prefixSpan.h"
#include "Predictor.h"
#include "Tester.h"
#include "PatternMiner.h"

void DebugMinerCorrectnessCheck();
void MinerCorrectnessCheck();
void FileIdenticalCheck(const string &source, const string &target);
void TfIdfTableCheck();
void DataGenerate();
void QueryTest();
void ChoiceTest();

void FifaDataProcess(const string &filename);