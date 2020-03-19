#pragma once
#include "Database.h"
#include "prefixSpan.h"
#include "Predictor.h"
#include "Tester.h"
#include "TestFunctions.h"
#include "PatternMiner.h"

void DebugCorrectnessCheck();
void CorrectnessCheck();
void FileIdenticalCheck(const string &source, const string &target);
void TfIdfTableCheck();
void DataGenerate();
void QueryTest();
void ChoiceTest();

void FifaDataProcess(const string &filename);