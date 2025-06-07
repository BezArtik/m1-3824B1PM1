#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void ComputeTermsSin(float* ArrTerms, int count, float x);
void ComputeTermsCos(float* ArrTerms, int count, float x);
void ComputeTermsExp(float* ArrTerms, int count, float x);
void ComputeTermsLn(float* ArrTerms, int count, float x);

float DirectSummation(float* ArrTerms, int count);
float ReverseSummation(float* ArrTerms, int count);
float PairwiseSummation(float* ArrTerms, int start, int end);
float PairwiseSummationRes(float* ArrTerms, int count);

typedef float (*SummationFunction)(float* x, int count);

typedef void (*ComputeTermsFunction)(float* terms, int count, float x);

float ComputeTaylorRow(float x, int count, ComputeTermsFunction ComputeTerms, SummationFunction SummationMethod);
void CompareAndOutputResults(int functionType, int summationMethod, float x, int count);
