#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

void Merge(float* arr1, int size1, float* arr2, int size2, float* arr12);

void MergeSort(float* arr, int size, float* arr0);

void MergeComplexityOfWork(float* arr1, int size1, float* arr2, int size2, float* arr12);

void MergeSortComplexityOfWork(float* arr, int size, float* arr0);

void SelectionSort(float* arr, int size);

void SelectionSortComplexityOfWork(float* arr, int size);

void CombSort(float* arr, int size);

void CombSortComplexityOfWork(float* arr, int size);

void Fcount(float* arr, unsigned int counter[], int n, int offset);

void RadixSort(float* arr, int size);

void RadixSortWithSign(float* arr, int size);

void FcountComplexityOfWork(float* arr, unsigned int counter[], int size, int offset);

void RadixSortComplexityOfWork(float* arr, int size);

void RadixSortWithSignComplexityOfWork(float* arr, int size);

void TimeOfWork(void (*f)(float* arr, int size), int size);

void TimeOfWorkMerge(void (*f)(float* arr, int size, float* arr0), int size);

void ComplexityOfWork(void (*f)(float* arr, int size), int size);

void ComplexityOfWorkMerge(void (*f)(float* arr, int size, float* mas0), int size);

void ComplexityOfWorkRadix(void (*f)(float* arr, int size), void (*f1)(float* arr, int size), int size);

int cmp(const void* a, const void* b);

void ConfirmationOfCorrectness(void (*f)(float* arr, int size), int size);

void ConfirmationOfCorrectnessMerge(void (*f)(float* arr, int size, float* arr0), int size);
