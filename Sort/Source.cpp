#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

long long count = 0;
long long count1 = 0;
//���������� ��������
void Merge(float* arr1, int size1, float* arr2, int size2, float* arr12) {
	int j = 0, k = 0; // j - ������ ��� ������� mas1, k - ������ ��� ������� mas2
	for (int i = 0; i < size1 + size2; i++) { // �������� �� ���� ��������� mas12
		if (j == size1) { // ���� ��� �������� �� mas1 ��� ��������� � mas12
			arr12[i] = arr2[k++]; // ��������� ���������� �������� �� mas2
		}
		else {
			if (k == size2) { // ���� ��� �������� �� mas2 ��� ��������� � mas12
				arr12[i] = arr1[j++]; // ��������� ���������� �������� �� mas1
			}
			else { // ���� ��� �������� �������� � ����� ��������
				if (arr1[j] < arr2[k]) { // ���������� ������� �������� mas1 � mas2
					arr12[i] = arr1[j++]; // ���� ������� mas1 ������, ��������� ��� � mas12
				}
				else {
					arr12[i] = arr2[k++]; // ����� ��������� ������� mas2 � mas12
				}
			}
		}
	}
}
void MergeSort(float* arr, int size, float* arr0) {
	if (size <= 1) { // ���� ������ �������� 1 �������, �� ��� ������������
		return;
	}

	MergeSort(&arr[0], size / 2, &arr0[0]);   // ���������� ��������� ������ �������� �������
	MergeSort(&arr[size / 2], size - size / 2, &arr0[size / 2]);  // ���������� ��������� ������ �������� �������
	Merge(&arr[0], size / 2, &arr[size / 2], size - size / 2, &arr0[0]);   // ������� ��� ��������������� �������� � ������ mas0

	for (int i = 0; i < size; i++) {
		arr[i] = arr0[i];
	}
}

void MergeComplexityOfWork(float* arr1, int size1, float* arr2, int size2, float* arr12) {
	int j = 0, k = 0;
	for (int i = 0; i < size1 + size2; i++) {
		if (j == size1) {
			arr12[i] = arr2[k++];
			count++;
		}
		else {
			if (k == size2) {
				count++;
				arr12[i] = arr1[j++];
			}
			else {
				if (arr1[j] < arr2[k]) {
					count++;
					arr12[i] = arr1[j++];
				}
				else {
					count++;
					arr12[i] = arr2[k++];
				}
			}
		}
	}
}
void MergeSortComplexityOfWork(float* arr, int size, float* arr0) {
	if (size <= 1) {
		return;
	}

	MergeSortComplexityOfWork(&arr[0], size / 2, &arr0[0]);
	MergeSortComplexityOfWork(&arr[size / 2], size - size / 2, &arr0[size / 2]);
	MergeComplexityOfWork(&arr[0], size / 2, &arr[size / 2], size - size / 2, &arr0[0]);

	for (int i = 0; i < size; i++) {
		arr[i] = arr0[i];
		count++;
	}
}

//���������� �������
void SelectionSort(float* arr, int size) {
	int i_min; // ������ ������������ �������� � ������� ��������
	float min; // ����������� �������� � ������� ��������

	for (int i = 0; i < size - 1; i++) {   //�������� �� ���� ��������� �������, ����� ����������
		min = arr[i]; // ������������, ��� ������� ������� - �����������
		i_min = i; // ���������� ������ �������� ������������ ��������


		for (int j = i + 1; j < size; j++) {   //���� ����������� ������� � ���������� ����� �������
			if (min > arr[j]) { // ���� ������ ������� ������ �������� ��������
				i_min = j; // ��������� ������ ������������ ��������
				min = arr[j]; // ��������� �������� ��������
			}
		}

		// ����� ����������: ������ ��������� ����������� ������� �� ����� ��������
		float t = arr[i];
		arr[i] = arr[i_min];
		arr[i_min] = t;
	}
}

void SelectionSortComplexityOfWork(float* arr, int size) {
	int i_min;
	float min;
	for (int i = 0; i < size - 1; i++) {
		min = arr[i];
		count++;
		i_min = i;
		count++;

		for (int j = i + 1; j < size; j++) {
			count++;
			if (min > arr[j]) {
				count++;
				i_min = j;
				count++;
				min = arr[j];
				count++;
			}
		}
		float t = arr[i];
		count++;
		arr[i] = arr[i_min];
		count++;
		arr[i_min] = t;
		count++;
	}
}

//���������� ���������
void CombSort(float* arr, int size) {
	int gap = size;  // ��������� ���
	float ratio = 1.27f;  // ����������� ���������� ����
	int flag = 1;  // ���� ��� ������������ �������

	while (gap > 1 || flag) {
		gap /= ratio;  // ��������� ���
		if (gap < 1) {
			gap = 1;
		}
		flag = 0;
		for (int i = 0; i < size - gap; i++) {  // ���������� �������� � ������ ��
			if (arr[i] > arr[i + gap]) {
				float t = arr[i];
				arr[i] = arr[i + gap];
				arr[i + gap] = t;
				flag = 1;  // ����� ���������
			}
		}
	}
}
void CombSortComplexityOfWork(float* arr, int size) {
	int gap = size;
	float ratio = 1.27f;
	int flag = 1;
	while (gap > 1 || flag) {
		gap /= ratio;
		if (gap < 1) {
			gap = 1;
		}
		flag = 0;
		for (int i = 0; i < size - gap; i++) {
			if (arr[i] > arr[i + gap]) {
				count++;
				float t = arr[i];
				arr[i] = arr[i + gap];
				arr[i + gap] = t;
				flag = 1;
			}
		}
	}
}

//����������� ����������
void Fcount(float* arr, unsigned int counter[], int size, int offset) {
	unsigned char* b = (unsigned char*)arr; // �������� ��������� �� ������ float � unsigned char
	memset(counter, 0, 256 * sizeof(unsigned int)); // �������� ������� ��� 256 ��������� �������� ������

	for (int i = 0; i < size; i++) {  // ������� ���������� ��������� ������� �����
		counter[b[i * sizeof(float) + offset]]++; // ����������� ������� ��� �������� �����
	}

	for (int i = 1; i < 256; i++) {  // ����������� ������� � ������� ��� ����������
		counter[i] += counter[i - 1]; // ��������� �������� ��� ��������� ��������
	}
}

void RadixSort(float* arr, int size) {
	unsigned int counter[256]; // ������ ��� �������� ��������� ������
	float* arr1 = (float*)malloc(size * sizeof(float)); // �������� ������ ��� ���������� �������
	if (arr1 == NULL) {
		return;
	}

	unsigned char* b = (unsigned char*)arr; // �������� ��������� �� ������ float � unsigned char

	for (int i = 0; i < sizeof(float); i++) {  // �������� ���� �� ������� ����� float (4 �����)
		Fcount(arr, counter, size, i); // �������� ������� Fcount ��� �������� ��������� ������

		for (int j = size - 1; j >= 0; j--) {   // ����������� ��������� � ��������������� ������ mas1
			arr1[counter[b[j * sizeof(float) + i]] - 1] = arr[j]; // ������������ �������� � ������ mas1 �� �������� �� counter
			counter[b[j * sizeof(float) + i]]--; // ��������� �������
		}
		// �������� ��������������� �������� ������� � ������������ ������ mas
		memcpy(arr, arr1, size * sizeof(float)); // �������� �������� �� mas1 ������� � mas

	}
	free(arr1); // ����������� ������, ���������� ��� ������ mas1
}
void RadixSortWithSign(float* arr, int size) {
	float* pos = (float*)malloc(size * sizeof(float)); // �������� ������ ��� ������� ������������� � ������� ������������� �����
	float* neg = (float*)malloc(size * sizeof(float));
	int posCount = 0, negCount = 0;
	if (pos == NULL || neg == NULL) {
		return;
	}
	// ��������� ������������� � ������������� �����
	for (int i = 0; i < size; i++) {
		if (arr[i] >= 0) {
			pos[posCount++] = arr[i];
		}
		else {
			neg[negCount++] = arr[i];
		}
	}

	// ��������� ������������� � ������������� ����� ��������
	RadixSort(neg, negCount);
	RadixSort(pos, posCount);

	// ���������� ��������������� �������
	for (int i = 0; i < negCount; i++) {
		arr[i] = neg[negCount - i - 1]; // �������� ������� ��� ������������� �����
	}

	for (int i = 0; i < posCount; i++) {
		arr[negCount + i] = pos[i]; // ������������� ����� ���� ����� �������������
	}

	free(pos);
	free(neg);
}

void FcountComplexityOfWork(float* arr, unsigned int counter[], int n, int offset) {
	unsigned char* b = (unsigned char*)arr;
	memset(counter, 0, 256 * sizeof(unsigned int));
	for (int i = 0; i < n; i++) {
		counter[b[i * sizeof(float) + offset]]++;
	}

	for (int i = 1; i < 256; i++) {
		counter[i] += counter[i - 1];
	}
}

void RadixSortComplexityOfWork(float* arr, int size) {
	unsigned int counter[256];
	float* mas1 = (float*)malloc(size * sizeof(float));
	if (mas1 == NULL) {
		return;
	}

	unsigned char* c = (unsigned char*)arr;

	for (int i = 0; i < sizeof(float); i++) {
		FcountComplexityOfWork(arr, counter, size, i);

		for (int j = size - 1; j >= 0; j--) {
			mas1[counter[c[j * sizeof(float) + i]] - 1] = arr[j];
			count++;
			counter[c[j * sizeof(float) + i]]--;
		}

		memcpy(arr, mas1, size * sizeof(float));
	}
	free(mas1);
}
void RadixSortWithSignComplexityOfWork(float* arr, int size) {
	float* pos = (float*)malloc(size * sizeof(float));
	float* neg = (float*)malloc(size * sizeof(float));
	int posCount = 0, negCount = 0;
	if (pos == NULL || neg == NULL) {
		return;
	}
	for (int i = 0; i < size; i++) {
		if (arr[i] >= 0) {
			pos[posCount++] = arr[i];
			count1++;
		}
		else {
			neg[negCount++] = arr[i];
			count1++;
		}
	}

	RadixSort(neg, negCount);
	RadixSort(pos, posCount);

	for (int i = 0; i < negCount; i++) {
		arr[i] = neg[negCount - i - 1];
		count1++;
	}

	for (int i = 0; i < posCount; i++) {
		arr[negCount + i] = pos[i];
		count1++;
	}

	free(pos);
	free(neg);
}


//������� ������� ������ ���������
void TimeOfWork(void (*f)(float* arr, int size), int size) {
	float* arr = (float*)malloc(sizeof(float) * size);
	if (arr == NULL) {
		return;
	}
	for (int i = 0; i < size; i++) {
		arr[i] = (float)rand() / (float)(RAND_MAX / (float)2000000) - 1000000;
	}
	clock_t start = clock();
	f(arr, size);
	clock_t finish = clock();
	double time_taken = ((double)(finish - start)) / CLOCKS_PER_SEC;
	printf("%f seconds\n", time_taken);

	free(arr);
}

void TimeOfWorkMerge(void (*f)(float* arr, int size, float* arr0), int size) {
	float* arr = (float*)malloc(sizeof(float) * size);
	float* arr0 = (float*)malloc(sizeof(float) * size);
	if (arr == NULL || arr0 == NULL) {
		return;
	}
	for (int i = 0; i < size; i++) {
		arr[i] = (float)rand() / (float)(RAND_MAX / (float)2000000) - 1000000;
		arr0[i] = 0;
	}
	clock_t start = clock();
	f(arr, size, arr0);
	clock_t finish = clock();
	double time_taken = ((double)(finish - start)) / CLOCKS_PER_SEC;
	printf("%f seconds\n", time_taken);

	free(arr);
	free(arr0);
}

//������� ��������� ������ ���������
void ComplexityOfWork(void (*f)(float* arr, int size), int size) {
	float* arr = (float*)malloc(sizeof(float) * size);
	if (arr == NULL) {
		return;
	}
	for (int i = 0; i < size; i++) {
		arr[i] = (float)rand() / (float)(RAND_MAX / (float)2000000) - 1000000;
	}

	f(arr, size);

	printf("%lli operation\n", count);
	free(arr);
}

void ComplexityOfWorkMerge(void (*f)(float* arr, int size, float* mas0), int size) {
	float* arr = (float*)malloc(sizeof(float) * size);
	float* arr0 = (float*)malloc(sizeof(float) * size);
	if (arr == NULL || arr0 == NULL) {
		return;
	}
	for (int i = 0; i < size; i++) {
		arr[i] = (float)rand() / (float)(RAND_MAX / (float)2000000) - 1000000;
		arr0[i] = 0;
	}

	f(arr, size, arr0);

	printf("%lli operation\n", count);
	free(arr);
	free(arr0);
}

void ComplexityOfWorkRadix(void (*f)(float* arr, int size), void (*f1)(float* arr, int size), int size) {
	float* arr = (float*)malloc(sizeof(float) * size);
	if (arr == NULL) {
		return;
	}
	for (int i = 0; i < size; i++) {
		arr[i] = (float)rand() / (float)(RAND_MAX / (float)2000000) - 1000000;
	}

	f(arr, size);
	f1(arr, size);

	printf("%lli operation\n", count + count1);
	free(arr);
}

int cmp(const void* a, const void* b) {
	return (*(float*)a > *(float*)b) - (*(float*)a < *(float*)b);
}

//������������� ������������ ����������
void ConfirmationOfCorrectness(void (*f)(float* arr, int size), int size) {
	float* arr = (float*)malloc(sizeof(float) * size);
	float* arr1 = (float*)malloc(sizeof(float) * size);
	int flag = 1;
	if (arr == NULL || arr1 == NULL) {
		return;
	}
	for (int i = 0; i < size; i++) {
		arr[i] = (float)rand() / (float)(RAND_MAX / (float)2000000) - 1000000;
		arr1[i] = arr[i];
	}
	f(arr, size);
	qsort(arr1, size, sizeof(float), cmp);
	for (int i = 0; i < size; i++) {
		if (arr[i] != arr1[i]) {
			flag = 0;
			break;
		}
	}
	if (flag == 1) {
		printf("The array has been sorted successfully!\n");
	}
	else {
		printf("The array wasn't sorted.\n");
	}
	free(arr);
	free(arr1);
}
void ConfirmationOfCorrectnessMerge(void (*f)(float* arr, int size, float* arr0), int size) {
	float* arr = (float*)malloc(sizeof(float) * size);
	float* arr0 = (float*)malloc(sizeof(float) * size);
	float* arr1 = (float*)malloc(sizeof(float) * size);
	int flag = 1;
	if (arr == NULL || arr1 == NULL || arr0 == NULL) {
		return;
	}
	for (int i = 0; i < size; i++) {
		arr[i] = (float)rand() / (float)(RAND_MAX / (float)2000000) - 1000000;
		arr1[i] = arr[i];
		arr0[i] = 0;
	}
	f(arr, size, arr0);
	qsort(arr1, size, sizeof(float), cmp);
	for (int i = 0; i < size; i++) {
		if (arr[i] != arr1[i]) {
			flag = 0;
			break;
		}
	}
	if (flag == 1) {
		printf("The array has been sorted successfully!\n");
	}
	else {
		printf("The array wasn't sorted.\n");
	}
	free(arr);
	free(arr0);
	free(arr1);
}