#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

long long count = 0;
long long count1 = 0;
//СОРТИРОВКА СЛИЯНИЕМ
void Merge(float* arr1, int size1, float* arr2, int size2, float* arr12) {
	int j = 0, k = 0; // j - индекс для массива mas1, k - индекс для массива mas2
	for (int i = 0; i < size1 + size2; i++) { // Проходим по всем элементам mas12
		if (j == size1) { // Если все элементы из mas1 уже добавлены в mas12
			arr12[i] = arr2[k++]; // Добавляем оставшиеся элементы из mas2
		}
		else {
			if (k == size2) { // Если все элементы из mas2 уже добавлены в mas12
				arr12[i] = arr1[j++]; // Добавляем оставшиеся элементы из mas1
			}
			else { // Если еще остались элементы в обоих массивах
				if (arr1[j] < arr2[k]) { // Сравниваем текущие элементы mas1 и mas2
					arr12[i] = arr1[j++]; // Если элемент mas1 меньше, добавляем его в mas12
				}
				else {
					arr12[i] = arr2[k++]; // Иначе добавляем элемент mas2 в mas12
				}
			}
		}
	}
}
void MergeSort(float* arr, int size, float* arr0) {
	if (size <= 1) { // Если массив содержит 1 элемент, он уже отсортирован
		return;
	}

	MergeSort(&arr[0], size / 2, &arr0[0]);   // Рекурсивно сортируем первую половину массива
	MergeSort(&arr[size / 2], size - size / 2, &arr0[size / 2]);  // Рекурсивно сортируем вторую половину массива
	Merge(&arr[0], size / 2, &arr[size / 2], size - size / 2, &arr0[0]);   // Сливаем две отсортированные половины в массив mas0

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

//СОРТИРОВКА ВЫБОРОМ
void SelectionSort(float* arr, int size) {
	int i_min; // Индекс минимального элемента в текущей итерации
	float min; // Минимальное значение в текущей итерации

	for (int i = 0; i < size - 1; i++) {   //проходим по всем элементам массива, кроме последнего
		min = arr[i]; // Предполагаем, что текущий элемент - минимальный
		i_min = i; // Запоминаем индекс текущего минимального элемента


		for (int j = i + 1; j < size; j++) {   //ищем минимальный элемент в оставшейся части массива
			if (min > arr[j]) { // Если найден элемент меньше текущего минимума
				i_min = j; // Обновляем индекс минимального элемента
				min = arr[j]; // Обновляем значение минимума
			}
		}

		// Обмен значениями: ставим найденный минимальный элемент на место текущего
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

//СОРТИРОВКА РАСЧЕСКОЙ
void CombSort(float* arr, int size) {
	int gap = size;  // Начальный шаг
	float ratio = 1.27f;  // Коэффициент уменьшения шага
	int flag = 1;  // Флаг для отслеживания обменов

	while (gap > 1 || flag) {
		gap /= ratio;  // Уменьшаем шаг
		if (gap < 1) {
			gap = 1;
		}
		flag = 0;
		for (int i = 0; i < size - gap; i++) {  // Сравниваем элементы и меняем их
			if (arr[i] > arr[i + gap]) {
				float t = arr[i];
				arr[i] = arr[i + gap];
				arr[i + gap] = t;
				flag = 1;  // Обмен произошел
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

//ПОРАЗРЯДНАЯ СОРТИРОВКА
void Fcount(float* arr, unsigned int counter[], int size, int offset) {
	unsigned char* b = (unsigned char*)arr; // Приводим указатель на массив float к unsigned char
	memset(counter, 0, 256 * sizeof(unsigned int)); // Обнуляем счетчик для 256 возможных значений байтов

	for (int i = 0; i < size; i++) {  // Подсчет количества вхождений каждого байта
		counter[b[i * sizeof(float) + offset]]++; // Увеличиваем счетчик для текущего байта
	}

	for (int i = 1; i < 256; i++) {  // Преобразуем счетчик в индексы для сортировки
		counter[i] += counter[i - 1]; // Суммируем значения для получения индексов
	}
}

void RadixSort(float* arr, int size) {
	unsigned int counter[256]; // Массив для подсчета вхождений байтов
	float* arr1 = (float*)malloc(size * sizeof(float)); // Выделяем память для временного массива
	if (arr1 == NULL) {
		return;
	}

	unsigned char* b = (unsigned char*)arr; // Приводим указатель на массив float к unsigned char

	for (int i = 0; i < sizeof(float); i++) {  // Основной цикл по каждому байту float (4 байта)
		Fcount(arr, counter, size, i); // Вызываем функцию Fcount для подсчета вхождений байтов

		for (int j = size - 1; j >= 0; j--) {   // Перемещение элементов в отсортированный массив mas1
			arr1[counter[b[j * sizeof(float) + i]] - 1] = arr[j]; // Распределяем элементы в массив mas1 по индексам из counter
			counter[b[j * sizeof(float) + i]]--; // Уменьшаем счетчик
		}
		// Копируем отсортированные элементы обратно в оригинальный массив mas
		memcpy(arr, arr1, size * sizeof(float)); // Копируем элементы из mas1 обратно в mas

	}
	free(arr1); // Освобождаем память, выделенную под массив mas1
}
void RadixSortWithSign(float* arr, int size) {
	float* pos = (float*)malloc(size * sizeof(float)); // Выделяем память для массива положительных и массива отрицательных чисел
	float* neg = (float*)malloc(size * sizeof(float));
	int posCount = 0, negCount = 0;
	if (pos == NULL || neg == NULL) {
		return;
	}
	// Разделяем положительные и отрицательные числа
	for (int i = 0; i < size; i++) {
		if (arr[i] >= 0) {
			pos[posCount++] = arr[i];
		}
		else {
			neg[negCount++] = arr[i];
		}
	}

	// Сортируем положительные и отрицательные числа отдельно
	RadixSort(neg, negCount);
	RadixSort(pos, posCount);

	// Объединяем отсортированные массивы
	for (int i = 0; i < negCount; i++) {
		arr[i] = neg[negCount - i - 1]; // Обратный порядок для отрицательных чисел
	}

	for (int i = 0; i < posCount; i++) {
		arr[negCount + i] = pos[i]; // Положительные числа идут после отрицательных
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


//ПОДСЧЕТ ВРЕМЕНИ РАБОТЫ ПРОГРАММЫ
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

//ПОДСЧЕТ СЛОЖНОСТИ РАБОТЫ АЛГОРИТМА
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

//ПОДТВЕРЖДЕНИЕ КОРРЕКТНОСТИ СОРТИРОВОК
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