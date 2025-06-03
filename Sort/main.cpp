#include "Header.h"
int main() {
	int select, size;

	printf("Choise number sort:\n 1 - RadixSort\n 2 - SelectionSort\n 3 - MergeSort\n 4 - CombSort\n");

	printf("Sort - "); scanf_s("%d", &select);

	printf("Array size - "); scanf_s("%d", &size);

	switch (select) {
	case 1:

		ConfirmationOfCorrectness(RadixSortWithSign, size);
		ComplexityOfWorkRadix(RadixSortComplexityOfWork, RadixSortWithSignComplexityOfWork, size);
		TimeOfWork(RadixSortWithSign, size);

		break;
	case 2:

		ConfirmationOfCorrectness(SelectionSort, size);
		ComplexityOfWork(SelectionSortComplexityOfWork, size);
		TimeOfWork(SelectionSort, size);

		break;
	case 3:

		ConfirmationOfCorrectnessMerge(MergeSort, size);
		ComplexityOfWorkMerge(MergeSortComplexityOfWork, size);
		TimeOfWorkMerge(MergeSort, size);

		break;
	case 4:

		ConfirmationOfCorrectness(CombSort, size);
		ComplexityOfWork(CombSortComplexityOfWork, size);
		TimeOfWork(CombSort, size);


		break;
	default:
		printf("In the section <<Sort - >> enter a number from 1 to 4.");
	}
	return 0;
}
