#include "HeaderTaylor.h"


int main() {
    int SelectFunction;
    int SelectMethod;
    float x;
    int count;

    printf("Enter the number function:\n 1 - sin(x)\n 2 - cos(x)\n 3 - exp(x)\n 4 - ln(1+x)\n");
    printf("Function - "); scanf_s("%d", &SelectFunction);

    printf("Enter the number method summation: \n 1 - Direct Summation\n 2 - Reverse Summation\n 3 - Pairwise Summation\n");
    printf("Method - "); scanf_s("%d", &SelectMethod);

    printf("Enter the value: "); scanf_s("%f", &x);

    printf("Enter the number of terms: "); scanf_s("%d", &count);

    CompareAndOutputResults(SelectFunction, SelectMethod, x, count);

    return 0;
}
