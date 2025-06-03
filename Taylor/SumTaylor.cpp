#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define PI M_PI
#define COUNT_FUNCTION_TYPE 4
#define COUNT_SUMMATION_METHOD 3

//sin(x)
void ComputeTermsSin(float* ArrTerms, int count, float x) {
    if (count <= 0) {
        exit(1);
    }
    x = fmodf(x, 2.0f * (float)PI);
    ArrTerms[0] = x;
    for (int n = 1; n < count; n++) {
        ArrTerms[n] = ArrTerms[n - 1] * (-x * x) / ((2 * n + 1) * (2 * n));
    }
}

//cos(x)
void ComputeTermsCos(float* ArrTerms, int count, float x) {
    if (count <= 0) {
        exit(1);
    }
    x = fmodf(x, 2.0f * (float)PI);
    ArrTerms[0] = 1.0f;
    for (int n = 1; n < count; n++) {
        ArrTerms[n] = ArrTerms[n - 1] * (-x * x) / ((2 * n) * (2 * n - 1));
    }
}

//exp(x)
void ComputeTermsExp(float* ArrTerms, int count, float x) {
    if (count <= 0) {
        exit(1);
    }
    ArrTerms[0] = 1.0f;
    for (int n = 1; n < count; n++) {
        ArrTerms[n] = ArrTerms[n - 1] * x / n;
    }
}

//ln(1 + x)
void ComputeTermsLn(float* ArrTerms, int count, float x) {
    if (count <= 0) {
        exit(1);
    }
    ArrTerms[0] = x;
    for (int n = 1; n < count; n++) {
        ArrTerms[n] = ArrTerms[n - 1] * (-x) * n / (n + 1);
    }
}

//опълне ясллхпнбюмхе
float DirectSummation(float* ArrTerms, int count) {
    float SumRes = 0.0f;
    for (int i = 0; i < count; i++) {
        SumRes += ArrTerms[i];
    }
    return SumRes;
}

//напюрмне ясллхпнбюмхе
float ReverseSummation(float* ArrTerms, int count) {
    float SumRes = 0.0f;
    for (int i = count - 1; i >= 0; i--) {
        SumRes += ArrTerms[i];
    }
    return SumRes;
}

//оноюпмне ясллхпнбюмхе
float PairwiseSummation(float* ArrTerms, int start, int end) {
    if (start == end) {
        return ArrTerms[start];
    }

    int mid = (start + end) / 2;
    return PairwiseSummation(ArrTerms, start, mid) + PairwiseSummation(ArrTerms, mid + 1, end);
}

float PairwiseSummationRes(float* ArrTerms, int count) {
    if (count <= 0) {
        exit(1);
    }
    return PairwiseSummation(ArrTerms, 0, count - 1);
}

typedef float (*SummationFunction)(float* x, int count);

typedef void (*ComputeTermsFunction)(float* terms, int count, float x);

//бшвхякемхе ясллш пъдю
float ComputeTaylorRow(float x, int count, ComputeTermsFunction ComputeTerms, SummationFunction SummationMethod) {
    float* ArrTerms = (float*)malloc(sizeof(float) * count);
    if (ArrTerms == NULL) {
        exit(1);
    }

    ComputeTerms(ArrTerms, count, x);

    float Res = SummationMethod(ArrTerms, count);

    free(ArrTerms);
    return Res;
}

//бшбнд пегскэрюрнб
void CompareAndOutputResults(int FunctionType, int SummationMethod, float x, int count) {

    if (FunctionType <= 0 || SummationMethod <= 0 || FunctionType > COUNT_FUNCTION_TYPE || SummationMethod > COUNT_SUMMATION_METHOD) {
        printf("Invalid type of function or method of summation!");
        return;
    }
    //люяяхб тсмйжхи дкъ ондяверю вкемнб пъдю
    ComputeTermsFunction ComputeTermsFunctions[] = {
        NULL,
        ComputeTermsSin,  //sin(x)
        ComputeTermsCos,  //cos(x)
        ComputeTermsExp,  //exp(x)
        ComputeTermsLn    //ln(1 + x)
    };

    //люяяхб тсмйжхи дкъ ясллхпнбюмхъ
    SummationFunction SummationMethods[] = {
        NULL,
        DirectSummation,      //опълне ясллхпнбюмхе
        ReverseSummation,     //напюрмне ясллхпнбюмхе
        PairwiseSummationRes  //оноюпмне ясллхпнбюмхе
    };

    //люяяхб тсмйжхи хг math.h
    float (*MathLibFunctions[])(float x) = {
        NULL,
        sinf,   //sin(x)
        cosf,   //cos(x)
        expf,   //exp(x)
        log1pf  //ln(1 + x)
    };

    float MyResult = ComputeTaylorRow(x, count, ComputeTermsFunctions[FunctionType], SummationMethods[SummationMethod]);

    float MathLibResult = MathLibFunctions[FunctionType](x);

    printf("Result my method: %e\n", MyResult);
    printf("Result of library <math.h>: %e\n", MathLibResult);
    printf("Calculation error: %e\n", fabsf(MyResult - MathLibResult));

}