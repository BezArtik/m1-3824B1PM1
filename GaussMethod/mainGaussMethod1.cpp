#include "HeaderGaussMethod1.h"

int main() {
    try {
        FinalRun();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}