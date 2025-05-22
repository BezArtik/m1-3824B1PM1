#include "HeaderGaussMethod.h"

enum TypeChoice {
    FloatType = 1,    
    DoubleType = 2
};

int main() {
    try {
        std::cout << "Select numeric type:\n";
        std::cout << "1. float (single precision)\n";
        std::cout << "2. double (double precision)\n";
        std::cout << "3. Exit\n";
        std::cout << "Your choice: ";
        size_t type_choice;
        std::cin >> type_choice;

        switch (type_choice) {
        case FloatType:
            Run<float>();
            break;
        case DoubleType:
            Run<double>();
            break;
        case Exit:
            std::cout << "Exiting program...\n";
            return 0;
        default:
            throw std::invalid_argument("Invalid choice");
        }        
    }
    catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}