#include "Bitmap.h"
#include <cmath>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
    if (argc == 1) {
        std::cerr << "no input image" << std::endl;
        exit(1);
    }
    Bitmap raw;
    raw.from_file(argv[1]);

    // Assignment 5
    Bitmap mean(raw);
    mean.filter(Bitmap::Kernel::mean);

    return 0;
}