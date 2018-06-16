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

    // Assignment 6
    Bitmap bilateral(raw);
    bilateral.bilaterial_filter(12, 0.2);
    bilateral.to_file("bilateral.bmp");

    return 0;
}