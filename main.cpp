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
    mean.mean_filter();
    mean.to_file("mean.bmp");

    Bitmap laplacian(raw);
    laplacian.laplacian_enhance(0.5);
    laplacian.to_file("laplacian.bmp");

    return 0;
}