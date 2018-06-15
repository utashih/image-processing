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

    // Assignment 3
    Bitmap grayscale(raw);
    grayscale.grayscale();
    grayscale.to_file("grayscale.bmp");

    Bitmap logarithmic(grayscale);
    logarithmic.logarithmic_operate();
    logarithmic.to_file("logarithmic.bmp");

    Bitmap histeq(grayscale);
    histeq.equalize_histogram();
    histeq.to_file("histeq.bmp");

    // Assignment 4
    Bitmap translated(raw);
    translated.translate(100, 200);
    translated.to_file("translated.bmp");

    Bitmap mirror_x(raw);
    mirror_x.mirror(Bitmap::Axis::x_axis);
    mirror_x.to_file("mirror_x.bmp");

    Bitmap mirror_y(raw);
    mirror_y.mirror(Bitmap::Axis::y_axis);
    mirror_y.to_file("mirror_y.bmp");

    Bitmap shrunk(raw);
    shrunk.scale(0.707);
    shrunk.to_file("shrunk.bmp");

    Bitmap enlarged(raw);
    enlarged.scale(1.414);
    enlarged.to_file("enlarged.bmp");

    Bitmap rotated(raw);
    rotated.rotate(M_PI / 6);
    rotated.to_file("rotated.bmp");

    Bitmap sheared_x(raw);
    sheared_x.shear(Bitmap::Axis::x_axis, 1.0);
    sheared_x.to_file("sheared_x.bmp");

    Bitmap sheared_y(raw);
    sheared_y.shear(Bitmap::Axis::y_axis, -0.5);
    sheared_y.to_file("sheared_y.bmp");
    return 0;
}