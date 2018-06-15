#ifndef BITMAP_H__
#define BITMAP_H__

#include "bmp.h"
#include <functional>
#include <vector>

typedef std::vector<std::pair<int, int>> StructureElement;
typedef std::function<std::pair<int, int>(const int, const int)> Map;

class Bitmap {
  public:
    Bitmap();
    Bitmap(const Bitmap &) = default;
    ~Bitmap() = default;
    void from_file(const char *);
    void to_file(const char *);
    void to_file_binary(const char *);
    int get_width() const;
    int get_height() const;
    uint8_t get_threshold() const;

    void grayscale();

    void binarize(const uint8_t);
    void dilate(const StructureElement &);
    void erode(const StructureElement &);
    void open(const StructureElement &);
    void close(const StructureElement &);

    void logarithmic_operate();
    void equalize_histogram();

    enum class Axis { x_axis, y_axis };
    void transform(const Map &, const Map &);
    void translate(const int, const int);
    void mirror(const Axis);
    void scale(const double);
    void rotate(const double);
    void shear(const Axis, const double);

    class Kernel {
      public:
        using type = std::vector<double>;
        static type mean;
    };
    void filter(const Kernel::type &);

  private:
    bool initialized;
    BITMAPFILEHEADER file_header;
    BITMAPINFOHEADER info_header;
    std::vector<RGBQUAD> data;
    static uint8_t get_brightness(const RGBQUAD);
    static void set_rgb(RGBQUAD &, uint8_t r, uint8_t g, uint8_t b);
    static void set_brightness(RGBQUAD &, uint8_t);
    static uint8_t clamp(const double);
    template <typename T> static T clamp(const T value, const T lb, const T ub);
    static RGBQUAD new_RGBQUAD(const uint8_t, const uint8_t, const uint8_t);
    static std::vector<std::pair<int, int>> dir4, dir8;
};

#endif