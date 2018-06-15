#include "Bitmap.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
using std::size_t;

Bitmap::Bitmap() : initialized(false) {}

uint8_t Bitmap::clamp(const double x) {
    if (x > 255.0)
        return 255;
    if (x < 0.0)
        return 0;
    return static_cast<uint8_t>(x);
}

template <typename T>
T Bitmap::clamp(const T value, const T lower_bound, const T upper_bound) {
    if (value >= upper_bound)
        return upper_bound;
    if (value <= lower_bound)
        return lower_bound;
    return value;
}

uint8_t Bitmap::get_brightness(const RGBQUAD quad) {
    uint8_t brightness = clamp(0.299 * quad.rgbRed + 0.587 * quad.rgbGreen +
                               0.114 * quad.rgbBlue);
    return brightness;
}

void Bitmap::set_rgb(RGBQUAD &quad, uint8_t r, uint8_t g, uint8_t b) {
    quad.rgbBlue = b;
    quad.rgbGreen = g;
    quad.rgbRed = r;
    quad.rgbReserved = 0;
}

void Bitmap::set_brightness(RGBQUAD &quad, uint8_t brightness) {
    quad.rgbBlue = brightness;
    quad.rgbGreen = brightness;
    quad.rgbRed = brightness;
    quad.rgbReserved = 0;
}

std::vector<std::pair<int, int>> Bitmap::dir4{{0, 1}, {1, 0}, {0, -1}, {-1, 0}};

std::vector<std::pair<int, int>> Bitmap::dir8{
    {0, 1}, {1, 0}, {0, -1}, {-1, 0}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1}};

void Bitmap::from_file(const char *filename) {
    if (initialized) {
        throw std::runtime_error("Bitmap already initialized");
    }
    initialized = true;
    data.clear();

    std::ifstream fin(filename, std::ios::binary);
    fin.read(reinterpret_cast<char *>(&file_header), 14);
    fin.read(reinterpret_cast<char *>(&info_header), 40);

    RGBQUAD temp;
    char alignment[4];
    int W = get_width();
    int H = get_height();
    int offset = (4 - (W * 3) % 4) % 4;
    size_t quadSize = info_header.biBitCount / 8;

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            fin.read(reinterpret_cast<char *>(&temp), quadSize);
            data.push_back(temp);
        }
        fin.read(alignment, offset);
    }
    fin.close();
}

void Bitmap::to_file(const char *filename) {
    if (!initialized) {
        throw std::runtime_error("Bitmap not initialized");
    }

    std::ofstream fout(filename);
    fout.write(reinterpret_cast<const char *>(&file_header), 14);
    fout.flush();
    fout.write(reinterpret_cast<const char *>(&info_header), 40);
    fout.flush();

    const char alignment[4] = {0, 0, 0, 0};
    int W = get_width();
    int H = get_height();
    int offset = (4 - (W * 3) % 4) % 4;
    size_t quadSize = info_header.biBitCount / 8;

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            fout.write(reinterpret_cast<char *>(&data[y * W + x]), quadSize);
        }
        fout.write(alignment, offset);
        fout.flush();
    }
    fout.close();
}

void Bitmap::to_file_binary(const char *filename) {
    if (!initialized) {
        throw std::runtime_error("Bitmap not initialized");
    }
    for (auto &quad : data) {
        set_brightness(quad, quad.rgbBlue == 255 ? 0 : 255);
    }
    to_file(filename);
    for (auto &quad : data) {
        set_brightness(quad, quad.rgbBlue == 255 ? 0 : 255);
    }
}

int Bitmap::get_width() const {
    if (!initialized) {
        throw std::runtime_error("Bitmap not initialized");
    }
    return info_header.biWidth;
}

int Bitmap::get_height() const {
    if (!initialized) {
        throw std::runtime_error("Bitmap not initialized");
    }
    return info_header.biHeight;
}

void Bitmap::grayscale() {
    for (auto &quad : data) {
        uint8_t brightness = get_brightness(quad);
        set_brightness(quad, brightness);
    }
}

uint8_t Bitmap::get_threshold() const {
    uint8_t min_bri = std::accumulate(
        data.begin(), data.end(), std::numeric_limits<uint8_t>::max(),
        [](uint8_t min_yet, const RGBQUAD &quad) -> uint8_t {
            return std::min(min_yet, get_brightness(quad));
        });

    uint8_t max_bri = std::accumulate(
        data.begin(), data.end(), std::numeric_limits<uint8_t>::min(),
        [](uint8_t max_yet, const RGBQUAD &quad) -> uint8_t {
            return std::max(max_yet, get_brightness(quad));
        });

    uint8_t threshold = min_bri;
    double threshold_metric = std::numeric_limits<double>::min();

    for (uint8_t bri = min_bri; bri <= max_bri; bri++) {
        double metric = [this](const uint8_t bri) -> double {
            int NFgrd = std::count_if(data.begin(), data.end(),
                                      [&](const RGBQUAD &quad) -> bool {
                                          return get_brightness(quad) >= bri;
                                      });
            int NBgrd = std::count_if(data.begin(), data.end(),
                                      [&](const RGBQUAD &quad) -> bool {
                                          return get_brightness(quad) < bri;
                                      });
            int N = data.size();
            double wf = static_cast<double>(NFgrd) / N;
            double wb = static_cast<double>(NBgrd) / N;

            double muf =
                std::accumulate(data.begin(), data.end(), 0.0,
                                [&](double sum, const RGBQUAD &quad) -> double {
                                    uint8_t brightness = get_brightness(quad);
                                    if (brightness >= bri)
                                        return sum + static_cast<double>(
                                                         get_brightness(quad));
                                    else
                                        return sum;
                                }) /
                NFgrd;

            double mub =
                std::accumulate(data.begin(), data.end(), 0.0,
                                [&](double sum, const RGBQUAD &quad) -> double {
                                    uint8_t brightness = get_brightness(quad);
                                    if (brightness < bri)
                                        return sum + static_cast<double>(
                                                         get_brightness(quad));
                                    else
                                        return sum;
                                }) /
                NBgrd;
            return wf * wb * (muf - mub) * (muf - mub);
        }(bri);
        if (metric > threshold_metric) {
            threshold = bri;
            threshold_metric = metric;
        }
        if (bri == 255)
            break;
    }
    return threshold;
}

void Bitmap::binarize(const uint8_t threshold) {
    for (auto &quad : data) {
        if (quad.rgbBlue >= threshold) {
            set_brightness(quad, 0);
        } else {
            set_brightness(quad, 255);
        }
    }
}

void Bitmap::dilate(const StructureElement &elem) {
    std::vector<RGBQUAD> original(data);

    for (auto &quad : data) {
        set_brightness(quad, 0);
    }

    int W = get_width();
    int H = get_height();

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            bool ok = false;
            for (auto [dx, dy] : elem) {
                if (y + dy >= 0 && y + dy < H && x + dx >= 0 && x + dx < W) {
                    ok |= (original[(y + dy) * W + (x + dx)].rgbBlue == 255);
                }
            }
            set_brightness(data[y * W + x], ok ? 255 : 0);
        }
    }
}

void Bitmap::erode(const StructureElement &elem) {
    std::vector<RGBQUAD> original(data);

    for (auto &quad : data) {
        set_brightness(quad, 0);
    }

    int W = get_width();
    int H = get_height();

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            bool ok = true;
            for (auto [dx, dy] : elem) {
                if (y + dy >= 0 && y + dy < H && x + dx >= 0 && x + dx < W) {
                    ok &= (original[(y + dy) * W + (x + dx)].rgbBlue == 255);
                }
            }
            set_brightness(data[y * W + x], ok ? 255 : 0);
        }
    }
}

void Bitmap::open(const StructureElement &elem) {
    erode(elem);
    dilate(elem);
}

void Bitmap::close(const StructureElement &elem) {
    dilate(elem);
    erode(elem);
}

void Bitmap::logarithmic_operate() {
    uint8_t max_bri = std::accumulate(
        data.begin(), data.end(), std::numeric_limits<uint8_t>::min(),
        [](uint8_t max_yet, const RGBQUAD &quad) -> uint8_t {
            return std::max(max_yet, get_brightness(quad));
        });
    for (auto &quad : data) {
        uint8_t brightness = get_brightness(quad);
        uint8_t new_brightness = clamp(255 * std::log(brightness / 255.0 + 1) /
                                       std::log(max_bri / 255.0 + 1));
        set_brightness(quad, new_brightness);
    }
}

void Bitmap::equalize_histogram() {
    std::vector<int> histogram(256, 0);
    for (auto &quad : data)
        histogram[get_brightness(quad)] += 1;
    std::partial_sum(histogram.begin(), histogram.end(), histogram.begin());
    int total_cnt = data.size();
    std::vector<uint8_t> brightness(256);
    for (int i = 0; i < 256; i++)
        brightness[i] =
            clamp(255 * static_cast<double>(histogram[i]) / total_cnt);
    for (auto &quad : data)
        set_brightness(quad, brightness[get_brightness(quad)]);
}

void Bitmap::transform(const Map &map, const Map &inv_map) {
    std::vector<RGBQUAD> original(std::move(data));
    int W = get_width();
    int H = get_height();
    int x_max = std::numeric_limits<int>::min(),
        x_min = std::numeric_limits<int>::max(),
        y_max = std::numeric_limits<int>::min(),
        y_min = std::numeric_limits<int>::max();
    for (int k = 0; k < 4; k++) {
        int x = k & 1 ? W - 1 : 0;
        int y = k & 2 ? H - 1 : 0;
        auto [x1, y1] = map(x, y);
        x_max = std::max(x_max, x1);
        x_min = std::min(x_min, x1);
        y_max = std::max(y_max, y1);
        y_min = std::min(y_min, y1);
    }
    int W1 = x_max - x_min + 1;
    int H1 = y_max - y_min + 1;
    info_header.biWidth = W1;
    info_header.biHeight = H1;
    data.resize(W1 * H1);
    int x_offset = -x_min;
    int y_offset = -y_min;
    for (int y1 = 0; y1 < H1; y1++) {
        for (int x1 = 0; x1 < W1; x1++) {
            RGBQUAD &quad = data[y1 * W1 + x1];
            auto [x, y] = inv_map(x1 - x_offset, y1 - y_offset);
            if (y >= 0 && y < H && x >= 0 && x < W)
                quad = original[y * W + x];
        }
    }
}

void Bitmap::translate(const int dx, const int dy) {
    if (dx < 0 || dy < 0) {
        throw std::range_error("cannot translate with negative offset");
    }
    std::vector<RGBQUAD> original(std::move(data));
    int W = get_width();
    int H = get_height();
    int W1 = W + dx;
    int H1 = H + dy;
    info_header.biWidth = W1;
    info_header.biHeight = H1;
    for (int y1 = 0; y1 < H1; y1++) {
        for (int x1 = 0; x1 < W1; x1++) {
            if (y1 >= dy && x1 >= dx) {
                data.push_back(original[(y1 - dy) * W + (x1 - dx)]);
            } else {
                data.push_back(RGBQUAD());
            }
        }
    }
}

void Bitmap::mirror(const Axis axis) {
    int W = get_width();
    int H = get_height();
    if (axis == Bitmap::Axis::x_axis) {
        for (int y = 0; y < H; y++) {
            for (int x = 0, z = W - 1; x < z; x++, z--) {
                std::swap(data[y * W + x], data[y * W + z]);
            }
        }
    } else if (axis == Bitmap::Axis::y_axis) {
        for (int y = 0, z = H - 1; y < z; y++, z--) {
            for (int x = 0; x < W; x++) {
                std::swap(data[y * W + x], data[z * W + x]);
            }
        }
    } else {
        throw std::runtime_error("unknown axis");
    }
}

void Bitmap::scale(const double ratio) {
    if (ratio <= 0.0) {
        throw std::runtime_error("non-positive scaling ratio");
    }
    std::vector<RGBQUAD> original(std::move(data));
    int W = get_width();
    int H = get_height();
    int W1 = static_cast<int>(W * ratio);
    int H1 = static_cast<int>(H * ratio);
    info_header.biWidth = W1;
    info_header.biHeight = H1;
    data.resize(W1 * H1);
    for (int y1 = 0; y1 < H1; y1++) {
        for (int x1 = 0; x1 < W1; x1++) {
            double yd = (y1 + 0.5) / ratio - 0.5;
            double xd = (x1 + 0.5) / ratio - 0.5;
            int yf = static_cast<int>(std::floor(yd));
            int yc = static_cast<int>(std::ceil(yd));
            int xf = static_cast<int>(std::floor(xd));
            int xc = static_cast<int>(std::ceil(xd));
            std::vector<std::pair<RGBQUAD, double>> accumulated;
            if (yf >= 0 && yf < H && xf >= 0 && xf < W) {
                accumulated.push_back(std::make_pair(original[yf * W + xf],
                                                     (xc - xd) * (yc - yd)));
            }
            if (yf >= 0 && yf < H && xc >= 0 && xc < W) {
                accumulated.push_back(std::make_pair(original[yf * W + xc],
                                                     (xd - xf) * (yc - yd)));
            }
            if (yc >= 0 && yc < H && xf >= 0 && xf < W) {
                accumulated.push_back(std::make_pair(original[yc * W + xf],
                                                     (xc - xd) * (yd - yf)));
            }
            if (yc >= 0 && yc < H && xc >= 0 && xc < W) {
                accumulated.push_back(std::make_pair(original[yc * W + xc],
                                                     (xd - xf) * (yd - yf)));
            }
            double Z = std::accumulate(
                accumulated.begin(), accumulated.end(), 0.0,
                [](const double sum, const std::pair<RGBQUAD, double> pair)
                    -> double { return sum + pair.second; });
            uint8_t r =
                clamp(std::accumulate(
                          accumulated.begin(), accumulated.end(), 0.0,
                          [](const double sum,
                             const std::pair<RGBQUAD, double> pair) -> double {
                              return sum + pair.first.rgbRed * pair.second;
                          }) /
                      Z);
            uint8_t g =
                clamp(std::accumulate(
                          accumulated.begin(), accumulated.end(), 0.0,
                          [](const double sum,
                             const std::pair<RGBQUAD, double> pair) -> double {
                              return sum + pair.first.rgbGreen * pair.second;
                          }) /
                      Z);
            uint8_t b =
                clamp(std::accumulate(
                          accumulated.begin(), accumulated.end(), 0.0,
                          [](const double sum,
                             const std::pair<RGBQUAD, double> pair) -> double {
                              return sum + pair.first.rgbBlue * pair.second;
                          }) /
                      Z);
            set_rgb(data[y1 * W1 + x1], r, g, b);
        }
    }
}

void Bitmap::rotate(const double theta) {
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);
    auto map = [=](const int x, const int y) -> std::pair<int, int> {
        int x1 = static_cast<int>(cos_theta * x - sin_theta * y);
        int y1 = static_cast<int>(sin_theta * x + cos_theta * y);
        return std::make_pair(x1, y1);
    };
    auto inv_map = [=](const int x, const int y) -> std::pair<int, int> {
        int x1 = static_cast<int>(cos_theta * x + sin_theta * y);
        int y1 = static_cast<int>(-sin_theta * x + cos_theta * y);
        return std::make_pair(x1, y1);
    };
    transform(map, inv_map);
}

void Bitmap::shear(const Axis axis, const double ratio) {
    if (axis == Bitmap::Axis::x_axis) {
        auto map = [=](const int x, const int y) -> std::pair<int, int> {
            int x1 = static_cast<int>(x + ratio * y);
            return std::make_pair(x1, y);
        };
        auto inv_map = [=](const int x, const int y) -> std::pair<int, int> {
            int x1 = static_cast<int>(x - ratio * y);
            return std::make_pair(x1, y);
        };
        transform(map, inv_map);
    } else if (axis == Bitmap::Axis::y_axis) {
        auto map = [=](const int x, const int y) -> std::pair<int, int> {
            int y1 = static_cast<int>(y + ratio * x);
            return std::make_pair(x, y1);
        };
        auto inv_map = [=](const int x, const int y) -> std::pair<int, int> {
            int y1 = static_cast<int>(y - ratio * x);
            return std::make_pair(x, y1);
        };
        transform(map, inv_map);
    } else {
        throw std::runtime_error("unknown axis");
    }
}

void Bitmap::filter(const Bitmap::Kernel::type &kernel) {
    std::vector<RGBQUAD> original(data);
    int W = get_width();
    int H = get_height();
    double Z = std::accumulate(kernel.begin(), kernel.end(), 0.0);
    for (int y = 1; y < H - 1; y++) {
        for (int x = 1; x < W - 1; x++) {
            double r = 0.0, g = 0.0, b = 0.0;
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    double w = kernel[dy * 3 + dx + 4];
                    RGBQUAD &quad = original[(y + dy) * W + (x + dx)];
                    r += w * quad.rgbRed;
                    g += w * quad.rgbGreen;
                    b += w * quad.rgbBlue;
                }
            }
            set_rgb(data[y * W + x], clamp(r / Z), clamp(g / Z), clamp(b / Z));
        }
    }
}

Bitmap::Kernel::type Bitmap::Kernel::mean(9, 1.0);