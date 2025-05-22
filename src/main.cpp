#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include <matplot/matplot.h>

#include <vector>
#include <complex>
#include <cmath>
#include <string>
#include <numeric>
#include <iomanip>

constexpr double PI = 3.14159265358979323846;

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace matplot;

// === Proste operacje matematyczne ===
int add(int i, int j) { return i + j; }
int subtract(int i, int j) { return i - j; }

// === Generowanie sygnałów ===
std::vector<double> generate_sin(double freq, double start, double end, size_t n_samples) {
    std::vector<double> signal; signal.reserve(n_samples);
    double dt = (end - start) / (n_samples - 1);
    for (size_t i = 0; i < n_samples; ++i) {
        double t = start + i * dt;
        signal.push_back(std::sin(2 * PI * freq * t));
    }
    return signal;
}

std::vector<double> generate_cos(double freq, double start, double end, size_t n_samples) {
    std::vector<double> signal; signal.reserve(n_samples);
    double dt = (end - start) / (n_samples - 1);
    for (size_t i = 0; i < n_samples; ++i) {
        double t = start + i * dt;
        signal.push_back(std::cos(2 * PI * freq * t));
    }
    return signal;
}

std::vector<double> generate_square(double freq, double start, double end, size_t n_samples) {
    std::vector<double> signal; signal.reserve(n_samples);
    double dt = (end - start) / (n_samples - 1);
    for (size_t i = 0; i < n_samples; ++i) {
        double t = start + i * dt;
        signal.push_back((std::sin(2 * PI * freq * t) >= 0) ? 1.0 : -1.0);
    }
    return signal;
}

std::vector<double> generate_sawtooth(double freq, double start, double end, size_t n_samples) {
    std::vector<double> signal; signal.reserve(n_samples);
    double dt = (end - start) / (n_samples - 1);
    for (size_t i = 0; i < n_samples; ++i) {
        double t = start + i * dt;
        double period = 1.0 / freq;
        signal.push_back(2.0 * (std::fmod(t, period) / period) - 1.0);
    }
    return signal;
}

// === DFT i IDFT ===
std::vector<std::complex<double>> dft(const std::vector<double>& x) {
    size_t N = x.size();
    std::vector<std::complex<double>> X(N);
    for (size_t k = 0; k < N; ++k) {
        std::complex<double> sum = 0;
        for (size_t n = 0; n < N; ++n) {
            double angle = -2 * PI * k * n / N;
            sum += x[n] * std::exp(std::complex<double>(0, angle));
        }
        X[k] = sum;
    }
    return X;
}

std::vector<double> idft(const std::vector<std::complex<double>>& X) {
    size_t N = X.size();
    std::vector<double> x(N);
    for (size_t n = 0; n < N; ++n) {
        std::complex<double> sum = 0;
        for (size_t k = 0; k < N; ++k) {
            double angle = 2 * PI * k * n / N;
            sum += X[k] * std::exp(std::complex<double>(0, angle));
        }
        x[n] = (sum.real() / static_cast<double>(N));
    }
    return x;
}

// === Filtracja 1D (IIR) ===
std::vector<double> filter1d(const std::vector<double>& x,
                             const std::vector<double>& b,
                             const std::vector<double>& a) {
    size_t N = x.size();
    size_t M = b.size();
    size_t L = a.size();
    std::vector<double> y(N, 0.0);
    for (size_t n = 0; n < N; ++n) {
        double acc = 0.0;
        for (size_t i = 0; i < M; ++i) {
            if (n >= i) acc += b[i] * x[n - i];
        }
        for (size_t j = 1; j < L; ++j) {
            if (n >= j) acc -= a[j] * y[n - j];
        }
        y[n] = acc / a[0];
    }
    return y;
}

// === Filtracja 2D (konwolucja) ===
std::vector<std::vector<double>> filter2d(
    const std::vector<std::vector<double>>& img,
    const std::vector<std::vector<double>>& kernel) {
    int h = img.size();
    int w = img[0].size();
    int kh = kernel.size();
    int kw = kernel[0].size();
    int ph = kh / 2;
    int pw = kw / 2;
    std::vector<std::vector<double>> out(h, std::vector<double>(w, 0.0));
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            double sum = 0.0;
            for (int m = 0; m < kh; ++m) {
                for (int n = 0; n < kw; ++n) {
                    int ii = i + m - ph;
                    int jj = j + n - pw;
                    if (ii >= 0 && ii < h && jj >= 0 && jj < w)
                        sum += img[ii][jj] * kernel[kh - 1 - m][kw - 1 - n];
                }
            }
            out[i][j] = sum;
        }
    }
    return out;
}

std::vector<std::vector<double>> gaussian_kernel(int size, double sigma) {
    std::vector<std::vector<double>> kernel(size, std::vector<double>(size));
    double sum = 0.0;
    int offset = size / 2;

    for (int i = -offset; i <= offset; ++i) {
        for (int j = -offset; j <= offset; ++j) {
            double value = std::exp(-(i * i + j * j) / (2 * sigma * sigma));
            kernel[i + offset][j + offset] = value;
            sum += value;
        }
    }

    // Normalizacja
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            kernel[i][j] /= sum;

    return kernel;
}

std::vector<std::vector<double>> gaussian_blur(const std::vector<std::vector<double>>& img, int kernel_size, double sigma) {
    auto kernel = gaussian_kernel(kernel_size, sigma);
    return filter2d(img, kernel);
}

// === Rysowanie: dowolny sygnał ===
void SimpleSignal(const std::vector<double>& x, const std::vector<double>& y) {
    plot(x, y);
    title("Signal");
    xlabel("X");
    ylabel("Y");
    show();
}

// === Wizualizacje sygnałów ===
void plot_sin(double freq, double start, double end, size_t n_samples) {
    auto t = std::vector<double>(n_samples);
    double dt = (end - start) / (n_samples - 1);
    for (size_t i = 0; i < n_samples; ++i) t[i] = start + i * dt;
    auto y = generate_sin(freq, start, end, n_samples);
    plot(t, y);
    title("Sygnał sinusoidalny");
    xlabel("Czas [s]");
    ylabel("Amplituda");
    show();
}

void plot_cos(double freq, double start, double end, size_t n_samples) {
    auto t = std::vector<double>(n_samples);
    double dt = (end - start) / (n_samples - 1);
    for (size_t i = 0; i < n_samples; ++i) t[i] = start + i * dt;
    auto y = generate_cos(freq, start, end, n_samples);
    plot(t, y);
    title("Sygnał cosinusoidalny");
    xlabel("Czas [s]");
    ylabel("Amplituda");
    show();
}

void plot_square(double freq, double start, double end, size_t n_samples) {
    auto t = std::vector<double>(n_samples);
    double dt = (end - start) / (n_samples - 1);
    for (size_t i = 0; i < n_samples; ++i) t[i] = start + i * dt;
    auto y = generate_square(freq, start, end, n_samples);
    plot(t, y);
    title("Sygnał prostokątny");
    xlabel("Czas [s]");
    ylabel("Amplituda");
    show();
}

void plot_sawtooth(double freq, double start, double end, size_t n_samples) {
    auto t = std::vector<double>(n_samples);
    double dt = (end - start) / (n_samples - 1);
    for (size_t i = 0; i < n_samples; ++i) t[i] = start + i * dt;
    auto y = generate_sawtooth(freq, start, end, n_samples);
    plot(t, y);
    title("Sygnał piłokształtny");
    xlabel("Czas [s]");
    ylabel("Amplituda");
    show();
}

// === Wizualizacja DFT i IDFT ===
void plot_dft(const std::vector<double>& signal, double sampling_rate) {
    auto X = dft(signal);
    size_t N = X.size();
    std::vector<double> freqs(N), mag(N);
    for (size_t k = 0; k < N; ++k) {
        freqs[k] = k * sampling_rate / static_cast<double>(N);
        mag[k] = std::abs(X[k]);
    }
    plot(freqs, mag);
    title("Widmo amplitudowe (DFT)");
    xlabel("Częstotliwość [Hz]");
    ylabel("|X(f)|");
    show();
}

void plot_idft(const std::vector<double>& signal, double sampling_rate) {
    auto X = dft(signal);
    auto x_rec = idft(X);
    size_t N = x_rec.size();
    std::vector<double> t(N);
    double dt = 1.0 / sampling_rate;
    for (size_t i = 0; i < N; ++i) t[i] = i * dt;
    plot(t, x_rec);
    title("Rekonstrukcja sygnału (IDFT)");
    xlabel("Czas [s]");
    ylabel("Amplituda");
    show();
}

// === Wizualizacja filtracji ===
void plot_filter1d(const std::vector<double>& signal,
                   const std::vector<double>& b,
                   const std::vector<double>& a,
                   double sampling_rate) {
    size_t N = signal.size();
    std::vector<double> t(N);
    double dt = 1.0 / sampling_rate;
    for (size_t i = 0; i < N; ++i) t[i] = i * dt;
    auto y = filter1d(signal, b, a);
    plot(t, signal); hold(on); plot(t, y, "--"); hold(off);
    title("Filtracja 1D");
    xlabel("Czas [s]");
    ylabel("Amplituda");
    legend({"oryginal", "filtrowany"});
    show();
}

#include <iomanip> // upewnij się, że to jest na początku pliku

void plot_filter2d(const std::vector<std::vector<double>>& img,
                   const std::vector<std::vector<double>>& kernel) {
    auto out = filter2d(img, kernel);

    // Wypisywanie przefiltrowanej macierzy
    std::cout << "Przefiltrowana macierz (Filtracja 2D):\n";
    for (const auto& row : out) {
        for (double val : row)
            std::cout << std::fixed << std::setprecision(4) << val << "\t";
        std::cout << "\n";
    }

    // Wizualizacja w formie spłaszczonej
    std::vector<double> x, y1, y2;
    for (const auto& row : img)
        for (double val : row) y1.push_back(val);
    for (const auto& row : out)
        for (double val : row) y2.push_back(val);

    x.resize(y1.size());
    std::iota(x.begin(), x.end(), 0);

    plot(x, y1, "-");
    hold(on);
    plot(x, y2, "--");
    hold(off);
    legend({"oryginal", "po filtracji"});
    title("Filtracja 2D (spłaszczone dane)");
    xlabel("Indeks próbki");
    ylabel("Wartość");
    show();
}


void plot_gaussian_blur(const std::vector<std::vector<double>>& img, int kernel_size, double sigma) {
    auto blurred = gaussian_blur(img, kernel_size, sigma);

    // Wypisywanie przefiltrowanej macierzy
    std::cout << "Przefiltrowana macierz (rozmycie Gaussa):\n";
    for (const auto& row : blurred) {
        for (double val : row)
            std::cout << std::fixed << std::setprecision(4) << val << "\t";
        std::cout << "\n";
    }

    // Przygotowanie danych do wykresu (spłaszczone)
    std::vector<double> x, y1, y2;
    for (const auto& row : img)
        for (double val : row) y1.push_back(val);

    for (const auto& row : blurred)
        for (double val : row) y2.push_back(val);

    x.resize(y1.size());
    std::iota(x.begin(), x.end(), 0);  // Indeksy: 0, 1, 2, ...

    using namespace matplot;
    plot(x, y1, "-");
    hold(on);
    plot(x, y2, "--");
    hold(off);
    legend({"oryginal", "po rozmyciu Gaussa"});
    title("Rozmycie Gaussa (2D spłaszczone)");
    xlabel("Indeks");
    ylabel("Wartość");
    show();
}


// === Rejestracja w module _core ===
PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Moduł C++ do generowania, analizowania (DFT/IDFT),
        filtracji (1D/2D) i wizualizacji sygnałów z pybind11 + matplot++
    )pbdoc";

    m.def("add", &add, "Dodaje dwie liczby całkowite");
    m.def("subtract", &subtract, "Odejmuje dwie liczby całkowite");

    m.def("generate_sin", &generate_sin, py::arg("freq"), py::arg("start"), py::arg("end"), py::arg("n_samples"));
    m.def("generate_cos", &generate_cos, py::arg("freq"), py::arg("start"), py::arg("end"), py::arg("n_samples"));
    m.def("generate_square", &generate_square, py::arg("freq"), py::arg("start"), py::arg("end"), py::arg("n_samples"));
    m.def("generate_sawtooth", &generate_sawtooth, py::arg("freq"), py::arg("start"), py::arg("end"), py::arg("n_samples"));

    m.def("SimpleSignal", &SimpleSignal, py::arg("x"), py::arg("y"));
    m.def("plot_sin", &plot_sin, py::arg("freq"), py::arg("start"), py::arg("end"), py::arg("n_samples"));
    m.def("plot_cos", &plot_cos, py::arg("freq"), py::arg("start"), py::arg("end"), py::arg("n_samples"));
    m.def("plot_square", &plot_square, py::arg("freq"), py::arg("start"), py::arg("end"), py::arg("n_samples"));
    m.def("plot_sawtooth", &plot_sawtooth, py::arg("freq"), py::arg("start"), py::arg("end"), py::arg("n_samples"));

    m.def("dft", &dft, py::arg("signal"));
    m.def("idft", &idft, py::arg("spectrum"));
    m.def("plot_dft", &plot_dft, py::arg("signal"), py::arg("sampling_rate"));
    m.def("plot_idft", &plot_idft, py::arg("signal"), py::arg("sampling_rate"));

    m.def("filter1d", &filter1d, py::arg("signal"), py::arg("b"), py::arg("a"));
    m.def("plot_filter1d", &plot_filter1d, py::arg("signal"), py::arg("b"), py::arg("a"), py::arg("sampling_rate"));
    m.def("filter2d", &filter2d, py::arg("img"), py::arg("kernel"));
    m.def("plot_filter2d", &plot_filter2d, py::arg("img"), py::arg("kernel"));

    m.def("gaussian_blur", &gaussian_blur, py::arg("img"), py::arg("kernel_size"), py::arg("sigma"));
    m.def("plot_gaussian_blur", &plot_gaussian_blur, py::arg("img"), py::arg("kernel_size"), py::arg("sigma"));


#ifdef VERSION_INFO
    m.attr("version") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("version") = "dev";
#endif
}
