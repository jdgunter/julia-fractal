#include <iostream>
#include <array>
#include <complex>
#include <vector>
#include <string>
#include <algorithm>

#include <cmath>
#include <cstdlib>

#include <omp.h>
#include <png++/png.hpp>


// Added recommendations given in
// https://codereview.stackexchange.com/questions/188095/julia-fractal-drawing-in-c


// defining a few constants
static constexpr auto max_iterations = 2000;
static constexpr auto width = 3840;
static constexpr auto height = 2160;

typedef std::complex<double> complex;

// code for this taken from
// https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
class InputParser {
public:
    InputParser(int argc, char **argv) {
	for (int i = 1; i < argc; ++i) {
	    tokens.emplace_back(std::string(argv[i]));
	}
    }

    const std::string& optionValue(const std::string& option) const {
	static const std::string not_found{};
	auto iter = std::find(tokens.begin(), tokens.end(), option);
	return iter != tokens.end() && ++iter != tokens.end() ? *iter : not_found;
    }

    std::vector<std::string> optionValues(const std::string& option, int n) const {
	static const std::vector<std::string> not_found{};
	auto iter = std::find(tokens.begin(), tokens.end(), option);
	if (std::distance(iter, tokens.end()) <= n) return not_found;

	std::vector<std::string> results;
	results.reserve(n);
	while (n--) results.push_back(*++iter);
	return results;
    }

    bool contains(const std::string& option) const {
	return std::find(tokens.begin(), tokens.end(), option) != tokens.end();
    }

private:
    std::vector<std::string> tokens;
};

// counts the number of iterations it takes for a complex function `f(z) = z^power + c` evaluated iteratively
// at an initial point `init` to grow greater than 2 in magnitude
// normalized to achieve smoother coloring, look at this webpage for details:
// http://linas.org/art-gallery/escape/escape.html
double normalized_iterations(complex z, complex c, int power)
{
    int iters;
    for (iters = 0; std::norm(z) <= 4; ++iters) {
	if (iters == max_iterations) return iters;
	z = std::pow(z,power) + c;
    }
    return iters + 1 - std::log(std::log(std::abs(z))) / std::log(power);
}

png::rgb_pixel linear_interpolation(const png::rgb_pixel& v, const png::rgb_pixel& u, double a)
{
    auto b = 1 - a;
    return png::rgb_pixel(b*v.red + a*u.red,
			  b*v.green + a*u.green,
			  b*v.blue + a*u.blue);
}

int main(int argc, char *argv[])
{
    static const auto usage =
	"Usage: -f <filename> [-p <power>] -c <real_part> <imag_part>"
        " [-origin <x> <y>] [-z <zoom>] [-verbose]\n"
        "Power defaults to 2, origin defaults to (0,0)\n";
    
    // Parsing command line arguments
    InputParser input(argc, argv);

    if (input.contains("-h") || input.contains("--help")) {
	std::cout << usage;
	return 0;
    }
    
    const bool verbose = input.contains("-verbose");
    const auto filename = input.optionValue("-f");
    const auto power_string = input.optionValue("-p");
    const auto complex_strings = input.optionValues("-c", 2);
    const auto origin_strings = input.optionValues("-origin", 2);
    const auto zoom_string = input.optionValue("-z");
    
    int power = 2;
    if (!power_string.empty()) {
	power = std::stoi(power_string);
    }
    
    if (complex_strings.empty()) {
	std::cerr << usage;
	return 1;
    }
    const complex complex_constant{std::stod(complex_strings[0]),
	                           std::stod(complex_strings[1])};

    double origin_x = 0.0, origin_y = 0.0;
    if (!origin_strings.empty()) {
        origin_x = std::stod(origin_strings[0]);
	origin_y = std::stod(origin_strings[1]);
    }

    double zoom = 1.0;
    if (!zoom_string.empty()) {
	zoom = std::stod(zoom_string);
    }

    std::ofstream outfile;
    std::ostream& out = filename.empty()
	? std::cout
	: (outfile.open(filename, std::ios::out|std::ios::trunc|std::ios::binary), outfile);

    if (!out) {
	perror(filename.c_str());
	return 1;
    }
    
    // setting up Julia parameters
    // computing C -> pixel mapping
    double im_start = origin_y + 1.08/zoom;
    double re_start = origin_x - 1.92/zoom;
    double delta_y = 2*1.08/zoom / height;
    double delta_x = 2*1.92/zoom / width;
    double im, re;

    if (verbose) {
	std::clog << "im_start = " << im_start << "\nre_start = " << re_start << std::endl;
	std::clog << "delta_y = " << delta_y << "\ndelta_x = " << delta_x << std::endl;
	std::clog << "zoom = " << zoom << std::endl;
	std::clog << "Running on " << omp_get_max_threads() << " threads" << std::endl;
    }

    // another thing that would be nice to add is allow the user to input a file
    // consisting of RGB triples to set up the color palette with
    static const std::vector<png::rgb_pixel> colors{
	{  0,   0,   0},
	{213,  67,  31},
	{251, 255, 121},
	{ 62, 223,  89},
	{ 43,  30, 218},
	{  0, 255, 247}
    };
    static const auto max_colors = colors.size()-1;

    png::image<png::rgb_pixel> image(width, height);
    
#pragma omp parallel for private(re) private(im)
    for (int y = 0; y < height; y++) {
	if (verbose)
#pragma omp critical
	{
	    std::clog << "Computing row " << y+1 << '/' << height << "...\n";
	}
	double im = im_start - y*delta_y;
	for (int x = 0; x < width; x++) {
	    re = re_start + x*delta_x;
	    double mu = normalized_iterations({re, im}, complex_constant, power);
	    // scale mu to be in the range of colours
	    mu *= static_cast<double>(max_colors)/max_iterations;
	    auto i_mu = static_cast<std::size_t>(mu);
	    auto color1 = colors[i_mu];
	    auto color2 = colors[std::min(i_mu+1, max_colors)];
	    double tmp;
	    image[y][x] = linear_interpolation(color1, color2, modf(mu, &tmp));
	}
    }
    
    image.write_stream(out);
    return 0;
}
