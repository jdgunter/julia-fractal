#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include <algorithm>

#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <omp.h>

#include <png++/png.hpp>

// defining a few constants
#define MAX_ITERATIONS 500
#define XPIXELS 3840
#define YPIXELS 2160

typedef std::complex<double> complex;

// code for this taken from
// https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
class InputParser {
public:
    InputParser(int& argc, char **argv) {
	for (int i = 1; i < argc; ++i) {
	    this->tokens.push_back(std::string(argv[i]));
	}
    }

    const std::string& getCmdOption(const std::string& option) const{
	std::vector<std::string>::const_iterator iter;
	iter = std::find(this->tokens.begin(), this->tokens.end(), option);
	if (iter != this->tokens.end() && ++iter != this->tokens.end()) {
	    return *iter;
	}
	static const std::string empty_string("");
	return empty_string;
    }

    std::vector<std::string> getCmdOptionN(const std::string& option, int n) const{
	std::vector<std::string>::const_iterator iter;
	std::vector<std::string> results;
	std::vector<std::string> empty_vec;
	iter = std::find(this->tokens.begin(), this->tokens.end(), option);
	if (iter == this->tokens.end()) {
	    return empty_vec;
	}
        for (int i = 1; i < n+1; i++) {
	    if (iter+i != this->tokens.end()) {
		results.push_back(*(iter+i));
	    } else {
		return empty_vec;
	    }
	}
	return results;
    }

    bool cmdOptionExists(const std::string& option) const{
	return std::find(this->tokens.begin(), this->tokens.end(), option) != this->tokens.end();
    }

private:
    std::vector<std::string> tokens;
};

// struct representing vectors in R^3
// along with some basic operations on them
struct Vec3 {
    double x;
    double y;
    double z;

    Vec3() {}
    Vec3(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}
    
    Vec3 operator+(const Vec3& v) {
	return Vec3(x+v.x, y+v.y, z+v.z);
    }

    Vec3 operator-(const Vec3& v) {
	return Vec3(x-v.x, y-v.y, z-v.z);
    }

    Vec3 operator*(double& a) {
	return Vec3(a*x, a*y, a*z);
    }

    Vec3 operator/(double& a) {
	return Vec3(x/a, y/a, z/a);
    }
    
    bool operator==(const Vec3& v) {
	return x == v.x && y == v.y && z == v.z;
    }

    void print() {
	std:: cout << "(" << x << ", " << y << ", " << z << ")\n";
    }
    
    double dist(const Vec3& u, const Vec3& v) {
	return sqrt((u.x - v.x)*(u.x - v.x) +
		    (u.y - v.y)*(u.y - v.y) +
		    (u.z - v.z)*(u.z - v.z));
    }
};

Vec3 operator*(double& a, const Vec3& v)
{
    return Vec3(a*v.x, a*v.y, a*v.z);
};

// counts the number of iterations it takes for a complex function `f(z) = z^power + c0` evaluated iteratively
// at an initial point `init` to grow greater than 2 in magnitude
// normalized to achieve smoother coloring, look at this webpage for details:
// http://linas.org/art-gallery/escape/escape.html
double normalized_iterations(complex init, complex c0, int power)
{
    complex z = init;
    int iters = 0;
    while ((abs(z) <= 2) && (iters < MAX_ITERATIONS)) {
	z = std::pow(z,power) + c0;
	iters += 1;
    }
    double mu = iters;
    if ( iters < MAX_ITERATIONS ) {
	mu = iters + 1 - log(log(abs(z))) / log(power);
    }
    return mu;
}

// computes v + t(u - v)
// t should be a value between 0 and 1
Vec3 linear_interpolation(Vec3& v, Vec3& u, double t)
{
    return v + t*(u - v);
}

// creates a linear gradient of SIZE colours, using RGB values from PTS
// interspersed evenly
std::vector<Vec3> linear_interpolated_gradient(std::vector<Vec3> pts, int size)
{
    std::vector<Vec3> pal;
    int to_travel = size;
    int lines_left = pts.size();
    int pts_to_color;
    for (int i = 0; i < pts.size()-1; i++) {
	if (to_travel % lines_left != 0) {
	    pts_to_color = (to_travel / lines_left)+1;
	} else {
	    pts_to_color = to_travel / lines_left;
	}
	to_travel = to_travel - pts_to_color;
	lines_left--;
	double scaling = 1.0 / pts_to_color;
	Vec3 delta_vec = scaling*(pts[i+1] - pts[i]);
	Vec3 next_color = pts[i];
	for (int j = 0; j < pts_to_color; j++) {
	    pal.push_back(next_color);
	    next_color = next_color + delta_vec;
	}
    }
    return pal;
}
    

int main(int argc, char *argv[])
{
    const std::string& usage = "Usage: -f <filename> [-p <power>] -c <real_part> <imag_part> [-origin <x> <y>] [-z <zoom>] [-verbose]\nPower defaults to 2, origin defaults to (0,0)\n";

    // Parsing command line arguments
    InputParser input(argc, argv);
    const std::string& filename = input.getCmdOption("-f");
    if (filename.empty()) {
	std::cout << usage;
	return 0;
    }
    
    const std::string& power_string = input.getCmdOption("-p");
    int power = 2;
    if (!power_string.empty()) {
	power = stoi(power_string);
    }
    
    const std::vector<std::string>& complex_strings = input.getCmdOptionN("-c", 2);
    if (complex_strings.empty()) {
	std::cout << usage;
	return 0;
    }
    const double real_part = stod(complex_strings[0]);
    const double imag_part = stod(complex_strings[1]);

    double origin_x = 0.0, origin_y = 0.0;
    const std::vector<std::string>& origin_strings = input.getCmdOptionN("-origin", 2);
    if (!origin_strings.empty()) {
        origin_x = stod(origin_strings[0]);
	origin_y = stod(origin_strings[1]);
    }

    double zoom = 1.0;
    const std::string& zoom_string = input.getCmdOption("-z");
    if (!zoom_string.empty()) {
	zoom = stod(zoom_string);
    }
    
    bool verbose = input.cmdOptionExists("-verbose");

    // Setting up parameters
    const complex complex_constant(real_part, imag_part);

    // computing C -> pixel mapping
    double im_start = origin_y + 1.08/zoom;
    double re_start = origin_x - 1.92/zoom;
    double delta_y = 2*std::abs(im_start) / YPIXELS, delta_x = 2*std::abs(re_start) / XPIXELS;
    double im, re;

    if (verbose) {
	std::cout << "im_start = " << im_start << "\nre_start = " << re_start << std::endl;
	std::cout << "delta_y = " << delta_y << "\ndelta_x = " << delta_x << std::endl;
	std::cout << "zoom = " << zoom << std::endl;
	std::cout << "Running on " << omp_get_max_threads() << " threads" << std::endl;
    }

    // another thing that would be nice to add is allow the user to input a file
    // consisting of RGB triples to set up the color palette with
    std::vector<Vec3> colors;
    colors.push_back(Vec3(9, 15, 113));
    colors.push_back(Vec3(233, 221, 236));
    colors.push_back(Vec3(242, 164, 58));  
    
    std::vector<Vec3> palette = linear_interpolated_gradient(colors, 100);
    png::image<png::rgb_pixel> image(XPIXELS, YPIXELS);
    
    #pragma omp parallel for private(re) private(im)
    for (int y = 0; y < YPIXELS; y++) {
	if (verbose) {
	    std::cout << "Computing row " << y+1 << '/' << YPIXELS << "...\n";
	}
	im = im_start - y*delta_y;
	for (int x = 0; x < XPIXELS; x++) {
	    re = re_start + x*delta_x;
	    complex init(re,im);
	    double mu = normalized_iterations(init, complex_constant, power);
	    // scale mu to be in the range of 1-100
	    mu *= 100.0/MAX_ITERATIONS;
	    double tmp;
	    Vec3 color1 = palette[(int)floor(mu)];
	    Vec3 color2 = palette[(int)ceil(mu)];
	    Vec3 color = linear_interpolation(color1, color2, modf(mu, &tmp));
	    image[y][x] = png::rgb_pixel(color.x, color.y, color.z);
	}
    }
    image.write(filename);
    return 0;
}
