////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <string>
////////////////////////////////////////////////////////////////////////////////

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;

double inline det(const Point &u, const Point &v) {
	Point p1 = u;
	Point p2 = v; 
	double detr = p1.real() * p2.real() - p1.imag()*p2.imag();

	return detr;
}

struct Compare {
	Point p0; // Leftmost point of the poly
	bool operator ()(const Point &p1, const Point &p2) {
		// TODO
		return true;
	}
};

bool inline salientAngle(Point &a, Point &b, Point &c) {
	// TODO
	return false;
}

////////////////////////////////////////////////////////////////////////////////

Polygon convex_hull(std::vector<Point> &points) {
	Compare order;
	// TODO
	order.p0 = Point(0, 0);
	std::sort(points.begin(), points.end(), order);
	Polygon hull;
	// TODO
	// use salientAngle(a, b, c) here
	return hull;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Point> load_xyz(const std::string &filename) {
	std::vector<Point> points;
	
	std::ifstream in(filename);
	int num;
	
	in >> num;
	
	
	
	float a, b, z;
	if (in.is_open()) {
    while (in>> a >> b >> z){
        Point p1 = Point(a, b);
        points.push_back(p1);
		
    }
	}
	return points;
}

void save_obj(const std::string &filename, Polygon &poly) {
	std::ofstream out(filename);
	if (!out.is_open()) {
		throw std::runtime_error("failed to open file " + filename);
	}
	out << std::fixed;
	for (const auto &v : poly) {
		out << "v " << v.real() << ' ' << v.imag() << " 0\n";
	}
	for (size_t i = 0; i < poly.size(); ++i) {
		out << "l " << i+1 << ' ' << 1+(i+1)%poly.size() << "\n";
	}
	out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
	
	if (argc <= 2) {
		std::cerr << "Usage: " << argv[0] << " points.xyz output.obj" << std::endl;
	}
	
  	std::vector<Point> points = load_xyz(argv[1]);
	for(int i=0; i < points.size(); i++)
   		std::cout << points.at(i) << ' ' <<'\n';
	//Polygon hull = convex_hull(points);
	//save_obj(argv[2], hull);
	return 0;
}
