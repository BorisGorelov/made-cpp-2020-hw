#pragma once
#include <cmath>
#include <vector>
#include <algorithm>

const double EPS = 1e-6;

double det(double a, double b, double c, double d) {
    // |a b|
    // |c d|
    return a * d - c * b;
}

struct Point {
    double x, y;
    Point(double x, double y): x(x), y(y) {}
    bool operator==(const Point& another) const {
        return abs(x - another.x) < EPS && abs(y - another.y) < EPS;
    }
};

class Line {
    // ax + by + c = 0
    double a;
    double b;
    double c;

 public:
    Line(Point first, Point second) {
        if (first.x == second.x) {
            // x = x1
            a = 1.0;
            b = 0.0;
            c = -first.x;
        } else {
            a = -(second.y - first.y) / (second.x - first.x);
            b = 1.0;
            c = (first.x * second.y - first.y * second.x) / \
            (second.x - first.x);
        }
    }

    // y = kx + b
    // -kx + y - b = 0
    Line(double angular_coefficient, double shift): \
    a(-angular_coefficient), b(1.0), c(-shift) {}

    // y - y1 = k(x - x1)
    // -kx + y - y1 + kx1 = 0
    Line(Point p, double angular_coefficient): \
    a(-angular_coefficient), \
    b(1.0), \
    c(angular_coefficient * p.x - p.y) {}

    bool operator==(Line line) {
        return abs(det(a, b, line.a, line.b)) < EPS && \
               abs(det(a, b, line.a, line.c)) < EPS && \
               abs(det(b, c, line.b, line.c)) < EPS;
    }

    bool operator!=(Line line) {
        return !(*this == line);
    }
};

class Shape {
 public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator==(const Shape& another) const = 0;
    virtual bool operator!=(const Shape& another) const {
        return !(*this == another);
    }
};

class Polygon: virtual public Shape {
    std::vector<Point> points;

 public:
    explicit Polygon(const std::vector<Point>& vec): points(vec) {}
    int verticesCount() { return points.size(); }

    std::vector<Point> getVertices() const {
        return points;
    }

    bool operator==(const Shape& another) const override {
        const Polygon* another_polygon = dynamic_cast<const Polygon*>(&another);
        if (!another_polygon) {
            std::cerr << "Bad dynamic cast (shape to polygon)\n";
            return false;
        }

        if (another_polygon->points.size() != this->points.size())
            return false;

        // buffer has duplicated vector of points of another polygon
        std::vector<Point> buf(another_polygon->points);
        buf.insert(std::end(buf), std::begin(another_polygon->points), \
        std::end(another_polygon->points));
        std::vector<Point> reversed_buf(buf);
        std::reverse(std::begin(reversed_buf), std::end(reversed_buf));
        int counter = 0;
        int max_counter = 0;
        int reversed_counter = 0;
        int reversed_max_counter = 0;
        for (int i = 0; i < buf.size(); ++i) {
            if (this->points[counter] == buf[i]) {
                counter += 1;
                if (counter > max_counter)
                    max_counter = counter;
                if (max_counter == this->points.size() - 1)
                    return true;
            } else {
                counter = 0;
            }
            if (this->points[reversed_counter] == reversed_buf[i]) {
                reversed_counter += 1;
                if (reversed_counter > reversed_max_counter)
                    reversed_max_counter = reversed_counter;
                if (reversed_max_counter == this->points.size() - 1)
                    return true;
            } else {
                reversed_counter = 0;
            }
        }
        return false;
    }
};
