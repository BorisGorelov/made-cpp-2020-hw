#pragma once
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <utility>

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

double dist(const Point& a, const Point& b) {
    return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y-b.y));
}

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

    double perimeter() const override {
        double ans = 0.0;
        for (int i = 1; i < this->points.size(); ++i)
            ans += dist(this->points[i - 1], this->points[i]);
        ans += dist(this->points[this->points.size() - 1], this->points[0]);
        return ans;
    }

    double area() const override {
        double ans = 0;
        // sholace formula
        // pivot - (0. 0)
        std::vector<Point> buf{this->points};
        buf.push_back(buf[0]);
        for (int i = 1; i < buf.size(); ++i)
            ans += (buf[i].y * buf[i-1].x - buf[i].x * buf[i-1].y);
        return fabs(ans) / 2;
    }
};

class Ellipse: virtual public Shape {
    Point focus1;
    Point focus2;
    double distance;
    double semi_major_axis;
    double semi_minor_axis;

 public:
    Ellipse(const Point& f1, const Point& f2, double d): \
    focus1(f1), focus2(f2), distance(d) {
        semi_major_axis = distance / 2;
        semi_minor_axis = semi_major_axis * \
                          sqrt(1 - pow(this->eccentricity(), 2));
    }

    std::pair<Point, Point> focuses() const {
        return std::make_pair(focus1, focus2);
    }

    double eccentricity() const {
        return dist(focus1, focus2) / distance;
    }

    Point center() const {
        return Point((focus1.x + focus2.x) / 2, (focus1.y + focus2.y) / 2);
    }

    bool operator==(const Shape& another) const override {
        const Ellipse* another_ellipse = dynamic_cast<const Ellipse*>(&another);
        if (!another_ellipse) {
            std::cerr << "Bad dynamic cast (shape to ellipse)\n";
            return false;
        }

        if (distance != another_ellipse->distance)
            return false;

        if (this->focuses() == another_ellipse->focuses() || \
            this->focuses() == std::make_pair(another_ellipse->focus2, \
            another_ellipse->focus1)) {
            return true;
        }
        return false;
    }

    double perimeter() const override {
        return 4.0 * semi_major_axis * \
               std::comp_ellint_2(this->eccentricity());
    }

    double area() const override {
        double pi = atan(1) * 4;
        return pi * semi_major_axis * semi_minor_axis;
    }
};

class Circle: public Ellipse {
    Point center;
    double radius;

 public:
    Circle(const Point& center, double radius): \
    center(center), radius(radius) {}

    double radius()

};
