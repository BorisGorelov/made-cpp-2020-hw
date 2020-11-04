#pragma once
#include <cmath>
#include <tr1/cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <utility>
#include <cassert>

const double EPS = 1e-7;

double det(double a, double b, double c, double d) {
  // |a b|
  // |c d|
  return a * d - c * b;
}

struct Point {
  double x, y;
  Point(): x(0), y(0) {}

  Point(double x, double y): x(x), y(y) {}

  bool operator==(const Point& another) const {
    return abs(x - another.x) < EPS && abs(y - another.y) < EPS;
  }

  bool operator!=(const Point& another) const {
    return !(*this == another);
  }

  friend std::ostream& operator<< (std::ostream &out, const Point &point);
};

std::ostream& operator<<(std::ostream &out, const Point &point) {
  out << "(" << point.x << ", " << point.y << "); ";
  return out;
}

Point operator+(const Point& a, const Point& b) {
  return {a.x + b.x, a.y + b.y};
}

Point operator-(const Point& a, const Point& b) {
  return {a.x - b.x, a.y - b.y};
}

Point points_center(const Point& a, const Point& b) {
  return {(a.x + b.x) / 2, (a.y + b.y) / 2};
}

double dist(const Point& a, const Point& b) {
  return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

class Line {
  // ax + by + c = 0
  double a;
  double b;
  double c;

 public:
  Line(): a(0), b(0), c(0) {}
  void normalize() {
    double z = sqrt(a * a + b * b);
    if (z > EPS) {
      a /= z;
      b /= z;
      c /= z;
    }
  }

  Line(const Point& first, const Point& second) {
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
    normalize();
  }

  // y = kx + b
  // -kx + y - b = 0
  Line(double angular_coefficient, double shift): \
  a(-angular_coefficient), b(1.0), c(-shift) {normalize();}

  // y - y1 = k(x - x1)
  // -kx + y - y1 + kx1 = 0
  Line(const Point& p, double angular_coefficient): \
    a(-angular_coefficient), \
    b(1.0), \
    c(angular_coefficient * p.x - p.y) {
    normalize();
  }

  bool operator==(const Line& line) {
    return abs(det(a, b, line.a, line.b)) < EPS && \
         abs(det(a, b, line.a, line.c)) < EPS && \
         abs(det(b, c, line.b, line.c)) < EPS;
  }

  bool operator!=(const Line& line) {
    return !(*this == line);
  }

  friend bool intersect(const Line&, const Line&, Point*);
  friend Line perpendicular(const Line&, const Point&);
};

bool intersect(const Line& m, const Line& n, Point* res) {
  double z = det(m.a, m.b, n.a, n.b);
  if (fabs(z) < EPS)
    return false;
  res->x = -det(m.c, m.b, n.c, n.b) / z;
  res->y = -det(m.a, m.c, n.a, n.c) / z;
  return true;
}

Line perpendicular(const Line& l, const Point& p) {
  Line ans;
  ans.a = -l.b;
  ans.b = l.a;
  ans.c = -l.a * p.y + l.b * p.x;
  return ans;
}

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
 protected:
  std::vector<Point> points;

 public:
  Polygon(): points() {}
  explicit Polygon(const std::vector<Point>& vec): points(vec) {}
  int verticesCount() { return points.size(); }

  void print_polygon() {
    for (auto i : points)
      std::cout << i;
    std::cout << '\n';
  }

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
    // shoelace formula
    std::vector<Point> buf{this->points};
    buf.push_back(buf[0]);
    for (int i = 1; i < buf.size(); ++i)
      ans += (buf[i].y * buf[i-1].x - buf[i].x * buf[i-1].y);
    return fabs(ans) / 2;
  }
};

class Ellipse: virtual public Shape {
 protected:
  Point focus1;
  Point focus2;
  double distance;
  double semi_major_axis;
  double semi_minor_axis;

 public:
  Ellipse(): focus1(0, 0), focus2(0, 0), distance(0), \
         semi_major_axis(0), semi_minor_axis(0) {}

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
          std::tr1::comp_ellint_2(this->eccentricity());
  }

  double area() const override {
    return M_PI * semi_major_axis * semi_minor_axis;
  }
};

class Circle: public Ellipse {
 protected:
  Point cent;
  double rad;

 public:
  Circle(const Point& center, double radius): \
  cent(center), rad(radius) {}

  double radius() const {return rad;}

  Point center() const {return cent;}

  bool operator==(const Shape& another) const override {
    const Circle* another_circle = dynamic_cast<const Circle*>(&another);
    if (!another_circle) {
      std::cerr << "Bad dynamic cast (shape to circle)\n";
      return false;
    }

    if (cent != another_circle->cent || rad != another_circle->rad)
      return false;
    else
      return true;
  }

  double perimeter() const override {
    return 2 * M_PI * rad;
  }

  double area() const override {
    return M_PI * rad * rad;
  }
};

class Rectangle: public Polygon {
 public:
  Rectangle(): Polygon() {}
  Rectangle(const Point& a, const Point& b, double ac_over_cb) {
    assert(ac_over_cb > EPS);
    if (a == b) {
      points = {a, a, a, a};
      return;
    }
    // alpha is slope of AB
    double alpha = atan2(b.y - a.y, b.x - a.x);

    ac_over_cb = std::min(ac_over_cb, 1 / ac_over_cb);
    // beta is slope of AC
    double beta = alpha + atan(ac_over_cb);

    double ab_magnitude = dist(a, b);
    double ac_magnitude = (ab_magnitude * ac_over_cb) / \
                          sqrt(ac_over_cb * ac_over_cb + 1);
    Point c, d;
    Point center = points_center(a, b);
    c.x = cos(beta) * ac_magnitude + a.x;
    c.y = sin(beta) * ac_magnitude + a.y;
    d.x = 2 * center.x - c.x;
    d.y = 2 * center.y - c.y;
    points = {a, c, b, d};
  }

  Point center() {
    Point a = points[0];
    Point b = points[2];
    return points_center(a, b);
  }

  std::pair<Line, Line> diagonals() {
    return {Line(points[0], points[2]), Line(points[1], points[3])};
  }
};

class Square: public Rectangle {
 public:
  Square(): Rectangle() {}
  Square(const Point& a, const Point& b) {
    if (a == b) {
      points = {a, a, a, a};
      return;
    }

    double alpha = atan2(b.y - a.y, b.x - a.x);
    double beta = alpha + M_PI / 4;
    double ab_magnitude = dist(a, b);
    double ac_magnitude = ab_magnitude / sqrt(2);
    Point c, d;
    Point center = {(a.x + b.x) / 2, (a.y + b.y) / 2};
    c.x = cos(beta) * ac_magnitude + a.x;
    c.y = sin(beta) * ac_magnitude + a.y;
    d.x = 2 * center.x - c.x;
    d.y = 2 * center.y - c.y;
    points = {a, c, b, d};
  }

  Circle circumscribedCircle() {
    return Circle(this->center(), sqrt(2) * dist(points[0], points[1]) / 2);
  }

  Circle inscribedCircle() {
    return Circle(this->center(), dist(points[0], points[1]) / 2);
  }
};

class Triangle: public Polygon {
 public:
  Triangle(): Polygon() {}
  Triangle(const Point& a, const Point& b, const Point& c) {
    points = {a, b, c};
  }

  Point centroid() {
    return {(points[0].x + points[1].x + points[2].x) / 3, \
            (points[0].y + points[1].y + points[2].y) / 3};
  }

  Circle circumscribedCircle() {
    Point a = points[0];
    Point b = points[1];
    Point c = points[2];
    Point* circumcircle_center = new Point();
    Point ab_middle = points_center(a, b);
    Point bc_middle = points_center(b, c);
    assert(intersect(perpendicular(Line(a, b), ab_middle), \
                     perpendicular(Line(b, c), bc_middle), \
                     circumcircle_center));
    double circumcircle_radius = (dist(a, b) * dist(b, c) * dist(c, a)) /\
                                 (4 * this->area());
    Circle ans = Circle(*circumcircle_center, circumcircle_radius);
    delete circumcircle_center;
    return ans;
  }

  Circle inscribedCircle() {
    // center - intersection of bisectors
    Point a = points[0];
    Point b = points[1];
    Point c = points[2];
    double a_bisector = (atan2(b.y - a.y, b.x - a.x) + \
                         atan2(c.y - a.y, c.x - a.x)) / 2;
    double b_bisector = (atan2(c.y - b.y, c.x - b.x) + \
                         atan2(a.y - b.y, a.x - b.x)) / 2;

    Point* incircle_center = new Point();
    assert(intersect(Line(a, tan(a_bisector)), \
                     Line(b, tan(b_bisector)), \
                     incircle_center));
    double incircle_radius = 2 * this->area() / this->perimeter();
    Circle ans = Circle(*incircle_center, incircle_radius);
    delete incircle_center;
    return ans;
  }

  Point orthocenter() {
    // OH = OA + OB + OC, O - circumcenter, H - orthocenter
    Point a = points[0];
    Point b = points[1];
    Point c = points[2];
    Point o = this->circumscribedCircle().center();
    a = a - o;
    b = b - o;
    c = c - o;
    return a + b + c + o;
  }

  Line EulerLine() {
    return Line(this->centroid(), this->orthocenter());
  }

  Circle ninePointsCircle() {
    Point a = points[0];
    Point b = points[1];
    Point c = points[2];
    return Triangle(points_center(a, b), \
                    points_center(b, c), \
                    points_center(c, a)).circumscribedCircle();
  }
};
