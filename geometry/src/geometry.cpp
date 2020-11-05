#pragma once
#include <cmath>
#include <tr1/cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <utility>
#include <cassert>
#include "./geometry.h"

const double EPS = 1e-6;

double det(double a, double b, double c, double d) {
  // |a b|
  // |c d|
  return a * d - c * b;
}

bool Point::operator==(const Point& another) const {
  return abs(x - another.x) < EPS && abs(y - another.y) < EPS;
}

bool Point::operator!=(const Point& another) const {
  return !(*this == another);
}

void Point::rotate_point(const Point& center, double angle) {
  // angle in radians
  double ao_dist = dist(*this, center);
  double alpha = atan2(this->y - center.y, this->x - center.x);
  this->x = cos(alpha + angle) * ao_dist + center.x;
  this->y = sin(alpha + angle) * ao_dist + center.y;
}

void Point::reflex_point_center(const Point& center) {
  Point ao = center - *this;
  *this = ao + ao + *this;
}

void Point::reflex_point_axis(const Line& axis) {
  Line perpendicular_to_axis = perpendicular(axis, *this);
  Point* intersection = new Point();
  assert(intersect(axis, perpendicular_to_axis, intersection));
  reflex_point_center(*intersection);
}

void Point::scale_point(const Point& center, double coefficient) {
  assert(fabs(coefficient) > EPS);
  Point ox = *this - center;
  Point ox_k = {ox.x * coefficient, ox.y * coefficient};
  *this = ox_k + center;
}

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

void Line::normalize() {
  double z = sqrt(a * a + b * b);
  if (z > EPS) {
    a /= z;
    b /= z;
    c /= z;
  }
}

Line::Line(const Point& first, const Point& second) {
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

// y - y1 = k(x - x1)
// -kx + y - y1 + kx1 = 0
Line::Line(const Point& p, double angular_coefficient) {
  a = -angular_coefficient;
  b = 1.0;
  c = angular_coefficient * p.x - p.y;
  normalize();
}

bool Line::operator==(const Line& line) {
  return abs(det(a, b, line.a, line.b)) < EPS && \
        abs(det(a, b, line.a, line.c)) < EPS && \
        abs(det(b, c, line.b, line.c)) < EPS;
}

bool Line::operator!=(const Line& line) {
  return !(*this == line);
}

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
  ans.normalize();
  return ans;
}

void Polygon::print_polygon() {
  for (auto i : points)
    std::cout << i;
  std::cout << '\n';
}

std::vector<Point> Polygon::getVertices() const {
  return points;
}

bool Polygon::operator==(const Shape& another) const {
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

double Polygon::perimeter() const {
  double ans = 0.0;
  for (int i = 1; i < this->points.size(); ++i)
    ans += dist(this->points[i - 1], this->points[i]);
  ans += dist(this->points[this->points.size() - 1], this->points[0]);
  return ans;
}

double Polygon::area() const {
  double ans = 0;
  // shoelace formula
  std::vector<Point> buf{this->points};
  buf.push_back(buf[0]);
  for (int i = 1; i < buf.size(); ++i)
    ans += (buf[i].y * buf[i-1].x - buf[i].x * buf[i-1].y);
  return fabs(ans) / 2;
}

void Polygon::rotate(Point center, double angle) {
  angle = angle * M_PI / 180;
  for (Point& i : points)
    i.rotate_point(center, angle);
}

void Polygon::reflex(Point center) {
  for (Point& i : points)
    i.reflex_point_center(center);
}

void Polygon::reflex(Line axis) {
  for (Point& i : points)
    i.reflex_point_axis(axis);
}

void Polygon::scale(Point center, double coefficient) {
  for (Point& i : points)
    i.scale_point(center, coefficient);
}

Ellipse::Ellipse(const Point& f1, const Point& f2, double d): \
focus1(f1), focus2(f2), distance(d) {
  semi_major_axis = distance / 2;
  semi_minor_axis = semi_major_axis * \
            sqrt(1 - pow(this->eccentricity(), 2));
}

std::pair<Point, Point> Ellipse::focuses() const {
  return std::make_pair(focus1, focus2);
}

double Ellipse::eccentricity() const {
  return dist(focus1, focus2) / distance;
}

Point Ellipse::center() const {
  return Point((focus1.x + focus2.x) / 2, (focus1.y + focus2.y) / 2);
}

bool Ellipse::operator==(const Shape& another) const {
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

double Ellipse::perimeter() const {
  return 4.0 * semi_major_axis * \
        std::tr1::comp_ellint_2(this->eccentricity());
}

double Ellipse::area() const {
  return M_PI * semi_major_axis * semi_minor_axis;
}

void Ellipse::rotate(Point center, double angle) {
  angle = angle * M_PI / 180;
  this->focus1.rotate_point(center, angle);
  this->focus2.rotate_point(center, angle);
}

void Ellipse::reflex(Point center) {
  this->focus1.reflex_point_center(center);
  this->focus2.reflex_point_center(center);
}

void Ellipse::reflex(Line axis) {
  this->focus1.reflex_point_axis(axis);
  this->focus2.reflex_point_axis(axis);
}

void Ellipse::scale(Point center, double coefficient) {
  this->focus1.scale_point(center, coefficient);
  this->focus2.scale_point(center, coefficient);
  this->distance *= coefficient;
  this->semi_major_axis *= coefficient;
  this->semi_minor_axis *= coefficient;
}

bool Circle::operator==(const Shape& another) const {
  const Circle* another_circle = dynamic_cast<const Circle*>(&another);
  if (!another_circle) {
    std::cerr << "Bad dynamic cast (shape to circle)\n";
    return false;
  }

  if (dist(cent, another_circle->cent) > EPS || \
      fabs(rad - another_circle->rad) > EPS)
    return false;
  else
    return true;
}

double Circle::perimeter() const {
  return 2 * M_PI * rad;
}

double Circle::area() const {
  return M_PI * rad * rad;
}

void Circle::rotate(Point center, double angle) {
  angle = angle * M_PI / 180;
  this->cent.rotate_point(center, angle);
}

void Circle::reflex(Point center) {
  this->cent.reflex_point_center(center);
}

void Circle::reflex(Line axis) {
  this->cent.reflex_point_axis(axis);
}

void Circle::scale(Point center, double coefficient) {
  this->cent.scale_point(center, coefficient);
  this->rad *= coefficient;
}

Rectangle::Rectangle(const Point& a, const Point& b, double ac_over_cb) {
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

Point Rectangle::center() {
  Point a = points[0];
  Point b = points[2];
  return points_center(a, b);
}

std::pair<Line, Line> Rectangle::diagonals() {
  return {Line(points[0], points[2]), Line(points[1], points[3])};
}

Square::Square(const Point& a, const Point& b) {
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

Circle Square::circumscribedCircle() {
  return Circle(this->center(), sqrt(2) * dist(points[0], points[1]) / 2);
}

Circle Square::inscribedCircle() {
  return Circle(this->center(), dist(points[0], points[1]) / 2);
}

Triangle::Triangle(const Point& a, const Point& b, const Point& c) {
  points = {a, b, c};
}

Point Triangle::centroid() {
  return {(points[0].x + points[1].x + points[2].x) / 3, \
          (points[0].y + points[1].y + points[2].y) / 3};
}

Circle Triangle::circumscribedCircle() {
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

Circle Triangle::inscribedCircle() {
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

Point Triangle::orthocenter() {
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

Line Triangle::EulerLine() {
  return Line(this->centroid(), this->orthocenter());
}

Circle Triangle::ninePointsCircle() {
  Point a = points[0];
  Point b = points[1];
  Point c = points[2];
  return Triangle(points_center(a, b), \
                  points_center(b, c), \
                  points_center(c, a)).circumscribedCircle();
}
