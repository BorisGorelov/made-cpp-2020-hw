#pragma once
#include <vector>
#include <utility>
#include <iostream>

class Line;
struct Point {
  double x, y;
  Point(): x(0), y(0) {}

  Point(double x, double y): x(x), y(y) {}

  bool operator==(const Point& another) const;

  bool operator!=(const Point& another) const;
  void rotate_point(const Point& center, double angle);
  void reflex_point_center(const Point& center);
  void reflex_point_axis(const Line& axis);
  void scale_point(const Point& center, double coefficient);

  friend std::ostream& operator<< (std::ostream &out, const Point &point);
  friend double dist(const Point& a, const Point& b);
  friend Point operator+(const Point& a, const Point& b);
  friend Point operator-(const Point& a, const Point& b);
  // friend bool intersect(const Line& m, const Line& n, Point* res);
  friend Line perpendicular(const Line& l, const Point& p);
};

class Line {
  // ax + by + c = 0
  double a;
  double b;
  double c;

 public:
  Line(): a(0), b(0), c(0) {}
  void normalize();

  Line(const Point& first, const Point& second);
  Line(double angular_coefficient, double shift): \
  a(-angular_coefficient), b(1.0), c(-shift) {normalize();}

  Line(const Point& p, double angular_coefficient);
  bool operator==(const Line& line);
  bool operator!=(const Line& line);
  friend bool intersect(const Line&, const Line&, Point*);
  friend Line perpendicular(const Line&, const Point&);
};

class Shape {
 public:
  virtual double perimeter() const = 0;
  virtual double area() const = 0;
  virtual bool operator==(const Shape& another) const = 0;
  virtual bool operator!=(const Shape& another) {
    return !(*this == another);
  }
  virtual void rotate(Point center, double angle) = 0;
  virtual void reflex(Point center) = 0;
  virtual void reflex(Line axis) = 0;
  virtual void scale(Point center, double coefficient) = 0;
  virtual ~Shape() {}
};

class Polygon: virtual public Shape {
 protected:
  std::vector<Point> points;

 public:
  Polygon(): points() {}
  explicit Polygon(const std::vector<Point>& vec): points(vec) {}
  int verticesCount() { return points.size(); }

  void print_polygon();
  std::vector<Point> getVertices() const;
  bool operator==(const Shape& another) const override;
  double perimeter() const override;
  double area() const override;

  void rotate(Point center, double angle);
  void reflex(Point center);
  void reflex(Line axis);
  void scale(Point center, double coefficient);
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

  Ellipse(const Point& f1, const Point& f2, double d);
  std::pair<Point, Point> focuses() const;
  double eccentricity() const;

  Point center() const;
  bool operator==(const Shape& another) const override;
  double perimeter() const override;
  double area() const override;
  void rotate(Point center, double angle);
  void reflex(Point center);
  void reflex(Line axis);
  void scale(Point center, double coefficient);
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

  bool operator==(const Shape& another) const override;

  double perimeter() const override;
  double area() const override;
  void rotate(Point center, double angle);
  void reflex(Point center);
  void reflex(Line axis);
  void scale(Point center, double coefficient);
};

class Rectangle: public Polygon {
 public:
  Rectangle(): Polygon() {}
  Rectangle(const Point& a, const Point& b, double ac_over_cb);
  Point center();
  std::pair<Line, Line> diagonals();
};

class Square: public Rectangle {
 public:
  Square(): Rectangle() {}
  Square(const Point& a, const Point& b);
  Circle circumscribedCircle();
  Circle inscribedCircle();
};

class Triangle: public Polygon {
 public:
  Triangle(): Polygon() {}
  Triangle(const Point& a, const Point& b, const Point& c);
  Point centroid();
  Circle circumscribedCircle();
  Circle inscribedCircle();
  Point orthocenter();
  Line EulerLine();
  Circle ninePointsCircle();
};
