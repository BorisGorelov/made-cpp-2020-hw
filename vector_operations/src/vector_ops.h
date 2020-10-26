#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>

namespace task {
const double EPS1 = 1e-7;

double operator*(const std::vector<double>& v1, \
                 const std::vector<double>& v2) {
  assert(v1.size() == v2.size());
  double buf = 0;
  for (int i = 0; i < v1.size(); ++i)
    buf += v1[i] * v2[i];
  return buf;
}

std::vector<double> operator+(const std::vector<double>& v1, \
                              const std::vector<double>& v2) {
  assert(v1.size() == v2.size());
  std::vector<double> buf;
  for (int i = 0; i < v1.size(); ++i)
      buf.push_back(v1[i] + v2[i]);
  return buf;
}

std::vector<double> operator-(const std::vector<double>& v1, \
                              const std::vector<double>& v2) {
  assert(v1.size() == v2.size());
  std::vector<double> buf;
  for (int i = 0; i < v1.size(); ++i)
    buf.push_back(v1[i] - v2[i]);
  return buf;
}

std::vector<double> operator-(const std::vector<double>& v) {
  std::vector<double> buf;
  for (int i = 0; i < v.size(); ++i)
    buf.push_back(-v[i]);
  return buf;
}

std::vector<double> operator+(const std::vector<double>& v) {
  std::vector<double> buf;
  for (int i = 0; i < v.size(); ++i)
    buf.push_back(v[i]);
  return buf;
}

std::vector<double> operator%(const std::vector<double>& v1, \
                              const std::vector<double>& v2) {
  assert(v1.size() == v2.size());
  assert(v1.size() == 3);
  std::vector<double> buf(3);
  buf[0] = v1[1] * v2[2] - v1[2] * v2[1];
  buf[1] = v1[2] * v2[0] - v1[0] * v2[2];
  buf[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return buf;
}

bool operator||(const std::vector<double>& v1, \
                const std::vector<double>& v2) {
  assert(v1.size() == v2.size());
  double v1_len = sqrt(v1 * v1);
  double v2_len = sqrt(v2 * v2);
  double dot_product = v1 * v2;
  return (fabs(fabs(dot_product) - fabs(v2_len * v1_len)) < EPS1);
}

bool operator&&(const std::vector<double>& v1, \
                const std::vector<double>& v2) {
  assert(v1.size() == v2.size());
  double dot_product = v1 * v2;
  return ((v1 || v2) && (dot_product >= 0));
}

void reverse(std::vector<double>& v) {
  for (int i = 0; i < v.size() / 2; ++i)
    std::swap(v[v.size() - 1 - i], v[i]);
}

std::vector<int> operator|(const std::vector<int>& v1, \
                           const std::vector<int>& v2) {
  assert(v1.size() == v2.size());
  std::vector<int> buf;
  for (int i = 0; i < v1.size(); ++i) {
    buf.push_back(v1[i] | v2[i]);
  }
  return buf;
}

std::vector<int> operator&(const std::vector<int>& v1, \
                           const std::vector<int>& v2) {
  assert(v1.size() == v2.size());
  std::vector<int> buf;
  for (int i = 0; i < v1.size(); ++i) {
    buf.push_back(v1[i] & v2[i]);
  }
  return buf;
}

std::ostream& operator<< (std::ostream &out, \
                          const std::vector<double>& v) {
  for (auto i : v)
    out << i << ' ';
  out << '\n';
  return out;                  
}

std::istream& operator>> (std::istream &in, \
                          std::vector<double>& v) {
  int size;
  double buf;
  v.resize(0);
  in >> size;
  if (size)
    for (int i = 0; i < size; ++i) {
      in >> buf;
      v.push_back(buf);
    }
  return in;
}
}  // namespace task
