#pragma once
#include <vector>
#include <iostream>
#include <cmath>
//double EPS = 1e-6;

namespace task {

double operator*(const std::vector<double>& v1, const std::vector<double>& v2)
{
    if (v1.size() != v2.size()){
        throw -1;
    }
    double buf = 0;
    for (int i = 0; i < v1.size(); ++i) {
        buf += v1[i] * v2[i];
    }
    return buf;
}

std::vector<double> operator+(const std::vector<double>& v1, const std::vector<double>& v2)
{
    if (v1.size() != v2.size()){
        throw -1;
    }
    std::vector<double> buf;
    for (int i = 0; i < v1.size(); ++i) {
        buf.push_back(v1[i] - v2[i]);
    }
    return buf;
}

std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2)
{
    if (v1.size() != v2.size()){
        throw -1;
    }
    std::vector<double> buf;
    for (int i = 0; i < v1.size(); ++i) {
        buf.push_back(v1[i] + v2[i]);
    }
    return buf;
}

std::vector<double> operator-(const std::vector<double>& v) {
    std::vector<double> buf;
    for (int i = 0; i < v.size(); ++i) {
        buf.push_back(-v[i]);
    }
    return buf;
}

std::vector<double> operator+(const std::vector<double>& v) {
    std::vector<double> buf;
    for (int i = 0; i < v.size(); ++i) {
        buf.push_back(v[i]);
    }
    return buf;
}

std::vector<double> operator%(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size() && v1.size() != 3)
        throw -1;
    std::vector<double> buf;
    buf[0] = v1[1] * v2[2] - v1[2] * v2[1]; 
    buf[1] = v1[2] * v2[0] - v1[0] * v2[2]; 
    buf[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return buf; 
}

bool operator||(const std::vector<double>& v1, const std::vector<double>& v2)
{
    if (v1.size() != v2.size()){
        throw -1;
    }
    double v1_len = v1 * v1;
    double v2_len = v2 * v2;
    double dot_product = v1 * v2;
    return (fabs(fabs(dot_product) - fabs(v2_len * v1_len)) < EPS); 
}

bool operator&&(const std::vector<double>& v1, const std::vector<double>& v2)
{
    if (v1.size() != v2.size()){
        throw -1;
    }
    double dot_product = v1 * v2;
    return ((v1 || v2) && (dot_product >= 0));
}

std::vector<double> reverse(const std::vector<double>& v) {
    std::vector<double> buf(v.size());
    for (int i = 0; i <= (int)v.size() / 2; ++i) {
        buf[i] = v[i];
        buf[v.size() - 1 - i] = v[v.size() - 1 - i];
    }
    if (v.size() % 2 == 1) {
        buf[v.size() / 2 + 1] = v[v.size() / 2 + 1];
    }
    return buf;
}

std::vector<int> operator|(const std::vector<int>& v1, const std::vector<int>& v2) {
    std::vector<int> buf;
    if (v1.size() != v2.size()){
        throw -1;
    }
    for (int i = 0; i < v1.size(); ++i){
        buf.push_back(v1[i] | v2[i]);
    }
    return buf;
}

std::vector<int> operator&(const std::vector<int>& v1, const std::vector<int>& v2) {
    std::vector<int> buf;
    if (v1.size() != v2.size()){
        throw -1;
    }
    for (int i = 0; i < v1.size(); ++i){
        buf.push_back(v1[i] & v2[i]);
    }
    return buf;
}
}  // namespace task