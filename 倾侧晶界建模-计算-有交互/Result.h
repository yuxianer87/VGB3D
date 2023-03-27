#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <string>
using namespace Eigen;

struct Result
{
    friend  std::ostream& operator<<(std::ostream& os, const Result  & r);
    void create_id(void);
    void to_string(std::string &temp, Vector3i const& a);
public:
    Vector3i a;
    Vector3i b;
    Vector3i c;

    Vector3i a1;
    Vector3i b1;
    Vector3i c1;
    std::string id;
};

