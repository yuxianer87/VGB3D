#include "Result.h"

std::ostream & operator<<(std::ostream & os,const Result  & r)
{
    return os
        << "id=" << r.id << std::endl
        << "a=" << r.a.transpose()<<std::endl
        << "b=" << r.b.transpose() << std::endl
        << "c=" << r.c.transpose() << std::endl
        << "a1=" << r.a1.transpose() << std::endl
        << "b1=" << r.b1.transpose() << std::endl
        << "c1=" << r.c1.transpose()<< std::endl;
}

void Result::create_id(void)
{
    std::string temp;
    to_string(temp, a);
    to_string(temp, b);
    to_string(temp, c);
    to_string(temp, a1);
    to_string(temp, b1);
    to_string(temp, c1);
    id = temp;
}

void Result::to_string(std::string &temp, Vector3i const& a)
{
    for (int i = 0; i < a.size(); i++)
    {
        temp.append(std::to_string(a[i]));
    }
}
