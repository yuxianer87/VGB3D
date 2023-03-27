#include <iostream>
#include <Eigen/Dense>
#include <list>
#include <cassert>
#include<cmath>
#include <algorithm>
#include <cmath>
#include "Result.h"
#include "log.h"
#include <vector>
#include <unordered_map>
#include <set>
#include<fstream>
#include <map>

using namespace std;
using namespace Eigen;


//https://eigen.tuxfamily.org/index.php?title=Main_Page
const float pai = 3.141592653589793;

//角度转换成弧度
float radian(float degree)
{
    //cout << "sin 60=" << sin(3.141592653589793 / 6) << endl;
    return degree * pai / 180;
}
float arc_radian(float degree)
{
    return degree * 180 / pai;
}
typedef Vector3i vector3;
typedef Vector2i point_i;
vector3 get_b(void)
{
    vector3 v;
    cout << "请输入b（3个整数bi in [0,1,2,3,4,5]，中间用空格隔开,）：";
    cin >> v(0) >> v(1) >> v(2);
    cout << "你输入的b=" << v.transpose() << endl;

    return v;
}
int get_N(vector3 const& b)
{
    return b.squaredNorm();
}
int get_segema(void)
{
    int segema = 0;
    cout << "请输入segema, segema in [0, 99]. segema=";
    cin >> segema;
    return segema;
}
std::ostream& operator<<(std::ostream& os, Vector2i const& p)
{
    return os << "[" << p[0] << "," << p[1] << "]";
}
std::vector<Vector2i> get_x0_y0(int segema, int N)
{
    cout << "遍历求解可能的x0([1,N]=[1," << N << "])和y0([1,x0]):" << endl;
    std::vector<Vector2i> result;
    Vector2i temp;
    temp << -1, -1;
    int x0 = 0;//x0 in [1, segema]; x0 in Z
    int y0 = 0;//y0 in [1, segema]; y0 in Z
    int start = 1;
    int last = segema;
    for (size_t i = start; i <= last; i++)
    {
        for (size_t j = start; j <= last; j++)
        {
            x0 = i;
            y0 = j;
            ///cout << "x0=" << x0 << ", y0=" << y0;
            if (x0*x0+y0*y0*N == segema)
            {
                temp[0] = x0;
                temp[1] = y0;
                //cout << "根据segema和N,找到一组符合条件的[x0,y0]=[" 
                //    << x0 << "," << y0 << "]" << endl;
               // cout<< " OK." << endl;
                result.push_back(temp);
            }
           /* else
            {
                cout << " not OK." << endl;
            }
            */
        }
    }
    cout << "可以使用的[x0,y0]为：" << endl;
    if (result.size() == 0)
    {
        cout << "没有可用的x0,y0" << endl;
    }
    for (auto& item : result)
    {
        cout << "x0=" << item[0] << ", y0=" << item[1] << endl;
    }
    return result;
}
bool is_x0_y0_valid(Vector2i x0y0)
{
    if (x0y0[0] == -1 || x0y0[1] == -1 || x0y0[1] == 0)
    {
        return false;
    }
    if (!(x0y0[0] > x0y0[1]))
    {
        return false;
    }
    return true;
}
float get_seta_degree(Vector2i x0y0, float N)
{
    pl("根据 N=" << N<<", x0="<<x0y0[0]<<", y0="<<x0y0[1]<<"开始计算seta:");
    auto y0_div_x0 = static_cast<float>( x0y0[1]) / x0y0[0];
    pl("y0/x0=" << y0_div_x0);
    auto t = y0_div_x0 * std::sqrt(N);
    pl("y0/x0 * std::sqrt(N) = " << t);
    float seata_degree = 2 * atan(t);
    return seata_degree;
}
bool my_abs_ok(Vector3i const& a)
{
    auto t1 = std::abs(a[0]);
    auto t2 = std::abs(a[1]);
    auto t3 = std::abs(a[2]);
    bool ret = (t1 < 10) && (t2 < 10) && (t3 < 10);
	return ret;
}
long gcd_i(long x, long y)
{
	auto x1 = x, y1 = y;
	if (x == 0)
	{
		return y;
	}
	if (y == 0)
	{
		return x;
	}
	long lRem = 0;
	while (y != 0)
	{
		lRem = x % y;
		x = y;
		y = lRem;
	}
	//pl("gcd_i " << x1 << ", " << y1 << "="<<x);
	return x;
}

int gcd_3i(Vector3i& result)
{
	auto copy = result;
	int a = result[0];
	for (size_t i = 0; i < result.size(); i++)
	{
		a = gcd_i(a, result[i]);
	}
	a = a < 0 ? a * (-1) : a;
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = result[i] / a;
	}
	//cout << "约去最大公约数： " << copy.transpose() << "=>" << result.transpose() << endl;
	return a;
}
void gcd(Result& result)
{
	gcd_3i(result.a);
	gcd_3i(result.b);
	gcd_3i(result.c);
	gcd_3i(result.a1);
	gcd_3i(result.b1);
	gcd_3i(result.c1);
}
std::map<std::string, Result> Processing(Eigen::Vector3d &b, float seata_degree, int &seagema)
{
    cout << "*************begin****************" << endl;
    cout << "b=" << b.transpose();
    cout << ", seagema=" << seagema;
    cout << ", seata_degree=" << seata_degree;
    cout << "b的单位向量[d,e,f]=" << endl;
    auto n = b.normalized();//normalized 单位向量
    auto d = n[0];

    auto e = n[1];
    auto f = n[2];
    cout << "d=" << d << ", e=" << e << ", f=" << f << endl;
    auto g = cos(radian(seata_degree));
    auto h = sin(radian(seata_degree));
    auto i = 1.0 - g;
    Matrix3f Rseata;
    Rseata << i * pow(d, 2) + g, i*d*e - h * f, i*d*f + h * e,
        i*d*e + h * f, i*pow(e, 2) + g, i*d*f - h * d,
        i*d*f - h * e, i*d*f + h * d, i*pow(f, 2) + g;
    cout << "Rseta=" << Rseata << endl;
    auto R = seagema * Rseata;
    cout << "R=" << R << endl;
    Matrix3i iR;


    //一直R b求a a1 c c1，使得
   /*
    已知：
    R*a = a1;
    R*b = b1;
    R*c = c1;
    求：a,c,a1,c1
    a,b,c 两两正交
    a,b,c 坐标[-10, 10]
    a,b,c 坐标约去公约数
    */


    for (size_t i = 0; i < R.rows(); i++)
    {
        for (int j = 0; j < R.cols(); ++j)
        {
            iR(i, j) = round(R(i, j));
        }
    }
    cout << "\nR的整数形式=iR=\n" << iR << endl;
    auto ib = b.cast<int>();
    auto b1 = iR * ib;
    cout << "根据三个等式，求出a,a1,b,b1,c,c1. 其中b，b1已知。三个等式为：\n";
    cout << " iR*a=a1\n";
    cout << " iR*b=b1\n";
    cout << " iR*c=c1\n";
    cout << " 计算出的结果，如果a=a1,c=c1则为对称晶界，否则为非对称晶界。";
    cout << "下面对a[0],a[1],a[2], c[0],c[1],c[2], 分别从[-9,9]尝试，看看上面三个等式是否成立：" << endl;
    cout << ib << endl;
    Vector3i a;
    Vector3i c;
    int length =10;

	std::map<std::string, Result> resultList;
    for (int a_1 = -length; a_1 <= length; a_1++)
    {
        a[0] = a_1;
        for (int a_2 = -length; a_2 <= length; a_2++)
        {
            a[1] = a_2;
            for (int a_3 = -length; a_3 <= length; a_3++)
            {
                a[2] = a_3;

                for (int c_1 = -length; c_1 <= length; c_1++)
                {
                    c[0] = c_1;
                    for (int c_2 = -length; c_2 <= length; c_2++)
                    {
                        c[1] = c_2;
                        for (int c_3 = -length; c_3 <= length; c_3++)
                        {
                            c[2] = c_3;

                            if (a.transpose() * ib == 0
                                && a.transpose() * c == 0
                                && ib.transpose() * c == 0
                                && a.norm() != 0
                                && c.norm() != 0
                                )
                            {
                                Vector3i a1 = iR * a;
                                auto c1 = iR * c;
                                //cout << "\n找到了一组a和c向量：\n";
								Result result;
                                result.a = a;
                                result.b = ib;
                                result.c = c;
                                result.a1 = a1;
                                result.b1 = b1;
                                result.c1 = c1;

                                gcd(result);
                                result.create_id();

                                if (my_abs_ok(result.a1) && my_abs_ok(result.b1) && my_abs_ok(result.c1))
								{
									resultList.insert(std::make_pair(result.id, result));
								}

                            }
                        }
                    }
                }
            }
        }
    }
    if (resultList.empty())
    {
        cout << "@@@@@@@@@@@ resultList is empty(). not found!\n";
    }
    cout << "*************end****************" << endl;
    return resultList;
}
bool is_ok(Result const& item)
{
    //判断输出结果
    if (item.a.cwiseAbs() == item.a1.cwiseAbs()
        && item.c.cwiseAbs() == item.c1.cwiseAbs())
    {
        return true;
    }
    else
    {
        return false;
    }
}
void use_x0y0(Eigen::Vector2i &x0y0, vector3 &b, int N, int &segema)
{
    //单位向量[d, e, f]
    ofstream out_txt1, out_txt2;
    auto d_e_f = b.cast<float>().normalized();
    cout << "向量b的单位向量[d,e,f]=" << d_e_f.transpose() << endl;
    auto seata_degree = get_seta_degree(x0y0, N);
    std::cout << "第1步，求出的seta值(弧度)为："
        << seata_degree << ", 角度为：" << arc_radian(seata_degree) << std::endl;
    seata_degree = arc_radian(seata_degree);
    Vector3d bd = b.cast<double>();
    auto result = Processing(bd, seata_degree, segema);
    char* path1;
    char* path2;
    string ns1 = "result\\矩阵\\"+ to_string(segema)+".txt";
    string ns2 = "result\\5个输入\\"+to_string(segema)+".txt";
    path1 = (char*)ns1.c_str();
    path2 = (char*)ns2.c_str();
   
    int count = 0;
	for (auto& item : result)
	{
        if (is_ok(item.second))
        {
            ++count;
        }
	}
    if (count > 0)
    {
        out_txt1.open(path1, ios::out | ios::trunc);
        out_txt2.open(path2, ios::out | ios::trunc);
        out_txt2 << fixed;
        out_txt2 << fixed;
        out_txt1 << endl;
        out_txt2 << endl;
        for (auto& item : result)
        {
            if (is_ok(item.second))
            {
                out_txt1 << item.second;
                out_txt2 << "segema:" << segema << "  " << "seata:" << seata_degree << " " << "| a | :" << sqrt(get_N(item.second.a)) << " " << " | b | : " << sqrt(get_N(item.second.b)) << "  " << "| c | : " << sqrt(get_N(item.second.c)) << endl;
            }
        }
        out_txt1 << endl;
        out_txt1 << "result size = " << count << endl;
        out_txt1.close();
        out_txt2.close();
    }
}

int main()
{

    //auto r = atan(1);
    //cout << r * 4 / pai << endl;
    auto b = get_b();
    auto N = get_N(b);
    cout << "N=" << N << endl;
    for (int k=1; k <=60; ++k)
    {
        int segema =k;
        cout <<"segama="<<k<<"时:" << endl;
        auto x0y0list = get_x0_y0(segema, N);
        for (size_t i = 0; i < x0y0list.size(); i++)
        {
            cout << "\n现在对找到的每一对x0y0，进行后续计算" << i + 1 << "/" << x0y0list.size() << "：" << endl;
            //判断有没有合适的x0,y0存在
            if (!is_x0_y0_valid(x0y0list[i]))
            {
                cout << "x0,y0不合适:" << x0y0list[i] << endl;
                continue;
            }
            use_x0y0(x0y0list[i], b, N, segema);
            cout << endl << endl;
        }
    }
    cin.get();
    string line;
    getline(cin, line);
	pl("请按q键盘退出！");


}