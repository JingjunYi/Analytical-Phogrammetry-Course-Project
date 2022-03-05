#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <string>
#include "Gpoint.h"
#include "Ppoint.h"
#include "SpaceResection.h"
#include "SpaceIntersection.h"
#include "InteriorOrientation.h";
#include "RelativeOrientation.h"
#include "AbsoluteOrientation.h"

using namespace std;
using namespace cv;

int main()
{
    string work_pattern = "";
    cout << "Choose work_patter: A(Space Resection), B(Space Intersection), C(Interior Orientation)" << endl;
    cout << "D(Relative Orientation), E(Absolute Orientation)" << endl;
    cin >> work_pattern;
    //Task1 Space Resection
    if (work_pattern == "A")
    {
        string file1 = "..\\space_resection.txt";
        string output_path1 = "..\\space_resection_result.txt";
        SpaceResection::space_resection(file1, output_path1);
    }

    //Task2 Space Intersection
    else if (work_pattern == "B")
    {
        string file2_1 = "..\\spcae_intersection319.txt";
        string file2_2 = "..\\spcae_intersection320.txt";
        string output_path2 = "..\\space_intersection_result.txt";
        /*string file2_1 = "..\\spcae_intersection_example1.txt";
        string file2_2 = "..\\spcae_intersection_example2.txt";
        string output_path2 = "..\\space_intersection_result_example.txt";*/
        SpaceIntersection::space_intersection_pointfactor(file2_1, file2_2, output_path2);
    }

    //Task3 Interior Orientation
    else if (work_pattern == "C")
    {
        string file3 = "..\\interior_orientation.txt";
        string output_path3 = "..\\interior_orientation_result.txt";
        vector<double> parameter;
        InteriorOrientation::interior_orientation(file3, output_path3, parameter);
        
        cout << "--------------------------------------------" << endl;
        while (1)
        {
            cout << "Input location(row and column) of pixel to get principal coordinates:(format - i j)" << endl;
            double i, j;
            cin >> i >> j;
            return 0;
            if (i != int(i) || j != int(j))
            {
                cout << "Wrong location, has to be int" << endl;
                break;
            }
            double m0 = parameter[0];
            double m1 = parameter[1];
            double m2 = parameter[2];
            double n0 = parameter[3];
            double n1 = parameter[4];
            double n2 = parameter[5];
            double xbias = parameter[6];
            double ybias = parameter[7];
            double x = m0 + m1 * i - m2 * j - xbias;
            double y = n0 + n1 * i - n2 * j - ybias;
            cout << "Principal coordinates of this point: " << x << " " << y << endl;
        }
    }

    //Task4 Relative Orientation
    else if (work_pattern == "D")
    {
        string file4 = "..\\relative_orientation.txt";
        string output_path4 = "..\\relative_orientation_result.txt";
        RelativeOrientation::relative_orientation(file4, output_path4);
    }

    //Task Absolute Orientation
    else if (work_pattern == "E")
    {
        string file5 = "..\\absolute_orientation.txt";
        string output_path5 = "..\\absolute_orientation_result.txt";
        AbsoluteOrientation::absolute_orientation(file5, output_path5);
    }
    return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件


