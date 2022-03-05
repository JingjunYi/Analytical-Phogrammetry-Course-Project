#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <string>
#include "Gpoint.h"
#include "Ppoint.h"

using namespace std;
using namespace cv;


class AbsoluteOrientation
{
public:
	static void read_file_for_absolute_orientation(string file, vector<string>& pname, vector<Gpoint>& Pmodel, vector<Gpoint>& Pspace);
	static void calculate_rotation_matrix(Mat_<double>& R, Mat_<double>& Para);
	static void calculate_L_matrix(int i, Mat_<double>& L, Mat_<double>& mtp,  Mat_<double>& mp, Mat_<double>& p0, Mat_<double>& Para);
	static void calculate_A_matrix(int i, Mat_<double>& A, vector<Gpoint>& Pmodel, Mat_<double>& Para);
	static void correct(Mat_<double>& Para, Mat_<double>& X);
	static void absolute_orientation(string file, string outputfile);
};

