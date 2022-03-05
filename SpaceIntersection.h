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


class SpaceIntersection
{
public:
	static void read_file_for_space_intersection(string file1, string file2, vector<double>& intri1, vector<double>& extri1, vector<Ppoint>& P1,
		vector<double>& intri2, vector<double>& extri2, vector<Ppoint>& P2, vector<int>& ids);
	static void calculate_rotation_matrix(Mat_<double>& R, const double fai, const double omega, const double kappa);
	static Mat coordinate_change(Mat_<double>& R, double x, double y, double f);
	static Gpoint  N1N2_intersection(vector<double> extri1, vector<double> extri2, double  N1, double N2, Mat_<double>& C, Mat_<double>& C_);
	static void space_intersection_pointfactor(string file1, string file2, string outfile);
};
