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

//#define x0 0
//#define y0 0
//#define f (153.24/1000.0)
//#define m 50000
//#define limit 0.01

class SpaceResection
{
public:
	static void read_file_for_space_resection(string file1, vector<Gpoint>& G, vector<Ppoint>& P);
	static void init(double& Xs, double& Ys, double& Zs, double& fai, double& omega, double& kappa, const vector<Gpoint>& G, double m, double f);
	static void calculate_rotation_matrix(Mat_<double>& R, const double fai, const double omega, const double kappa);
	static void calculate_A_matrix(Mat_<double>& A, Mat_<double>& R, vector<Ppoint>& P, double Z, int i, double fai, double omega, double kappa, double f);
	static void calculate_L_matrix(Mat_<double>& L, const vector<Ppoint>& P, const vector<Ppoint>& A, int i);
	static void correction(double& Xs, double& Ys, double& Zs, double& fai, double& omega, double& kappa, Mat_<double>& X);
	static bool ifadjustmentend(Mat_<double>& X, double limit);
	static void space_resection(string file1, string file2);
};
