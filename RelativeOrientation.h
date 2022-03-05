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


class RelativeOrientation
{
public:
	static void read_file_for_relative_orientation(string file, string img1, string img2, vector<Ppoint>& P1, vector<Ppoint>& P2, vector<double>& intri1, vector<double>& intri2, vector<int>& pname);
	static void calculate_relarotation_matrix(Mat_<double>& R2, Mat_<double>& Para);
	static void calculate_A_matrix(int i, Mat_<double>& A, Gpoint& P, double N1, double N2, double Bx);
	static void calculate_L_matrix(int i, Mat_<double>& L, double Q);
	static void correct(Mat_<double>& Para, Mat_<double>& X);
	static void relative_orientation(string file, string outputfile);
};

