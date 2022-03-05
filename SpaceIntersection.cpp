#include "SpaceIntersection.h"

using namespace std;
using namespace cv;


void SpaceIntersection::read_file_for_space_intersection(string file1, string file2, vector<double>& intri1, vector<double>& extri1, vector<Ppoint>& P1,
    vector<double>& intri2, vector<double>& extri2, vector<Ppoint>& P2, vector<int>& ids)
{
    ifstream infile1;
    infile1.open(file1, ios::in);
    if (!infile1)
    {
        cout << "File1 doesn't exisit" << endl;
        return;
    }
    string a = "";
    int i = 0;
    cout << "--------------------------------------------" << endl;
    cout << "Input parameters and image coordinates" << endl;
    while (!infile1.eof())
    {
        getline(infile1, a, '\n');
        i += 1;
        cout << a << endl;
        istringstream str(a);
        string intrinsics[4];
        string extrinsics[6];
        string split[3];
        if (i == 1)
        {
            str >> intrinsics[0] >> intrinsics[1] >> intrinsics[2] >> intrinsics[3];
            double f = atof(intrinsics[0].c_str()) / 1000.0;
            double x0 = atof(intrinsics[1].c_str()) / 1000.0;
            double y0 = atof(intrinsics[2].c_str()) / 1000.0;
            double m = atof(intrinsics[3].c_str());
            intri1.push_back(f);
            intri1.push_back(x0);
            intri1.push_back(y0);
            intri1.push_back(m);
        }
        else if (i == 2)
        {
            str >> extrinsics[0] >> extrinsics[1] >> extrinsics[2] >> extrinsics[3] >> extrinsics[4] >> extrinsics[5];
            //double Xs = atof(extrinsics[0].c_str());
            //double Ys = atof(extrinsics[1].c_str());
            double Ys = atof(extrinsics[0].c_str());
            double Xs = atof(extrinsics[1].c_str());
            double Zs = atof(extrinsics[2].c_str());
            double fai = atof(extrinsics[3].c_str());
            double omega = atof(extrinsics[4].c_str());
            //double omega = atof(extrinsics[3].c_str());
            //double fai = atof(extrinsics[4].c_str());
            double kappa = atof(extrinsics[5].c_str());
            fai = fai * CV_PI / 180;
            omega = omega * CV_PI / 180;
            kappa = kappa * CV_PI / 180;
            extri1.push_back(Xs);
            extri1.push_back(Ys);
            extri1.push_back(Zs);
            extri1.push_back(fai);
            extri1.push_back(omega);
            extri1.push_back(kappa);
        }
        else
        {
            while (str >> split[0] >> split[1] >> split[2])
            {
                int id = atoi(split[0].c_str());
                double ix = atof(split[1].c_str()) / 1000.0;
                double iy = atof(split[2].c_str()) / 1000.0;
                P1.push_back(Ppoint(ix, iy));
                ids.push_back(id);
            }
        }
    }
    cout << P1.size() << " Points" << endl;

    ifstream infile2;
    infile2.open(file2, ios::in);
    if (!infile2)
    {
        cout << "File2 doesn't exisit" << endl;
        return;
    }
    string a_ = "";
    int i_ = 0;
    while (!infile2.eof())
    {
        getline(infile2, a_, '\n');
        i_ += 1;
        cout << a_ << endl;
        istringstream str_(a_);
        string intrinsics_[4];
        string extrinsics_[6];
        string split_[3];
        if (i_ == 1)
        {
            str_ >> intrinsics_[0] >> intrinsics_[1] >> intrinsics_[2] >> intrinsics_[3];
            double f_ = atof(intrinsics_[0].c_str()) / 1000.0;
            double x0_ = atof(intrinsics_[1].c_str()) / 1000.0;
            double y0_ = atof(intrinsics_[2].c_str()) / 1000.0;
            double m_ = atof(intrinsics_[3].c_str());
            intri2.push_back(f_);
            intri2.push_back(x0_);
            intri2.push_back(y0_);
            intri2.push_back(m_);
        }
        else if (i_ == 2)
        {
            str_ >> extrinsics_[0] >> extrinsics_[1] >> extrinsics_[2] >> extrinsics_[3] >> extrinsics_[4] >> extrinsics_[5];
            //double Xs_ = atof(extrinsics_[0].c_str());
            //double Ys_ = atof(extrinsics_[1].c_str());
            double Ys_ = atof(extrinsics_[0].c_str());
            double Xs_ = atof(extrinsics_[1].c_str());
            double Zs_ = atof(extrinsics_[2].c_str());
            double fai_ = atof(extrinsics_[3].c_str());
            double omega_ = atof(extrinsics_[4].c_str());
            //double omega_ = atof(extrinsics_[3].c_str());
            //double fai_ = atof(extrinsics_[4].c_str());
            double kappa_ = atof(extrinsics_[5].c_str());
            fai_ = fai_ * CV_PI / 180;
            omega_ = omega_ * CV_PI / 180;
            kappa_ = kappa_ * CV_PI / 180;
            extri2.push_back(Xs_);
            extri2.push_back(Ys_);
            extri2.push_back(Zs_);
            extri2.push_back(fai_);
            extri2.push_back(omega_);
            extri2.push_back(kappa_);
        }
        else
        {
            while (str_ >> split_[0] >> split_[1] >> split_[2])
            {
                int id_ = atoi(split_[0].c_str());
                double ix_ = atof(split_[1].c_str()) / 1000.0;
                double iy_ = atof(split_[2].c_str()) / 1000.0;
                P2.push_back(Ppoint(ix_, iy_));
            }
        }
    }
    cout << P2.size() << " Points" << endl;
}

void SpaceIntersection::calculate_rotation_matrix(Mat_<double>& R, const double fai, const double omega, const double kappa)
{
    R.at<double>(0, 0) = cos(fai) * cos(kappa) - sin(fai) * sin(omega) * sin(kappa);
    R.at<double>(0, 1) = -cos(fai) * sin(kappa) - sin(fai) * sin(omega) * cos(kappa);
    R.at<double>(0, 2) = -sin(fai) * cos(omega);
    R.at<double>(1, 0) = cos(omega) * sin(kappa);
    R.at<double>(1, 1) = cos(omega) * cos(kappa);
    R.at<double>(1, 2) = -sin(omega);
    R.at<double>(2, 0) = sin(fai) * cos(kappa) + cos(fai) * sin(omega) * sin(kappa);
    R.at<double>(2, 1) = -sin(fai) * sin(kappa) + cos(fai) * sin(omega) * cos(kappa);
    R.at<double>(2, 2) = cos(fai) * cos(omega);
}

Mat SpaceIntersection::coordinate_change(Mat_<double>& R, double x, double y, double f)
{
    Mat_<double> P(3, 1);
    P.at<double>(0, 0) = x;
    P.at<double>(1, 0) = y;
    P.at<double>(2, 0) = -f;
    Mat_<double> Aux = R * P;
    return Aux;
}

Gpoint  SpaceIntersection::N1N2_intersection(vector<double> extri1, vector<double> extri2, double  N1, double N2, Mat_<double>& C, Mat_<double>& C_)
{
    double tmp[3] = {};
    double tmp2[3] = {};
    double final_point[3] = {};
    for (int i = 0; i < 3; i++)
    {
        tmp[i] = extri1[i] + N1 * C.at<double>(i, 0);
        tmp2[i] = extri2[i] + N2 * C_.at<double>(i, 0);
        //For X,Y,Z tak average as Y should be 1/2(Y1+Y2)
        final_point[i] = (tmp[i] + tmp2[i]) / 2;
    }
    return Gpoint(final_point[0], final_point[1], final_point[2]);
}

void SpaceIntersection::space_intersection_pointfactor(string file1, string file2, string outputfile)
{
    vector<double> intri1, intri2;
    vector<double> extri1, extri2;
    vector<Ppoint> P1, P2;
    vector<Gpoint> G;
    vector<int> ids;
    SpaceIntersection::read_file_for_space_intersection(file1, file2, intri1, extri1, P1, intri2, extri2, P2, ids);
    //Define and Initiate parameters
    double f = intri1[0];
    double x0 = intri1[1];
    double y0 = intri1[2];
    double m = intri1[3];
    double Xs = extri1[0];
    double Ys = extri1[1];
    double Zs = extri1[2];
    double fai = extri1[3];
    double omega = extri1[4];
    double kappa = extri1[5];
    double f_ = intri2[0];
    double x0_ = intri2[1];
    double y0_ = intri2[2];
    double m_ = intri2[3];
    double Xs_ = extri2[0];
    double Ys_ = extri2[1];
    double Zs_ = extri2[2];
    double fai_ = extri2[3];
    double omega_ = extri2[4];
    double kappa_ = extri2[5];
    //Rotation Matrixs
    Mat_<double> R = Mat::zeros(3, 3, CV_32F);
    Mat_<double> R_ = Mat::zeros(3, 3, CV_32F);
    //Baseline
    double Bx = Xs_ - Xs;
    double By = Ys_ - Ys;
    double Bz = Zs_ - Zs;
    //Image auxiliary coordinates
    Mat_<double> C = Mat::zeros(3, 1, CV_32F);
    Mat_<double> C_ = Mat::zeros(3, 1, CV_32F);

    ////Calculate
    //Rotation matrixs
    calculate_rotation_matrix(R, fai, omega, kappa);
    calculate_rotation_matrix(R_, fai_, omega_, kappa_);
    int point_num = P1.size();
    for (int i = 0; i < point_num; i++)
    {
        //Transfer to auxiliary coordinates
        C = coordinate_change(R, P1[i].x, P1[i].y, f);
        C_ = coordinate_change(R_, P2[i].x, P2[i].y, f_);
        //Point factor
        double N1 = (Bx * C_.at<double>(2, 0) - Bz * C_.at<double>(0, 0)) / (C.at<double>(0, 0) * C_.at<double>(2, 0) - C.at<double>(2, 0) * C_.at<double>(0, 0));
        double N2 = (Bx * C.at<double>(2, 0) - Bz * C.at<double>(0, 0)) / (C.at<double>(0, 0) * C_.at<double>(2, 0) - C.at<double>(2, 0) * C_.at<double>(0, 0));
        //Intersection based on N1, N2
        Gpoint A = N1N2_intersection(extri1, extri2, N1, N2, C, C_);
        G.push_back(A);
    }

    ////Output result
    cout << "--------------------------------------------" << endl;
    cout << "Space Intersection result of point factor N1 N2 method" << endl;
    for (int i = 0; i < G.size(); i++)
    {
        cout << fixed << setprecision(5) << ids[i] << " " << G[i].y << " " << G[i].x << " " << G[i].z << endl;
    }
    cout << endl;
    ofstream outfile;
    outfile.open(outputfile, ios::out);
    for (int i = 0; i < G.size(); i++)
    {
        outfile << fixed << setprecision(5) << ids[i] << " " << G[i].y << " " << G[i].x << " " << G[i].z << endl;
    }
}

//void space_intersection_bundle(string file1, string file2)
//{
//    vector<double> intri1, intri2;
//    vector<double> extri1, extri2;
//    vector<Ppoint> P1, P2;
//    read_file_for_space_intersection(file1, file2, intri1, extri1, P1, intri2, extri2, P2);
//    //Define and Initiate parameters
//    double f = intri1[0];
//    double x0 = intri1[1];
//    double y0 = intri1[2];
//    double m = intri1[3];
//    double Xs = extri1[0];
//    double Ys = extri1[1];
//    double Zs = extri1[2];
//    double fai = extri1[3];
//    double omega = extri1[4];
//    double kappa = extri1[5];
//    double f_ = intri1[0];
//    double x0_ = intri1[1];
//    double y0_ = intri1[2];
//    double m_ = intri1[3];
//    double Xs_ = extri1[0];
//    double Ys_ = extri1[1];
//    double Zs_ = extri1[2];
//    double fai_ = extri1[3];
//    double omega_ = extri1[4];
//    double kappa_ = extri1[5];
//    //Rotation Matrixs
//    Mat_<double> R = Mat::zeros(3, 3, CV_32F);
//    Mat_<double> R_ = Mat::zeros(3, 3, CV_32F);
//    calculate_rotation_matrix(R, fai, omega, kappa);
//    calculate_rotation_matrix(R_, fai, omega, kappa);
//    //Baseline
//    double Bx = Xs_ - Xs;
//    double By = Ys_ - Ys;
//    double Bz = Zs_ - Zs;
//    //Image auxiliary coordinates
//    Mat_<double> C1 = Mat::zeros(3, 1, CV_32F);
//    Mat_<double> C2 = Mat::zeros(3, 1, CV_32F);
//
//    ////Calculate
//}