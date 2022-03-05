#include "InteriorOrientation.h"

using namespace std;
using namespace cv;


void InteriorOrientation::read_file_for_interior_orientation(string file, vector<Ppoint>& Pm, vector<Ppoint>& Pi)
{
    ifstream infile;
    infile.open(file, ios::in);
    if (!infile)
    {
        cout << "File doesn't exisit" << endl;
        return;
    }
    string a = "";
    while (!infile.eof())
    {
        getline(infile, a, '\n');
        //cout << a << endl;
        istringstream str(a);
        string split[4];
        while (str >> split[0] >> split[1] >> split[2] >> split[3])
        {
            double mx = atof(split[0].c_str());
            double my = atof(split[1].c_str());
            double ix = atof(split[2].c_str());
            double iy = atof(split[3].c_str());
            Pm.push_back(Ppoint(mx, my));
            Pi.push_back(Ppoint(ix, iy));
        }
    }
    cout << Pm.size() << " Points" << endl;
}

void InteriorOrientation::calculate_A_matrix(Mat_<double>& A, vector<Ppoint>& Pi, double pixel)
{
    int point_num = Pi.size();
    for (int i = 0; i < point_num; i++)
    {
        A.at<double>(i * 2, 0) = 1;
        A.at<double>(i * 2, 1) = Pi[i].x * pixel;
        A.at<double>(i * 2, 2) = -Pi[i].y * pixel;
        A.at<double>(i * 2, 3) = 0;
        A.at<double>(i * 2, 4) = 0;
        A.at<double>(i * 2, 5) = 0;

        A.at<double>(i * 2 + 1, 0) = 0;
        A.at<double>(i * 2 + 1, 1) = 0;
        A.at<double>(i * 2 + 1, 2) = 0;
        A.at<double>(i * 2 + 1, 3) = 1;
        A.at<double>(i * 2 + 1, 4) = Pi[i].x * pixel;
        A.at<double>(i * 2 + 1, 5) = -Pi[i].y * pixel;
    }
}

void InteriorOrientation::calculate_L_matrix(Mat_<double>& L, vector<Ppoint>& Pm)
{
    int point_num = Pm.size();
    for (int i = 0; i < point_num; i++)
    {
        L.at<double>(i * 2, 0) = Pm[i].x;
        L.at<double>(i * 2 + 1, 0) = Pm[i].y;
    }
}

void InteriorOrientation::interior_orientation(string file, string outputfile, vector<double>& parameter)
{
    ////Define and Initiate parameters
    vector<Ppoint> Pm;
    vector<Ppoint> Pi;
    InteriorOrientation::read_file_for_interior_orientation(file, Pm, Pi);
    double pixel = 0.021;
    double width = 11129;
    double height = 11263;
    double x0 = 0.011;
    double y0 = 0.002;
    double x_offset = (width * pixel) / 2;
    double y_offset = -(height * pixel) / 2;
    int point_num = Pm.size();
    vector<Ppoint> result;
    Mat_<double> A = Mat::zeros(point_num * 2, 6, CV_32F);
    Mat_<double> L = Mat::zeros(point_num * 2, 1, CV_32F);
    Mat_<double> Para = Mat::zeros(6, 1, CV_32F);
    
    ////Form norm equation to calculate parameters
    calculate_A_matrix(A, Pi, pixel);
    //cout << A << endl;
    calculate_L_matrix(L, Pm);
    //cout << L << endl;
    //X(Para) = inv(A' * A) * A' * L
    Para = (A.t() * A).inv() * A.t() * L;
    Mat_<double> V = A * Para - L;
    Mat_<double> V_ = V.t() * V;
    double accuracy = sqrt(V_.at<double>(0, 0) / (point_num * 2 - 6));

    double m0 = Para.at<double>(0, 0);
    double m1 = Para.at<double>(1, 0);
    double m2 = Para.at<double>(2, 0);
    double n0 = Para.at<double>(3, 0);
    double n1 = Para.at<double>(4, 0);
    double n2 = Para.at<double>(5, 0);

    ////Output result
    cout << "--------------------------------------------" << endl;
    cout << "Interior Orientation Result" << endl;
    cout << "m0 = " << fixed << setprecision(5) << m0 << endl;
    cout << "m1 = " << fixed << setprecision(5) << m1 << endl;
    cout << "m2 = " << fixed << setprecision(5) << m2 << endl;
    cout << "n0 = " << fixed << setprecision(5) << n0 << endl;
    cout << "n1 = " << fixed << setprecision(5) << n1 << endl;
    cout << "n2 = " << fixed << setprecision(5) << n2 << endl;
    cout << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
    ofstream outfile;
    outfile.open(outputfile, ios::out);
    outfile << "m0 = " << fixed << setprecision(5) << m0 << endl;
    outfile << "m1 = " << fixed << setprecision(5) << m1 << endl;
    outfile << "m2 = " << fixed << setprecision(5) << m2 << endl;
    outfile << "n0 = " << fixed << setprecision(5) << n0 << endl;
    outfile << "n1 = " << fixed << setprecision(5) << n1 << endl;
    outfile << "n2 = " << fixed << setprecision(5) << n2 << endl;
    outfile << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
    parameter.push_back(m0);
    parameter.push_back(m1);
    parameter.push_back(m2);
    parameter.push_back(n0);
    parameter.push_back(n1);
    parameter.push_back(n2);
    parameter.push_back(x0);
    parameter.push_back(y0);
}