#include "SpaceResection.h"

using namespace std;
using namespace cv;


void SpaceResection::read_file_for_space_resection(string file1, vector<Gpoint>& G, vector<Ppoint>& P)
{
    ifstream infile;
    infile.open(file1, ios::in);
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
        string split[6];
        while (str >> split[0] >> split[1] >> split[2] >> split[3] >> split[4] >> split[5])
        {
            int id = atoi(split[0].c_str());
            double ix = atof(split[1].c_str()) / 1000.0;
            double iy = atof(split[2].c_str()) / 1000.0;
            double gx = atof(split[3].c_str());
            double gy = atof(split[4].c_str());
            double gz = atof(split[5].c_str());
            P.push_back(Ppoint(ix, iy));
            G.push_back(Gpoint(gx, gy, gz));
        }
    }
    cout << P.size() << " Points" << endl;
}

void SpaceResection::init(double& Xs, double& Ys, double& Zs, double& fai, double& omega, double& kappa, const vector<Gpoint>& G, double m, double f)
{
    double Z = 0;
    int point_num = G.size();
    for (vector<Gpoint>::const_iterator it = G.cbegin(); it != G.cend(); it++)
    {
        Xs += (*it).x;
        Ys += (*it).y;
        Z += (*it).z;
    }
    Xs /= point_num;
    Ys /= point_num;
    Zs /= m * f + Z / point_num;
    fai = 0;
    omega = 0;
    kappa = 0;
}

void SpaceResection::calculate_rotation_matrix(Mat_<double>& R, const double fai, const double omega, const double kappa)
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

void SpaceResection::calculate_A_matrix(Mat_<double>& A, Mat_<double>& R, vector<Ppoint>& P, double Z, int i, double fai, double omega, double kappa, double f)
{
    A.at<double>(i * 2, 0) = (R.at<double>(0, 0) * f + R.at<double>(0, 2) * P[i].x) / Z;//a11
    A.at<double>(i * 2, 1) = (R.at<double>(1, 0) * f + R.at<double>(1, 2) * P[i].x) / Z;//a12
    A.at<double>(i * 2, 2) = (R.at<double>(2, 0) * f + R.at<double>(2, 2) * P[i].x) / Z;//a13
    A.at<double>(i * 2, 3) = P[i].y * sin(omega) - (P[i].x * (P[i].x * cos(kappa) - P[i].y * sin(kappa)) / f + f * cos(kappa)) * cos(omega);//a14
    A.at<double>(i * 2, 4) = -f * sin(kappa) - P[i].x * (P[i].x * sin(kappa) + P[i].y * cos(kappa)) / f;//a15
    A.at<double>(i * 2, 5) = P[i].y;//a16
    A.at<double>(i * 2 + 1, 0) = (R.at<double>(0, 1) * f + R.at<double>(0, 2) * P[i].y) / Z;//a21
    A.at<double>(i * 2 + 1, 1) = (R.at<double>(1, 1) * f + R.at<double>(1, 2) * P[i].y) / Z;//a22
    A.at<double>(i * 2 + 1, 2) = (R.at<double>(2, 1) * f + R.at<double>(2, 2) * P[i].y) / Z;//a23
    A.at<double>(i * 2 + 1, 3) = -P[i].x * sin(omega) - (P[i].y * (P[i].x * cos(kappa) - P[i].y * sin(kappa)) / f - f * sin(kappa)) * cos(omega);//a24
    A.at<double>(i * 2 + 1, 4) = -f * cos(kappa) - (P[i].y * (P[i].x * sin(kappa) + P[i].y * cos(kappa))) / f;//a25
    A.at<double>(i * 2 + 1, 5) = -P[i].x;//a26
}

void SpaceResection::calculate_L_matrix(Mat_<double>& L, const vector<Ppoint>& P, const vector<Ppoint>& A, int i)
{
    L.at<double>(i * 2, 0) = P[i].x - A[i].x;
    L.at<double>(i * 2 + 1, 0) = P[i].y - A[i].y;
}

void SpaceResection::correction(double& Xs, double& Ys, double& Zs, double& fai, double& omega, double& kappa, Mat_<double>& X)
{
    Xs += X.at<double>(0, 0);
    Ys += X.at<double>(1, 0);
    Zs += X.at<double>(2, 0);
    fai += X.at<double>(3, 0);
    omega += X.at<double>(4, 0);
    kappa += X.at<double>(5, 0);
}

bool SpaceResection::ifadjustmentend(Mat_<double>& X, double limit)
{
    if (fabs(X.at<double>(3, 0)) < limit && fabs(X.at<double>(4, 0)) < limit && fabs(X.at<double>(5, 0)) < limit)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void SpaceResection::space_resection(string file1, string file2)
{
    ////Define, Initiate and points input
    double f = (153.24 / 1000.0);
    double m = 50000;
    double limit = 0.01;
    int iteration = 0;
    double Xs = 0, Ys = 0, Zs = 0;
    double fai = 0, omega = 0, kappa = 0;
    Mat_<double> R = Mat::zeros(3, 3, CV_32F);
    Mat_<double> X = Mat::zeros(6, 1, CV_32F);
    vector<Ppoint> P;
    vector<Gpoint> G;
    SpaceResection::read_file_for_space_resection(file1, G, P);
    int point_num = P.size();
    vector<Ppoint> approx(point_num);
    Mat_<double> A = Mat::zeros(point_num * 2, 6, CV_32F);
    Mat_<double> L = Mat::zeros(point_num * 2, 1, CV_32F);
    //Get initial parameters
    init(Xs, Ys, Zs, fai, omega, kappa, G, m, f);
    
    ////Adjustment process with iterations
    do
    {
        //Calculate rotation matrix
        calculate_rotation_matrix(R, fai, omega, kappa);
        //Calculate approximation of x,y of each point
        //Form norm equation of each point
        for (int i = 0; i < point_num; i++)
        {
            //Approximation
            approx[i].x = -f * (R.at<double>(0, 0) * (G[i].x - Xs) + R.at<double>(1, 0) * (G[i].y - Ys) + R.at<double>(2, 0) * (G[i].z - Zs))
                / (R.at<double>(0, 2) * (G[i].x - Xs) + R.at<double>(1, 2) * (G[i].y - Ys) + R.at<double>(2, 2) * (G[i].z - Zs));
            approx[i].y = -f * (R.at<double>(0, 1) * (G[i].x - Xs) + R.at<double>(1, 1) * (G[i].y - Ys) + R.at<double>(2, 1) * (G[i].z - Zs))
                / (R.at<double>(0, 2) * (G[i].x - Xs) + R.at<double>(1, 2) * (G[i].y - Ys) + R.at<double>(2, 2) * (G[i].z - Zs));
            //Norm equation
            double Z = R.at<double>(0, 2) * (G[i].x - Xs) + R.at<double>(1, 2) * (G[i].y - Ys) + R.at<double>(2, 2) * (G[i].z - Zs);
            calculate_A_matrix(A, R, P, Z, i, fai, omega, kappa, f);
            calculate_L_matrix(L, P, approx, i);
        }
        //Calculate correction according to norm equation
        //X = inv(A' * A) * A' * L
        X = (A.t() * A).inv() * A.t() * L;
        correction(Xs, Ys, Zs, fai, omega, kappa, X);
        iteration += 1;
    } while (!ifadjustmentend(X, limit));
    Mat_<double> V = A * X - L;
    Mat_<double> V_ = V.t() * V;
    double accuracy = sqrt(V_.at<double>(0, 0) / (point_num * 2 - 6));

    ////Output final result
    cout << "--------------------------------------------" << endl;
    cout << "Space Resection Result" << endl;
    cout << "Iteration: " << iteration << endl;
    iteration = 0;
    cout << "Xs = " << fixed << setprecision(2) << Xs << endl;
    cout << "Ys = " << fixed << setprecision(2) << Ys << endl;
    cout << "Zs = " << fixed << setprecision(2) << Zs << endl;
    cout << "fai = " << fixed << setprecision(5) << fai << endl;
    cout << "omega = " << fixed << setprecision(5) << omega << endl;
    cout << "kappa = " << fixed << setprecision(5) << kappa << endl;
    cout << "Rotation Matrix:" << endl;
    cout << fixed << setprecision(5) << R << endl;
    cout << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
    cout << endl;

    ofstream outfile;
    outfile.open(file2, ios::out);
    outfile << "Xs = " << fixed << setprecision(2) << Xs << endl;
    outfile << "Ys = " << fixed << setprecision(2) << Ys << endl;
    outfile << "Zs = " << fixed << setprecision(2) << Zs << endl;
    outfile << "fai = " << fixed << setprecision(5) << fai << endl;
    outfile << "omega = " << fixed << setprecision(5) << omega << endl;
    outfile << "kappa = " << fixed << setprecision(5) << kappa << endl;
    outfile << "Rotation Matrix:" << endl;
    outfile << fixed << setprecision(5) << R << endl;
    outfile << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
    outfile << endl;
}