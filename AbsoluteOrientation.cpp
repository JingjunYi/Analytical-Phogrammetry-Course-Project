#include "AbsoluteOrientation.h"

using namespace std;
using namespace cv;


void AbsoluteOrientation::read_file_for_absolute_orientation(string file, vector<string>& pname, vector<Gpoint>& Pmodel, vector<Gpoint>& Pspace)
{
    ifstream infile;
    infile.open(file, ios::in);
    if (!infile)
    {
        cout << "File doesn't exisit" << endl;
        return;
    }
    string a = "";
    int i = 0;
    while (!infile.eof())
    {
        getline(infile, a, '\n');
        i += 1;
        //cout << a << endl;
        istringstream str(a);
        string split[7];
        while (str >> split[0] >> split[1] >> split[2] >> split[3] >> split[4] >> split[5] >> split[6])
        {
            string name = split[0];
            double xmodel = atof(split[1].c_str());
            double ymodel = atof(split[2].c_str());
            double zmodel = atof(split[3].c_str());
            double xspace = atof(split[4].c_str());
            double yspace = atof(split[5].c_str());
            double zspace = atof(split[6].c_str());
            pname.push_back(name);
            Pmodel.push_back(Gpoint(xmodel, ymodel, zmodel));
            Pspace.push_back(Gpoint(xspace, yspace, zspace));
        }
    }
}

void AbsoluteOrientation::calculate_rotation_matrix(Mat_<double>& R, Mat_<double>& Para)
{
    double fai = Para.at<double>(3, 0);
    double omega = Para.at<double>(4, 0);
    double kappa = Para.at<double>(5, 0);
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

void AbsoluteOrientation::calculate_L_matrix(int i, Mat_<double>& L, Mat_<double>& mtp, Mat_<double>& mp, Mat_<double>& p0, Mat_<double>& Para)
{
    Mat_<double> jL = Mat::zeros(3, 1, CV_32F);
    double s = Para.at<double>(6, 0);
    jL = mtp - s * mp - p0;
    L(3 * i, 0) = jL.at<double>(0, 0);
    L(3 * i + 1, 0) = jL.at<double>(1, 0);
    L(3 * i + 2, 0) = jL.at<double>(2, 0);
}
void AbsoluteOrientation::calculate_A_matrix(int i, Mat_<double>& A, vector<Gpoint>& Pmodel, Mat_<double>& Para)
{
    double fai = Para.at<double>(3, 0);
    double omega = Para.at<double>(4, 0);
    double kappa = Para.at<double>(5, 0);
    double s = Para.at<double>(6, 0);
    Mat_<double> jA = Mat::zeros(3, 7, CV_32F);
    jA.at<double>(0, 0) = 1; jA.at<double>(0, 1) = 0; jA.at<double>(0, 2) = 0; jA.at<double>(0, 3) = Pmodel[i].x;
    jA.at<double>(0, 4) = -s * Pmodel[i].z; //fai
    jA.at<double>(0, 5) = -s * Pmodel[i].y * sin(fai); //omega
    jA.at<double>(0, 6) = -s * Pmodel[i].y * cos(fai) * cos(omega) - s * Pmodel[i].z * sin(omega); //kappa;
    jA.at<double>(1, 0) = 0; jA.at<double>(1, 1) = 1; jA.at<double>(1, 2) = 0; jA.at<double>(1, 3) = Pmodel[i].y;
    jA.at<double>(1, 4) = 0; //fai
    jA.at<double>(1, 5) = s * Pmodel[i].x * sin(fai) - s * Pmodel[i].z * cos(fai); //omega
    jA.at<double>(1, 6) = s * Pmodel[i].x * cos(fai) * cos(omega) + s * Pmodel[i].z * sin(fai) * cos(omega); //kappa
    jA.at<double>(2, 0) = 0; jA.at<double>(2, 1) = 0; jA.at<double>(2, 2) = 1; jA.at<double>(2, 3) = Pmodel[i].z;
    jA.at<double>(2, 4) = s * Pmodel[i].x; //fai
    jA.at<double>(2, 5) = s * Pmodel[i].y * cos(fai); //omega
    jA.at<double>(2, 6) = s * Pmodel[i].x * sin(omega) - s * Pmodel[i].y * sin(fai) * cos(omega); //kappa
    for (int j = 0; j < 7; j++)
    {
        A.at<double>(3 * i, j) = jA.at<double>(0, j);
        A.at<double>(3 * i+1, j) = jA.at<double>(1, j);
        A.at<double>(3 * i+2, j) = jA.at<double>(2, j);
    }
}

void AbsoluteOrientation::correct(Mat_<double>& Para, Mat_<double>& X)
{
    Para.at<double>(0, 0) += X.at<double>(0, 0); //x
    Para.at<double>(1, 0) += X.at<double>(1, 0); //y
    Para.at<double>(2, 0) += X.at<double>(2, 0); //z
    Para.at<double>(6, 0) += X.at<double>(3, 0); //s
    Para.at<double>(3, 0) += X.at<double>(4, 0); //fai
    Para.at<double>(4, 0) += X.at<double>(5, 0); //omega
    Para.at<double>(5, 0) += X.at<double>(6, 0); //kappa
}

void AbsoluteOrientation::absolute_orientation(string file, string outputfile)
{
    ////Define and Initiate parameters
    vector<string> pname;
    vector<Gpoint> Pmodel;
    vector<Gpoint> Pspace;
    AbsoluteOrientation::read_file_for_absolute_orientation(file, pname, Pmodel, Pspace);
    int point_num = Pmodel.size();
    Mat_<double> Para = Mat::zeros(7, 1, CV_32F); //x, y, z, fai, omega, kappa, s
    Para.at<double>(6, 0) = 1; // s = 1
    Mat_<double> A = Mat::zeros(3 * point_num, 7, CV_32F);
    Mat_<double> L = Mat::zeros(3 * point_num, 1, CV_32F);
    Mat_<double> X = Mat::zeros(7, 1, CV_32F);
    //Calculate Scale
    double d1 = sqrt(pow((Pmodel[0].x - Pmodel[1].x),2) + pow((Pmodel[0].y - Pmodel[1].y), 2) + pow((Pmodel[0].z - Pmodel[1].z), 2));
    double d2 = sqrt(pow((Pspace[0].x - Pspace[1].x), 2) + pow((Pspace[0].y - Pspace[1].y), 2) + pow((Pspace[0].z - Pspace[1].z), 2));
    double m = d2 / d1;
    //m = 1;
    for (int i = 0; i < point_num; i++)
    {
        Pspace[i].x /= m;
        Pspace[i].y /= m;
        Pspace[i].z /= m;
    }
    ////Calculate Barycentric coordinates
    //Calculate Center coordinates
    double  xmodel_center = 0;
    double  ymodel_center = 0;
    double  zmodel_center = 0;
    double  xspace_center = 0;
    double  yspace_center = 0;
    double  zspace_center = 0;
    for (int i = 0; i < point_num; i++)
    {
        xmodel_center += Pmodel[i].x;
        ymodel_center += Pmodel[i].y;
        zmodel_center += Pmodel[i].z;
        xspace_center += Pspace[i].x;
        yspace_center += Pspace[i].y;
        zspace_center += Pspace[i].z;
    }
    xmodel_center /= point_num;
    ymodel_center /= point_num;
    zmodel_center /= point_num;
    xspace_center /= point_num;
    yspace_center /= point_num;
    zspace_center /= point_num;
    Gpoint model_center(xmodel_center, ymodel_center, zmodel_center);
    Gpoint space_center(xspace_center, yspace_center, zspace_center);
    //Barycenterization
    for (int i = 0; i < point_num; i++)
    {
        Pmodel[i].x -= model_center.x;
        Pmodel[i].y -= model_center.y;
        Pmodel[i].z -= model_center.z;
        Pspace[i].x -= space_center.x;
        Pspace[i].y -= space_center.y;
        Pspace[i].z -= space_center.z;
    }

    ////Calculate with Iteration
    int iteration = 0;
    while (true)
    {
        iteration += 1;
        cout << "Iteration: " << iteration << endl;
        Mat_<double> R = Mat::zeros(3, 3, CV_32F);
        AbsoluteOrientation::calculate_rotation_matrix(R, Para);
        //Calculate coefficient of each point
        for (int i = 0; i < point_num; i++)
        {
            Mat_<double> mtp = Mat::zeros(3, 1, CV_32F);
            Mat_<double> mp = Mat::zeros(3, 1, CV_32F);
            Mat_<double> t = Mat::zeros(3, 1, CV_32F);
            Mat_<double> p0 = Mat::zeros(3, 1, CV_32F);
            Gpoint P_(0, 0, 0);
            mtp.at<double>(0, 0) = Pspace[i].x;
            mtp.at<double>(1, 0) = Pspace[i].y;
            mtp.at<double>(2, 0) = Pspace[i].z;
            t.at<double>(0, 0) = Pmodel[i].x;
            t.at<double>(1, 0) = Pmodel[i].y;
            t.at<double>(2, 0) = Pmodel[i].z;
            mp = R * t;
            P_.x = mp.at<double>(0, 0);
            P_.y = mp.at<double>(1, 0);
            P_.z = mp.at<double>(2, 0);
            p0.at<double>(0, 0) = Para.at<double>(0, 0);
            p0.at<double>(1, 0) = Para.at<double>(1, 0);
            p0.at<double>(2, 0) = Para.at<double>(2, 0);
            //Calculate A, L to Form norm equation
            AbsoluteOrientation::calculate_L_matrix(i, L, mtp, mp, p0, Para);
            AbsoluteOrientation::calculate_A_matrix(i, A, Pmodel, Para);
        }
        //Solve norm equation
        X = (A.t() * A).inv() * A.t() * L;
        //cout << A << endl;
        //cout << L << endl;
        //cout << X << endl;
        //if (iteration == 3)
        //{
            //return;
        //}
        AbsoluteOrientation::correct(Para, X);
        //Adjust if achieve convergence
        if (abs(X(0, 0)) < 0.00003 && abs(X(1, 0)) < 0.00003 && abs(X(2, 0)) < 0.00003 && abs(X(3, 0)) < 0.00003 && abs(X(4, 0)) < 0.00003 && abs(X(5, 0)) < 0.00003 && abs(X(6, 0)) < 0.00003)
        {
            ofstream outfile;
            outfile.open(outputfile, ios::out);
            cout << "Convergency!!!" << endl;
            outfile << "Convergency!!!" << endl;
            cout << "Correction:" << endl;
            outfile << "Correction:" << endl;
            cout << X << endl;
            outfile << X << endl;
            cout << "--------------------------------------------" << endl;
            outfile << "--------------------------------------------" << endl;
            cout << "Absolute Orientation Result: " << endl;
            outfile << "Absolute Orientation Result: " << endl;
            cout << "Iteration: " << iteration << endl;
            outfile << "Iteration: " << iteration << endl;
            cout << "Residual:" << endl;
            outfile << "Residual:" << endl;
            Mat_<double>V = A * X - L;
            cout << V << endl;
            outfile << V << endl;
            cout << "Seven Parameters of Relative Orientation(x y z fai omega kappa s): " << endl;
            outfile << "Seven Parameters of Relative Orientation(x y z fai omega kappa s): " << endl;
            cout << Para.at<double>(0, 0) + space_center.x << " " << Para.at<double>(1, 0) + space_center.y << " " << Para.at<double>(2, 0) + space_center.z << " " << Para.at<double>(3, 0)
                << "  " << Para.at<double>(4, 0) << " " << Para.at<double>(5, 0) << " " << Para.at<double>(6, 0) * m << endl;
            outfile << Para.at<double>(0, 0) + space_center.x << " " << Para.at<double>(1, 0) + space_center.y << " " << Para.at<double>(2, 0) + space_center.z << " " << Para.at<double>(3, 0)
                << "  " << Para.at<double>(4, 0) << " " << Para.at<double>(5, 0) << " " << Para.at<double>(6, 0) * m << endl;
            Mat_<double>V_ = V.t() * V;
            double accuracy = sqrt(V_.at<double>(0, 0) / (3 * point_num - 7));
            cout << "RMS Error£º" << accuracy << endl;
            outfile << "RMS Error£º" << accuracy << endl;
            Mat_<double> Rfinal = Mat::zeros(3, 3, CV_32F);
            Mat_<double> G = Mat::zeros(3, 1, CV_32F);
            Mat_<double> result = Mat::zeros(3, 1, CV_32F);
            AbsoluteOrientation::calculate_rotation_matrix(Rfinal, Para);
            G.at<double>(0, 0) = space_center.x * m;
            G.at<double>(1, 0) = space_center.y * m;
            G.at<double>(2, 0) = space_center.z * m;
            cout << "Coordinate of ground points: " << endl;
            outfile << "Coordinate of ground points: " << endl;
            for (int i = 0; i < point_num; i++)
            {
                Mat_<double> m0 = Mat::zeros(3, 1, CV_32F);
                m0.at<double>(0, 0) = Pmodel[i].x;
                m0.at<double>(1, 0) = Pmodel[i].y;
                m0.at<double>(2, 0) = Pmodel[i].z;
                result = Para.at<double>(6, 0) * m * Rfinal * m0 + G;
                cout << pname[i] << " " << result << endl;
                outfile << pname[i] << " " << result << endl;
            }
            break;
        }
    }
}