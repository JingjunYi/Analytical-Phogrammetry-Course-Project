#include "RelativeOrientation.h"

using namespace std;
using namespace cv;


void RelativeOrientation::read_file_for_relative_orientation(string file, string img1, string img2, vector<Ppoint>& P1, vector<Ppoint>& P2, vector<double>& intri1, vector<double>& intri2, vector<int>& pname)
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
        if (i == 1)
        {
            string name[2];
            while (str >> name[0] >> name[1])
            {
                img1 = name[0];
                img2 = name[1];
            }
        }
        else if (i == 2)
        {
            string inp[4];
            while (str >> inp[0] >> inp[1] >> inp[2] >> inp[3])
            {
                intri1.push_back(atof(inp[0].c_str()) / 1000.0);
                intri1.push_back(atof(inp[1].c_str()) / 1000.0);
                intri1.push_back(atof(inp[2].c_str()) / 1000.0);
                intri1.push_back(atof(inp[3].c_str()));
            }
        }
        else if (i == 3)
        {
            string inp[4];
            while (str >> inp[0] >> inp[1] >> inp[2] >> inp[3])
            {
                intri2.push_back(atof(inp[0].c_str()) / 1000.0);
                intri2.push_back(atof(inp[1].c_str()) / 1000.0);
                intri2.push_back(atof(inp[2].c_str()) / 1000.0);
                intri2.push_back(atof(inp[3].c_str()));
            }
        }
        else
        {
            string split[5];
            while (str >> split[0] >> split[1] >> split[2] >> split[3] >> split[4])
            {
                int id = atoi(split[0].c_str());
                double x1 = atof(split[1].c_str()) / 1000.0;
                double y1 = atof(split[2].c_str()) / 1000.0;
                double x2 = atof(split[3].c_str()) / 1000.0;
                double y2 = atof(split[4].c_str()) / 1000.0;
                P1.push_back(Ppoint(x1, y1));
                P2.push_back(Ppoint(x2, y2));
                pname.push_back(id);
            }
        }
    }
    cout << "Model: " << img1 << "-" << img2 << endl;
    cout << P1.size() << " Points" << endl;
}

void RelativeOrientation::calculate_relarotation_matrix(Mat_<double>& R2, Mat_<double>& Para)
{
    double fai = Para.at<double>(0,0);
    double omega = Para.at<double>(1, 0);
    double kappa = Para.at<double>(2, 0);
    R2.at<double>(0, 0) = cos(fai) * cos(kappa) - sin(fai) * sin(omega) * sin(kappa);
    R2.at<double>(0, 1) = -cos(fai) * sin(kappa) - sin(fai) * sin(omega) * cos(kappa);
    R2.at<double>(0, 2) = -sin(fai) * cos(omega);
    R2.at<double>(1, 0) = cos(omega) * sin(kappa);
    R2.at<double>(1, 1) = cos(omega) * cos(kappa);
    R2.at<double>(1, 2) = -sin(omega);
    R2.at<double>(2, 0) = sin(fai) * cos(kappa) + cos(fai) * sin(omega) * sin(kappa);
    R2.at<double>(2, 1) = -sin(fai) * sin(kappa) + cos(fai) * sin(omega) * cos(kappa);
    R2.at<double>(2, 2) = cos(fai) * cos(omega);
}

void RelativeOrientation::calculate_A_matrix(int i, Mat_<double>& A, Gpoint& P2, double N1, double N2, double Bx)
{
    double a[5] = {};
    a[0] = -P2.x * P2.y * N2 / P2.z; //q
    a[1] = -(P2.z + P2.y * P2.y / P2.z) * N2; //w
    a[2] = P2.x * N2; //k
    a[3] = Bx; //u
    a[4] = -P2.y * Bx / P2.z; //v
    for (int j = 0; j < 5; j++)
    {
        A.at<double>(i, j) = a[j];
    }
}

void RelativeOrientation::calculate_L_matrix(int i, Mat_<double>& L, double Q)
{
    L.at<double>(i, 0) = Q;
}

void RelativeOrientation::correct(Mat_<double>& Para, Mat_<double>& X)
{
    for (int i = 0; i < 5; i++)
    {
        Para.at<double>(i, 0) += X.at<double>(i, 0);
    }
}

void RelativeOrientation::relative_orientation(string file, string outputfile)
{
    ////Define and Initiate parameter
    int iteration = 0;
    string img1, img2;
    vector<Ppoint> P1, P2;
    vector<double> intri1, intri2;
    vector<int> pname;
    RelativeOrientation::read_file_for_relative_orientation(file, img1, img2, P1, P2, intri1, intri2, pname);
    double f, f_;
    double x0, x0_, y0, y0_;
    double m1, m2;
    f = intri1[0];
    f_ = intri2[0];
    x0 = intri1[1];
    x0_ = intri2[1];
    y0 = intri1[2];
    y0_ = intri2[2];
    m1 = intri1[3];
    m2 = intri1[3];
    //Calculate scale
    double d1 = sqrt(pow((P1[0].x - P1[1].x), 2) + pow((P1[0].y - P1[1].y), 2));
    double d2 = sqrt(pow((P2[0].x - P2[1].x), 2) + pow((P2[0].y - P2[1].y), 2));
    double m = d2 / d1;
    //Initiate
    int point_num = P1.size();
    Mat_<double> Para = Mat::zeros(5, 1, CV_32F); //fai, omega, kappa, u, v 
    Mat_<double> X = Mat::zeros(5, 1, CV_32F); //dfai, domega, dkappa, du, dv 
    Mat_<double> A = Mat::zeros(point_num, 5, CV_32F);
    Mat_<double> L = Mat::zeros(point_num, 1, CV_32F);
    vector<Gpoint> Aux1(point_num);
    vector<Gpoint> Aux2(point_num);
    vector<Gpoint> model(point_num);
    double Bx = P1[0].x - P2[0].x;
    double By = 0, Bz = 0;

    ////Calculate parameter with Iteration
    while (true)
    {
        iteration += 1;
        cout << "Iteration: " << iteration << endl;
        model.clear();
        Mat_<double> R2 = Mat::zeros(3, 3, CV_32F);
        RelativeOrientation::calculate_relarotation_matrix(R2, Para);
        By = Bx * Para.at<double>(3, 0);
        Bz = Bx * Para.at<double>(4, 0);
        //Form norm equation of each point
        for (int i = 0; i < point_num; i++)
        {
            //Calculate Aux coordinate
            Aux1[i].x = P1[i].x;
            Aux1[i].y = P1[i].y;
            Aux1[i].z = -f;
            Mat_<double> tmp = Mat::zeros(3, 1, CV_32F);
            tmp.at<double>(0, 0) = P2[i].x;
            tmp.at<double>(1, 0) = P2[i].y;
            tmp.at<double>(2, 0) = -f_;
            Mat_<double> auxi2 = R2 * tmp;
            Aux2[i].x = auxi2.at<double>(0, 0);
            Aux2[i].y = auxi2.at<double>(1, 0);
            Aux2[i].z = auxi2.at<double>(2, 0);
            //Calculate N1, N2
            double N1 = 0, N2 = 0;
            N1 = (Bx * Aux2[i].z - Bz * Aux2[i].x) / (Aux1[i].x * Aux2[i].z - Aux2[i].x * Aux1[i].z);
            N2 = (Bx * Aux1[i].z - Bz * Aux1[i].x) / (Aux1[i].x * Aux2[i].z - Aux2[i].x * Aux1[i].z);
            //Calculate Q of each point pair
            double Q = N1 * Aux1[i].y - N2 * Aux2[i].y - By;
            //Calculate A, L
            RelativeOrientation::calculate_A_matrix(i, A, Aux2[i], N1, N2, Bx);
            RelativeOrientation::calculate_L_matrix(i, L, Q);
            //Calculate model coordinates
            double mx = m * N1 * Aux1[i].x;
            double my = 0.5 * m * (N1 * Aux1[i].y + N2 * Aux2[i].y + By);
            double mz = m * f + m * N1 * (-f);
            model.push_back(Gpoint(mx, my, mz));
        }
        //Calculate correct value of Para
        X = (A.t() * A).inv() * A.t() * L;
        /*cout << A << endl;
        cout << L << endl;
        cout << X << endl;
        return;*/
        RelativeOrientation::correct(Para, X);

        ////Ouput Result
        //Adjust if reach convergency
        //3e-5
        if (abs(X.at<double>(0, 0)) < 0.00003&& abs(X.at<double>(1, 0)) < 0.00003 && abs(X.at<double>(2, 0)) < 0.00003 && abs(X.at<double>(3, 0)) < 0.00003 && abs(X.at<double>(4, 0)) < 0.00003)
        {
            cout << "Convergency!!!" << endl;
            cout << "Correction:" << endl;
            cout << X << endl;
            cout << "--------------------------------------------" << endl;
            cout << "Relative Orientation Result: " << endl;
            cout << "Iteration: " << iteration << endl;
            cout << "Residual:" << endl;
            Mat_<double>V = A * X - L;
            cout << V << endl;
            cout << "Five Parameters of Relative Orientation(fai, omega, kappa, u, v): " << endl;
            cout << Para.at<double>(0, 0) << " " << Para.at<double>(1, 0) << " " << Para.at<double>(2, 0) << " " << Para.at<double>(3, 0) << "  " << Para.at<double>(4, 0) << endl;
            Mat_<double>V_ = V.t() * V;
            double accuracy = sqrt(V_.at<double>(0, 0) / (point_num - 5));
            cout << "RMS Error£º" << accuracy << endl;
            cout << "Model coordinates(x y z):" << endl;
            for (int i = 0; i < point_num; i++)
            {
                cout << pname[i] << " " << model[i].x << " " << model[i].y << " " << model[i].z << endl;
            }
            ofstream outfile;
            outfile.open(outputfile, ios::out);
            outfile << "Residual:" << endl;
            outfile << V << endl;
            outfile << "Five Parameters of Relative Orientation(fai, omega, kappa, u, v): " << endl;
            outfile << Para.at<double>(0, 0) << " " << Para.at<double>(1, 0) << " " << Para.at<double>(2, 0) << " " << Para.at<double>(3, 0) << "  " << Para.at<double>(4, 0) << endl;
            outfile << "RMS Error£º" << accuracy << endl;
            outfile << "Iteration: " << iteration << endl;
            outfile << "Model coordinates(x y z):" << endl;
            for (int i = 0; i < point_num; i++)
            {
                outfile << pname[i] << " " << model[i].x << " " << model[i].y << " " << model[i].z << endl;
            }
            break;
        }
    }
}