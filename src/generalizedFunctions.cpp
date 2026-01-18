#include "controller/generalizedFunctions.hpp"

Eigen::Matrix3d crossMatrix(Eigen::Vector3d v){
    Eigen::Matrix3d A = Eigen::Matrix3d::Zero(3,3); 
    A <<    0, -v(2),  v(1),
         v(2),     0, -v(0), 
        -v(1),  v(0),     0;
    return A;
}

Eigen::MatrixXd velocityMatrix(Eigen::Matrix4d T){
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(6,6);
    Eigen::Matrix3d R = T.block(0,0,3,3).transpose(); //Rotation matrix
    Eigen::Vector3d p = T.block(0,3,3,1); //position vector
    X.block(0,0,3,3) = R;
    X.block(3,0,3,3) = -R*crossMatrix(p);
    X.block(3,3,3,3) = R;
    return X;
}

Eigen::Matrix4d inverseTransformationMatrix(Eigen::Matrix4d T){
    Eigen::Matrix4d Tinv = Eigen::Matrix4d::Zero(4,4);
    Tinv.block(0,0,3,3) = T.block(0,0,3,3).transpose();
    Tinv.block(0,3,3,1) = (-T.block(0,0,3,3).transpose())*T.block(0,3,3,1);
    Tinv(3,3) = 1;
    return Tinv;
}

Eigen::MatrixXd spatialCrossMatrix(Eigen::VectorXd v){
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(6,6);
    m.block(0,0,3,3) = crossMatrix(v.segment(0,3));
    m.block(3,0,3,3) = crossMatrix(v.segment(3,3));
    m.block(3,3,3,3) = crossMatrix(v.segment(0,3));
    return m;
}

Eigen::MatrixXd spatialCrossMatrixForce(Eigen::VectorXd v){
    Eigen::MatrixXd f = Eigen::MatrixXd::Zero(6,6);
    f = -spatialCrossMatrix(v).transpose();
    return f;
}

Eigen::Matrix3d matrixAngularVelToEulerDot(Eigen::Vector3d eta){
    //eta = [roll,pitch,yaw]
    Eigen::Matrix3d Omega = Eigen::Matrix3d::Zero(3,3);
    Omega << std::cos(eta(2))/std::cos(eta(1)),  std::sin(eta(2))/std::cos(eta(1)), 0,
                             -std::sin(eta(2)),                   std::cos(eta(2)), 0,
              std::cos(eta(2))*std::tan(eta(1)), std::sin(eta(2))*std::tan(eta(1)), 1;
    return Omega;
}

Eigen::Matrix3d eulerAnglesToSO3(const Eigen::Vector3d& rpy)
{
    const double roll  = rpy(0);
    const double pitch = rpy(1);
    const double yaw   = rpy(2);

    const double cr = std::cos(roll);
    const double sr = std::sin(roll);
    const double cp = std::cos(pitch);
    const double sp = std::sin(pitch);
    const double cy = std::cos(yaw);
    const double sy = std::sin(yaw);

    Eigen::Matrix3d R;

    R <<  cy * cp,  cy * sp * sr - sy * cr,  cy * sp * cr + sy * sr,
          sy * cp,  sy * sp * sr + cy * cr,  sy * sp * cr - cy * sr,
          -sp,      cp * sr,                 cp * cr;

    return R;
}

Eigen::Vector3d rotMatToAxisAngle(const Eigen::Matrix3d& R)
{
    // Trace and clamp for numerical robustness
    double tr = R.trace();
    double c = std::max(-1.0, std::min(1.0, (tr - 1.0) / 2.0));
    double phi = std::acos(c);

    // Skew-symmetric part
    Eigen::Matrix3d Skew = R - R.transpose();

    // vee operator
    Eigen::Vector3d v;
    v << Skew(2,1),
         Skew(0,2),
         Skew(1,0);

    Eigen::Vector3d r;
    if (phi < 1e-6){
        // Small-angle approximation: sin(phi) ≈ phi
        r = 0.5 * v;
    }
    else{
        // General case
        r = (phi / (2.0 * std::sin(phi))) * v;
    }

    return r;
}

Eigen::VectorXd findPolyCoeff(Eigen::MatrixXd pos,
    Eigen::MatrixXd vel,
    Eigen::MatrixXd acc)
{
    int posRows = pos.rows();
    int posCols = pos.cols();

    int velRows = vel.rows();
    int velCols = vel.cols();

    int accRows = acc.rows();
    int accCols = acc.cols();

    if(posCols!=2 || velCols!=2 || accCols!=2){
        std::cout<< "Number of columns should be 2 for Position Velocity and Acceleration" << std::endl;
    }
    int n = posRows + velRows + accRows; //number of coefficients
    Eigen::VectorXd Coeff = Eigen::VectorXd::Zero(n);

    Eigen::MatrixXd A(n, n);
    Eigen::VectorXd b(n); //Linear system to solve for the coefficients

    //Fill the Matrix A and vector b
    int row = 0; //current row of A and b

    for(int i=0; i<posRows; i++){
        double tPow = 1;
        for(int j=0; j<n; j++){
            A(row,j) = tPow;
            tPow *= pos(i,0);
        }
        b(row) = pos(i,1);
        row++;
    }

    for(int i=0; i<velRows; i++){
        double tPow = 1;
        A(row,0) = 0;
        for(int j=1; j<n; j++){
            A(row,j) = j*tPow;
            tPow *= vel(i,0);
        }
        b(row) = vel(i,1);
        row++;
    }

    for(int i=0; i<accRows; i++){
        double tPow = 1;
        A(row,0) = 0;
        A(row,1) = 0;
        for(int j=2; j<n; j++){
            A(row,j) = j*(j-1)*tPow;
            tPow *= acc(i,0);
        }
        b(row) = acc(i,1);
        row++;
    }

    Coeff = A.colPivHouseholderQr().solve(b);
    return Coeff;
}

double polyval(const Eigen::VectorXd& poly, double x)
{
    double value = 0;
    double xPow = 1;   

    for (int i=0;i<poly.size();i++){
        value += poly(i)*xPow;
        xPow *= x;
    }

    return value;
}

const Eigen::VectorXd polyder(const Eigen::VectorXd& poly)
{
    int n = poly.size();

    if (n <= 1) {
        return Eigen::VectorXd::Zero(1); // derivative of constant
    }

    Eigen::VectorXd polyD(n-1);

    for (int i = 0; i < n - 1; ++i) {
        polyD(i) = (i + 1) * poly(i + 1);
    }

    return polyD;
}