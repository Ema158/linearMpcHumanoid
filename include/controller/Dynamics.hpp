#pragma once
#include "controller/Robot.hpp"
#include "controller/generalizedFunctions.hpp"
#include "controller/linkInertia.hpp"
#include "controller/invKinematics.hpp"
#include <Eigen/Dense>
#include <vector>

class Dynamics{
    public:
    Dynamics() = default;

    void computeAll(const Robot& robot);

    const Eigen::VectorXd& getC() const {return C_;}
    const Eigen::VectorXd& getCg() const {return Cg_;}
    const Eigen::MatrixXd& getM() const {return M_;}
    const Eigen::MatrixXd& getAG() const {return AG_;}
    const Eigen::VectorXd& getAGpqp() const {return AGpqp_;}
    
    private:
    void allSpatialInertiaMatrices(const Robot& robot); //computes the spatial inertia 6x6 matrices of all frames

    Eigen::VectorXd computeC(const Robot& robot,
        bool isGravity); //..computes the Coriolles,centrigula,gravitational force vector.
                        //it is the forward Newton-Euler with vD=0

    void computeM(const Robot& robot);//Composite rigid body algorithm

    void centroidalMatrixAndBias(const Robot& robot); //computes the centroidal matrix of the centroidal model

    Eigen::MatrixXd spatialInertiaMatrix(const linkInertia& link); //compues the spatial inertia 6x6 matrix of a frame

    void forwardNewtonEuler(const Robot& robot,
        const Eigen::VectorxXd qD,
        const Eigen::VectorXd qDD,
        std::vector<Eigen::VectorXd>& vel,
        std::vector<Eigen::VectorXd>& acc,
        std::vector<Eigen::VectorXd>& f);

    Eigen::VectorXd C_;
    Eigen::VectorXd Cg_; //Vector C without gravity, needed for centroidal dynamics
    Eigen::MatrixXd M_;
    Eigen::MatrixXd AG_;
    Eigen::VectorXd AGpqp_;
    std::vector<Eigen::MatrixXd> I_;
};