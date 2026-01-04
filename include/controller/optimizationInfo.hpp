#pragma once

class optimizationInfo{
public:
    //==================Whole-Body weights
    double wCoML = 0; //linear momentum weight
    double wCoMK = 10000; //angular momentum rate weight
    double wBasePos = 10; //base position
    double wBaseAng = 10; //base attitude
    double wJoints = 1; //rotational joints
    double wForce = 1; //reaction forces
    double wFoot = 100000; //position and orientation of both feet
};