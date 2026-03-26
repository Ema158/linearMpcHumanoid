#pragma once 

enum class Task{
    Stand,
    Walk,
    Jump
};

struct GaitParameters{
    int numSteps = 1;
    double timePerStep = 0.5;
    double currentXPos = 0;
    double currentYPos = 0;
    double futureXPos = 0.00;
    double futureYPos = -0.05;
    double stepHeight = 0.02;
    double stepLength = 0.025;
    double stepWidth = 0.1;
    double timeInitialDS = 0.4;
};
    