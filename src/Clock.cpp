#include "linearMpcHumanoid/general/Clock.hpp"

Clock::Clock(double dt, double simulationTime)
    :
    dt_(dt),
    T_(simulationTime)
{

}