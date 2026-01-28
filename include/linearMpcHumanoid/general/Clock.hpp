#pragma once

class Clock{
    public:
        Clock(double dt, double simulationTime);

        double getTime() const {return t_;}
        double getTimeStep() const {return dt_;}
        double getSimulationTime() const {return T_;}
        
        void step() {t_ += dt_;}

    private:
        double t_ = 0;
        double dt_;
        double T_;

};
