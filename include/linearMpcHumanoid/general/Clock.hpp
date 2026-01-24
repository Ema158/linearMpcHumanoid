#pragma once

class Clock{
    public:
        double getTime() const {return t_;}
        double getTimeStep() const {return dt_;}
        
        void step() {t_ += dt_;}

    private:
        double t_ = 0;
        double dt_ = 0.01;

};
