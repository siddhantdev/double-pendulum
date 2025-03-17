#pragma once

#include <cassert>
#include <cmath>
#include <vector>

class RungeKutta {
   public:
    RungeKutta();
    RungeKutta(double r1, double m1, double r2, double m2);

    double get_theta1();
    double get_theta2();

    double get_omega1();
    double get_omega2();

    // state is of the format : {theta1, theta2, omega1, omega2}
    std::vector<double> calc(std::vector<double>& state);
    std::vector<double> get_next(std::vector<double>& current_state);
    std::vector<double> get_current_state();

    void set_state(std::vector<double> state);

    double get_dt();

   private:
    double _g;

    double _theta1;
    double _theta2;

    double _omega1;
    double _omega2;

    double _r1;
    double _r2;

    double _m1;
    double _m2;

    double _dt;
    double _ct;

    std::vector<double> adjust(std::vector<double>& curr,
                               std::vector<double>& values, double factor);
};
