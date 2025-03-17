#include "RungeKutta.h"

RungeKutta::RungeKutta(double r1, double m1, double r2, double m2) {
    this->_r1 = r1;
    this->_m1 = m1;

    this->_r2 = r2;
    this->_m2 = m2;

    this->_g = 9.8;
    this->_dt = 0.05;
    this->_ct = 0;
}

double RungeKutta::get_theta1() {
    return this->_theta1;
}

double RungeKutta::get_theta2() {
    return this->_theta2;
}

double RungeKutta::get_omega1() {
    return this->_omega1;
}

double RungeKutta::get_omega2() {
    return this->_omega2;
}

std::vector<double> RungeKutta::calc(std::vector<double>& state) {
    std::vector<double> res;

    double theta1 = state[0];
    double theta2 = state[1];
    double omega1 = state[2];
    double omega2 = state[3];

    res.emplace_back(omega1);
    res.emplace_back(omega2);

    const double delta = theta1 - theta2;
    const double common_denom = 2 * _m1 + _m2 * (1 - cos(2 * delta));

    double omega1_prime =
        -1.0 * _g * (2 * _m1 + _m2) * sin(theta1) -
        _m2 * _g * sin(delta - theta2) -
        2 * sin(delta) * _m2 *
            (_r2 * omega2 * omega2 + _r1 * omega1 * omega1 * cos(delta));
    omega1_prime /= (_r1 * common_denom);

    double omega2_prime =
        2 * sin(delta) * (_r1 * omega1 * omega1 * (_m1 + _m2)) +
        cos(theta1) * _g * (_m1 + _m2) +
        _r2 * _m2 * omega2 * omega2 * cos(delta);
    omega2_prime /= (_r2 * common_denom);

    res.emplace_back(omega1_prime);
    res.emplace_back(omega2_prime);

    return res;
}

std::vector<double> RungeKutta::adjust(std::vector<double>& curr,
                                       std::vector<double>& values,
                                       double factor = 1.0) {
    std::vector<double> res;

    for (int i = 0; i < 4; ++i) {
        res.emplace_back(curr[i] + factor * values[i]);
    }

    return res;
}

std::vector<double> RungeKutta::get_next() {
    std::vector<double> current_state = get_current_state();
    std::vector<double> k1 = calc(current_state);
    k1 = adjust(current_state, k1, _dt / 2.0);
    std::vector<double> k2 = calc(k1);
    k2 = adjust(current_state, k2, _dt / 2.0);
    std::vector<double> k3 = calc(k2);
    k3 = adjust(current_state, k3, _dt);
    std::vector<double> k4 = calc(k3);

    k4 = adjust(k4, k3, 2);
    k4 = adjust(k4, k2, 2);
    k4 = adjust(k4, k1);

    _ct += _dt;
    std::vector<double> res = adjust(current_state, k4, _dt / 6.0);
    set_state(res);
    return res;
}

std::vector<double> RungeKutta::get_current_state() {
    return {_theta1, _theta2, _omega1, _omega2};
}

void RungeKutta::set_state(std::vector<double> state) {
    this->_theta1 = state[0];
    this->_theta2 = state[1];
    this->_omega1 = state[2];
    this->_omega2 = state[3];
}

double RungeKutta::get_dt() {
    return this->_dt;
}
