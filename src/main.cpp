#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include <vector>
#include "RungeKutta.h"

int main() {
    auto window = sf::RenderWindow{{800, 800}, "Double Pendulum"};
    RungeKutta rk(1.0, 2.0, 1.0, 2.0);
    rk.set_state({M_PI / 2.0, 0, 0, 0});

    double ct = 0;
    for (int i = 0; i < 100; ++i) {
        std::vector<double> curr = rk.get_current_state();

        double x1 = cos(curr[0]);
        double y1 = sin(curr[0]);

        double x2 = x1 + cos(curr[0] + curr[1]);
        double y2 = y1 + sin(curr[0] + curr[1]);

        std::cout << x1 << ' ' << y1 << ' ' << x2 << ' ' << y2 << '\n';

        rk.set_state(rk.get_next(0, curr));
        ct += rk.get_dt();
    }

    while (window.isOpen()) {
        for (auto event = sf::Event{}; window.pollEvent(event);) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        window.display();
    }
}
