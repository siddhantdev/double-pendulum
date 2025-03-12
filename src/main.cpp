#include <SFML/Graphics.hpp>
#include <SFML/System/Sleep.hpp>
#include <cmath>
#include <vector>
#include "RungeKutta.h"
#include "const.h"

// TODO: make an adjust function which takes the co-ordinates and shifts them
std::pair<double, double> adjust(double x, double y) {
    int cx = Const::WIDTH / 2;
    int cy = Const::HEIGHT / 2;
    x = cx + (x * Const::WIDTH / 6.0);
    y = cy - (y * Const::HEIGHT / 6.0);

    return {x, y};
}

int main() {
    sf::RenderWindow window(sf::VideoMode(Const::WIDTH, Const::HEIGHT),
                            "Double Pendulum Simulation", sf::Style::Default,
                            sf::ContextSettings(0, 0, 8));
    window.setFramerateLimit(10);
    RungeKutta rk(1.0, 2.0, 1.0, 2.0);
    std::vector<double> curr = {M_PI / 2.0, 0, 0, 0};
    rk.set_state(curr);

    double ct = 0;

    while (window.isOpen()) {
        for (auto event = sf::Event{}; window.pollEvent(event);) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }
        window.clear();

        curr = rk.get_next(ct, curr);
        ct += rk.get_dt();

        double x1 = cos(curr[0]);
        double y1 = sin(curr[0]);

        double x2 = (x1 + cos(curr[0] + curr[1]));
        double y2 = (y1 + sin(curr[0] + curr[1]));

        std::tie(x1, y1) = adjust(x1, y1);
        std::tie(x2, y2) = adjust(x2, y2);

        sf::CircleShape circle1(Const::RADIUS);
        circle1.setOrigin(Const::RADIUS, Const::RADIUS);
        circle1.setPosition(x1, y1);
        circle1.setFillColor(sf::Color(100, 0, 0));

        sf::CircleShape circle2(Const::RADIUS);
        circle2.setOrigin(Const::RADIUS, Const::RADIUS);
        circle2.setPosition(x2, y2);
        circle2.setFillColor(sf::Color(0, 0, 100));

        sf::Vertex line1[2];
        line1[0].position =
            sf::Vector2f(Const::WIDTH / 2.0, Const::HEIGHT / 2.0);
        line1[0].color = sf::Color::White;
        line1[1].position = sf::Vector2f(x1, y1);
        line1[1].color = sf::Color::White;

        sf::Vertex line2[2];
        line2[0].position = sf::Vector2f(x1, y1);
        line2[0].color = sf::Color::White;
        line2[1].position = sf::Vector2f(x2, y2);
        line2[1].color = sf::Color::White;

        window.draw(line1, 2, sf::Lines);
        window.draw(line2, 2, sf::Lines);

        window.draw(circle1);
        window.draw(circle2);

        window.display();
    }
}
