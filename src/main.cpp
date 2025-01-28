#include <SFML/Graphics.hpp>

int main() {
    auto window = sf::RenderWindow{{800, 800}, "Double Pendulum"};
    while (window.isOpen()) {
        for (auto event = sf::Event{}; window.pollEvent(event);) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        window.display();
    }
}
