#include <SDL2/SDL.h>
#include <cstdlib>
#include <ctime>
#include <emscripten.h>
#include <vector.hpp>
#include <vector>

using std::cout;
using std::endl;

// Doesn't really matter if it's not the real G
// Can be used to tweak gravitational power
const double G = 0.01;

const int NUM_PLANETS = 50;
const int MAX_MASS = 700;

// Tweak radius of planets
const double SIZE_CONST = 1;

const int CANVAS_W = 1000;
const int CANVAS_H = 700;

class Planet {
  private:
	double mass;
	Vector velocity;
	Vector acceleration;

  public:
	// Constructor
	Planet(double mass, Vector position, Vector velocity)
	    : mass(mass), position(position), velocity(velocity),
	      radius(cbrt(mass) * SIZE_CONST) {}
	Planet(double mass, Vector position)
	    : Planet(mass, position, Vector(0, 0)) {}

	double radius;
	Vector position;

	/**
	 * Calculates gravitational force on this planet due to another and
	 * updates the acceleration
	 */
	void gravitate(const Planet &other) {
		// Force = (G * m1 * m2) / r^2
		//    => Acc = (G * m2) / r^2
		double r_squared = pow(position.distance(other.position), 2);
		double delta_acc_mag = (G * other.mass) / r_squared;

		Vector direction_vector = other.position - position;
		double delta_acc_dir = direction_vector.direction();

		Vector delta_acc;
		delta_acc.set_as_polar(delta_acc_mag, delta_acc_dir);

		acceleration = acceleration + delta_acc;
	}

	/**
	 * Checks for a collision between this planet and another. If there's a hit,
	 * merge the other planet into this one and return true, else return false
	 */
	bool checkCollision(const Planet &other) {
		// Check spherical collision condition
		if (position.distance(other.position) <= radius + other.radius) {
			// Merge while conserving momentum
			Vector momentum = velocity * mass;
			Vector other_momentum = other.velocity * other.mass;
			Vector new_momentum = momentum + other_momentum;

			// Find center of mass
			double new_mass = mass + other.mass;
			Vector new_position =
			    ((position * mass) + (other.position * other.mass)) / new_mass;

			velocity = new_momentum / new_mass;
			radius = cbrt(new_mass) * SIZE_CONST;
			mass = new_mass;
			position = new_position;

			return true;
		}
		return false;
	}

	/**
	 * tick..
	 */
	void update() {
		velocity = velocity + acceleration;
		position = position + velocity;

		// Reset, to calculate net force again on the next iteration
		acceleration = Vector(0, 0);
	}
};

struct Universe {
	std::vector<Planet> planets;
	SDL_Renderer *renderer;
};

void NewDrawCircle(SDL_Renderer *renderer, Vector position, int radius) {
	for (int i = 0; i < radius * 2; ++i) {
		for (int j = 0; j < radius * 2; ++j) {
			int di = radius - i;
			int dj = radius - j;
			if (pow(di, 2) + pow(dj, 2) <= pow(radius, 2)) {
				SDL_RenderDrawPoint(renderer, position.x + di, position.y + dj);
			}
		}
	}
}

// Helper method to draw a circle
void DrawCircle(SDL_Renderer *renderer, Vector position, int radius) {
	const int diameter = (radius * 2);

	int x = (radius - 1);
	int y = 0;
	int tx = 1;
	int ty = 1;
	int error = (tx - diameter);
	double cx = position.x;
	double cy = position.y;

	while (x >= y) {
		//  Each of the following renders an octant of the circle
		SDL_RenderDrawPoint(renderer, cx + x, cy - y);
		SDL_RenderDrawPoint(renderer, cx + x, cy + y);
		SDL_RenderDrawPoint(renderer, cx - x, cy - y);
		SDL_RenderDrawPoint(renderer, cx - x, cy + y);
		SDL_RenderDrawPoint(renderer, cx + y, cy - x);
		SDL_RenderDrawPoint(renderer, cx + y, cy + x);
		SDL_RenderDrawPoint(renderer, cx - y, cy - x);
		SDL_RenderDrawPoint(renderer, cx - y, cy + x);

		if (error <= 0) {
			++y;
			error += ty;
			ty += 2;
		}

		if (error > 0) {
			--x;
			tx += 2;
			error += (tx - diameter);
		}
	}
}

void mainloop(void *arg) {
	Universe *universe = static_cast<Universe *>(arg);
	SDL_Renderer *renderer = universe->renderer;

	auto &planets = universe->planets;

	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
	SDL_RenderClear(renderer);

	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

	// Calculate G forces
	for (size_t i = 0; i < planets.size(); ++i) {
		for (size_t j = 0; j < planets.size(); ++j) {
			if (i == j)
				continue;
			planets[i].gravitate(planets[j]);
		}
	}

	// Update positions
	for (auto &planet : planets) {
		planet.update();
	}

	// Check for collisions
	for (size_t i = 0; i < planets.size() - 1; ++i) {
		for (size_t j = i + 1; j < planets.size(); j++) {
			if (planets[i].checkCollision(planets[j])) {
				planets.erase(planets.begin() + j);
				j -= 1;
			}
		}
	}

	for (auto &planet : planets) {
		DrawCircle(renderer, planet.position, planet.radius);
	}

	SDL_RenderPresent(renderer);
}

int main() {
	srand(time(0));
	SDL_Init(SDL_INIT_VIDEO);
	SDL_Window *window;
	SDL_Renderer *renderer;
	SDL_CreateWindowAndRenderer(CANVAS_W, CANVAS_H, 0, &window, &renderer);

	const int HALF_W = CANVAS_W / 2;
	const int HALF_H = CANVAS_H / 2;

	std::vector<Planet> planets;
	for (int i = 0; i < NUM_PLANETS; ++i) {
		double mass = rand() % MAX_MASS;
		Vector position = Vector((rand() % HALF_W) + HALF_W / 2,
		                         (rand() % HALF_H) + HALF_H / 2);
		int VEL = 200;
		Vector velocity =
		    Vector(((rand() % VEL) / (double)100) - (VEL / (double)200),
		           ((rand() % VEL) / (double)100) - (VEL / (double)200));
		planets.push_back(Planet(mass, position, velocity));
	}

	Universe universe;
	universe.renderer = renderer;
	universe.planets = planets;

	emscripten_set_main_loop_arg(mainloop, &universe, -1, 1);

	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();

	return EXIT_SUCCESS;
}
