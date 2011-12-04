#ifndef OFFENSE_H
#define OFFENSE_H

#include <cstdlib>

struct Bot;
struct State;

struct Offense {
	Offense(Bot & b);
	~Offense();

	void init();
	void run();

	void update_hill_distance(uint hillidx);
	void update_hills();
	uint closesthill();

	struct Data;

	Bot & bot;
	State & state;
	Data & d;
};

#endif // OFFENSE_H
