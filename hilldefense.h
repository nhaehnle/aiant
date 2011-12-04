#ifndef HILLDEFENSE_H
#define HILLDEFENSE_H

#include <cstdlib>

struct Bot;
struct Location;
struct State;

struct HillDefense {
	HillDefense(Bot & b);
	~HillDefense();

	void init();
	void run();

	bool hilldestroyed();
	uint getnrhills();
	const Location & gethill(uint idx);

	void update_hills();
	void update_hill_distances(uint hillidx);

	struct Data;

	Bot & bot;
	Data & d;
	State & state;
};

#endif // HILLDEFENSE_H
