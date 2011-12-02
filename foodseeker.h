#ifndef FOODSEEKER_H
#define FOODSEEKER_H

#include <cstdlib>

struct Bot;
struct Location;
struct State;

struct FoodSeeker {
	FoodSeeker(Bot & b);
	~FoodSeeker();

	void init();
	void run();

	uint foodidx_at(const Location & pos);
	void assign_food();

	struct Data;

	Bot & bot;
	State & state;
	Data & d;
};

#endif // FOODSEEKER_H
