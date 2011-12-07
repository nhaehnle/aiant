#ifndef TACTICAL_H
#define TACTICAL_H

#include <stdlib.h>

#include "module.h"

struct Bot;
struct Location;
struct Outcome;
struct PlayerMove;
template<typename T>
struct BaseSubmap;
struct Submap;
struct State;

struct Tactical : Module {
	Tactical(Bot & bot_);
	~Tactical();

	void init();
	void run();

	void gensubmap(Submap & sm, const Location & center);
	void gensubmap_field(Submap & sm, const Location & local, const Location & global);

	void compute_foodpotential(const Submap & sm, BaseSubmap<uint> & potential);
	void make_moves(const Location & center);

	void generate_initmoves();
	void compute_deaths(PlayerMove & pm, PlayerMove & enemymove, int & mydeaths, int & enemydeaths);
	int evaluate_ant_positions(PlayerMove & pm, bool myperspective);
	bool evaluate_new_moves();
	void improve(uint myidx, uint enemyidx);

	Bot & bot;
	State & state;

	struct Data;
	Data & d;
};

#endif // TACTICAL_H
