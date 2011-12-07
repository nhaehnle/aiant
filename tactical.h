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
	struct Data;
	struct ShadowAnt;
	struct Theater;

	Tactical(Bot & bot_);
	~Tactical();

	void init();
	void run();

	bool timeover();
	void generate_theater(const Location & center);
	void run_theater(uint theateridx);
	void make_moves();

	void gensubmap(Submap & sm, const Location & center);
	void gensubmap_field(Submap & sm, const Location & local, const Location & global);
	void compute_foodpotential(const Submap & sm, BaseSubmap<uint> & potential);

	void compute_deaths(Theater & th, PlayerMove & pm, PlayerMove & enemymove, int & mydeaths, int & enemydeaths);
	int evaluate_ant_positions(Theater & th, PlayerMove & pm, bool myperspective);
	void evaluate_new_moves(uint theateridx);
	void improve(uint theateridx, uint myidx, uint enemyidx);

	Bot & bot;
	State & state;

	Data & d;
};

#endif // TACTICAL_H
