#ifndef TACTICAL_H
#define TACTICAL_H

#include <stdlib.h>

struct Bot;
struct Location;
struct Outcome;
struct PlayerMove;
struct Scenarios;
struct Submap;
struct State;

struct Tactical {
	Tactical(Bot & bot_);

	void gensubmap(Submap & sm, const Location & center);
	void gensubmap_field(Submap & sm, const Location & local, const Location & global);

	int evaluate(const Submap & sm);
	void make_moves(const Location & center);

	void generate_moves(Scenarios & scn, const Location & offset, bool predefined);
	bool evaluate_new_moves(Scenarios & scn);
	void improve(Scenarios & scn, uint myidx, uint enemyidx);
	void apply_moves(Submap & sm, const PlayerMove & me, const PlayerMove & enemy, Outcome & outcome);

	Bot & bot;
	State & state;
};

#endif // TACTICAL_H
