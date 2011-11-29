#ifndef TACTICAL_H
#define TACTICAL_H

#include <stdlib.h>

struct Bot;
struct Location;
struct Outcome;
struct PlayerMove;
template<typename T>
struct BaseSubmap;
struct Submap;
struct State;

struct Tactical {
	Tactical(Bot & bot_);
	~Tactical();

	void gensubmap(Submap & sm, const Location & center);
	void gensubmap_field(Submap & sm, const Location & local, const Location & global);

	void compute_foodpotential(const Submap & sm, BaseSubmap<uint> & potential);
	void evaluate_food(const Submap & sm, int & myvalue, int & enemyvalue);
	int evaluate(const Submap & sm, int & myvalue, int & enemyvalue);
	void make_moves(const Location & center);

	void generate_initmoves();
	void complete_move(PlayerMove & pm);
	bool evaluate_new_moves();
	void improve(uint myidx, uint enemyidx);
	int try_collect_food(const Location & antpos);
	void apply_moves(Submap & sm, const PlayerMove & me, const PlayerMove & enemy, Outcome & outcome);

	Bot & bot;
	State & state;

	struct Data;
	Data & d;

private:
	Tactical(const Tactical &);
	Tactical & operator=(const Tactical &);
};

#endif // TACTICAL_H
