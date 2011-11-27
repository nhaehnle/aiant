#ifndef TACTICAL_H
#define TACTICAL_H

struct Bot;
struct Location;
struct Submap;
struct State;

static const unsigned int TacticalProximity = 6;

struct Tactical {
	Tactical(Bot & bot_);

	void gensubmap(Submap & sm, const Location & center);
	void gensubmap_field(Submap & sm, const Location & local, const Location & global);

	int evaluate(const Submap & sm);
	void make_moves(const Location & center);

	Bot & bot;
	State & state;
};

#endif // TACTICAL_H
