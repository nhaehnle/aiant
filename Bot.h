#ifndef BOT_H_
#define BOT_H_

#include "State.h"

struct FoodSeeker;
struct HillDefense;
struct Offense;
struct OpportunisticAttack;
struct Scout;
struct SymmetryFinder;
struct Zoc;

struct Ant {
	Location where;
	int direction;
	const int * dirperm;

	bool hastactical;

	Ant() : dirperm(0) {reset();}

	void reset();
};

/*
    This struct represents your bot in the game of Ants
*/
struct Bot
{
	State state;

	Bot();
	~Bot();

	void playGame();    //plays a single game of Ants

	void makeMoves();   //makes moves for a single turn
	void endTurn();     //indicates to the engine that it has made its moves

	void update_ants();
	uint myantidx_at(const Location & pos);

	bool try_rotate_move(uint antidx);
	void make_moves();

	struct Data;

	Data & d;
	Zoc & m_zoc;
	SymmetryFinder & m_symmetry;
	FoodSeeker & m_foodseeker;
	Scout & m_scout;
	HillDefense & m_hilldefense;
	OpportunisticAttack & m_opportunisticattack;
	Offense & m_offense;
	std::vector<Ant> m_ants;
};

#endif //BOT_H_
