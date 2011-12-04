#ifndef BOT_H_
#define BOT_H_

#include "State.h"
#include "map.h"

struct FoodSeeker;
struct HillDefense;
struct OpportunisticAttack;
struct Scout;
struct Zoc;

struct Ant {
	Location where;

	bool hastactical;

	int direction;

	Ant() : hastactical(false), direction(-1) {}
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

	uint myantidx_at(const Location & pos);

	bool try_rotate_move(uint antidx, const Map<bool> & claims);
	void make_moves();

	Zoc & m_zoc;
	FoodSeeker & m_foodseeker;
	Scout & m_scout;
	HillDefense & m_hilldefense;
	OpportunisticAttack & m_opportunisticattack;
	std::vector<Ant> m_ants;
};

#endif //BOT_H_
