#ifndef BOT_H_
#define BOT_H_

#include "State.h"

struct Diffusion;
struct FoodSeeker;
struct HillDefense;
struct Module;
struct Offense;
struct OpportunisticAttack;
struct Scout;
struct SymmetryFinder;
struct Zoc;

struct Ant {
	Location where;
	bool assigneddirection;
	int direction;
	const int * dirperm;

	Ant() : dirperm(0) {reset();}

	void reset();
};

/*
    This struct represents your bot in the game of Ants
*/
struct Bot
{
	State state;

	Bot(int argc, char * * argv);
	~Bot();

	float getargfloat(const std::string & key, float def);

	void playGame();    //plays a single game of Ants

	void makeMoves();   //makes moves for a single turn
	void endTurn();     //indicates to the engine that it has made its moves

	void update_ants();
	uint myantidx_at(const Location & pos);

	bool try_rotate_move(uint antidx);
	void make_moves();

	bool timeover();

	struct Data;

	Data & d;
	Zoc & m_zoc;
	SymmetryFinder & m_symmetry;
	FoodSeeker & m_foodseeker;
	Scout & m_scout;
	HillDefense & m_hilldefense;
	OpportunisticAttack & m_opportunisticattack;
	Diffusion & m_diffusion;
	Offense & m_offense;
	Module & m_tactical;
	std::vector<Ant> m_ants;
};

#endif //BOT_H_
