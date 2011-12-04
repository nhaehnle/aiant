#ifndef OPPORTUNISTICATTACK_H
#define OPPORTUNISTICATTACK_H

struct Bot;
struct State;

struct OpportunisticAttack {
	OpportunisticAttack(Bot & b);
	~OpportunisticAttack();

	void init();
	void run();

	Bot & bot;
	State & state;
};

#endif // OPPORTUNISTICATTACK_H
