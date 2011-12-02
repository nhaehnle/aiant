#ifndef SCOUT_H

struct Bot;
struct State;

struct Scout {
	Scout(Bot & b);
	~Scout();

	void init();
	void run();

	void recompute_maps();
	void update_nrscouts();

	struct Data;

	Data & d;
	Bot & bot;
	State & state;

private:
	Scout(const Scout & o);
	Scout & operator=(const Scout & o);
};

#endif // SCOUT_H
