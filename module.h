#ifndef MODULE_H
#define MODULE_H

struct Module {
	Module() {}
	virtual ~Module() {}

	virtual void init() = 0;
	virtual void run() = 0;
	virtual void learn() {}

private:
	Module(const Module &);
	Module & operator=(const Module &);
};

#endif // MODULE_H
