#!/usr/bin/env python

import argparse
import json
import os

import trueskill

class Bot:
	def __init__(self, idx, name):
		self.idx = idx
		self.name = name
		self.games = []
		self.status = {}
		self.results = []
		self.skill = (trueskill.INITIAL_MU, trueskill.INITIAL_SIGMA)

	def addgame(self, gameidx, botid, status):
		self.games.append((gameidx, botid, status))
		if status in self.status:
			self.status[status] += 1
		else:
			self.status[status] = 1

	def addmatch(self, vs, result):
		if vs >= len(self.results):
			self.results += [ [0] * 3 for i in range(vs - len(self.results) + 1) ]
		self.results[vs][result+1] += 1

	def getmatches(self, vs):
		"""
		Return [losses, ties, wins] against the vs player
		"""
		if vs >= len(self.results):
			self.results += [ [0] * 3 for i in range(vs - len(self.results) + 1) ]
		return self.results[vs]

NrGames = 0
Bots = []

def getbotidx(name):
	for idx in range(len(Bots)):
		if Bots[idx].name == name:
			return idx
	Bots.append(Bot(len(Bots), name))
	return len(Bots)-1

def process_replay(replayidx, replay):
	global NrGames

	if "error" in replay:
		return False

	botidx = [getbotidx(name) for name in replay["playernames"]]

	for gameidx in range(len(botidx)):
		bot = Bots[botidx[gameidx]]
		bot.addgame(replayidx, gameidx, replay["status"][gameidx])
		myscore = replay["score"][gameidx]
		myrank = replay["rank"][gameidx]
		for othergameidx in range(len(botidx)):
			otheridx = botidx[othergameidx]
			if botidx[gameidx] == otheridx:
				continue

			otherscore = replay["score"][othergameidx]
			otherrank = replay["rank"][othergameidx]

			if myscore > otherscore:
				assert myrank < otherrank
				bot.addmatch(otheridx, +1)
			elif myscore == otherscore:
				assert myrank == otherrank
				bot.addmatch(otheridx, 0)
			else:
				assert myrank > otherrank
				bot.addmatch(otheridx, -1)

	bots = [Bots[botidx[i]] for i in range(len(botidx))]
	for i in range(len(bots)):
		bots[i].rank = replay["rank"][i]
	trueskill.AdjustPlayers(bots)

	NrGames += 1
	return True

def loadreplays():
	replayindices = [int(name.split('.')[0]) for name in os.listdir("game_logs") if name.endswith(".replay")]
	replayindices.sort()
	for replayidx in replayindices:
		with open("game_logs/{0}.replay".format(replayidx), 'r') as filp:
			try:
				replay = json.load(filp)
			except ValueError:
				if replayidx == replayindices[-1]:
					continue # happens if we look at statistics while the games are running
				print "Error parsing replay", replayidx
				continue

			if not process_replay(replayidx, replay):
				print "Error in replay", replayidx

def main():
	trueskill.SetParameters(beta=25.0/3.0, draw_probability=0.3)
	loadreplays()

	maxnamelen = max([len(bot.name) for bot in Bots])
	status = set()
	for bot in Bots:
		status.update(bot.status.iterkeys())

	botmap = [i for i in range(len(Bots))]
	botmap.sort(cmp=lambda i,j: cmp(Bots[i].name, Bots[j].name))

	print "{0} games played".format(NrGames)
	print

	if args.bot:
		for bot in Bots:
			if bot.name == args.bot:
				print "Details for {0}".format(bot.name)
				details = {}
				for gameidx,botid,s in bot.games:
					if not s in details:
						details[s] = []
					details[s].append((gameidx, botid))
				for s in details.iterkeys():
					print "  {0}:".format(s),
					print ', '.join(["{0}:{1}".format(gameidx,botid) for gameidx,botid in details[s]])
				print

	statusfields = [(s, max(10, len(s))) for s in status]
	print "Bots", ' ' * maxnamelen, ' '.join(["{0:>{1}}".format(s, l) for s, l in statusfields])
	for idx in range(len(Bots)):
		bot = Bots[botmap[idx]]
		print "{0:>3}: {bot.name:<{maxnamelen}}".format(idx+1, **locals()),
		for s, l in statusfields:
			count = 0
			if s in bot.status:
				count = bot.status[s]
			print "{0:>{1}}".format(
				"{0} ({1}%)".format(count, int(100 * count / len(bot.games))),
				l),
		print
	print

	draws = 0
	total = 0

	print "    ", ' '.join(["{0:>3}".format(i+1) for i in range(len(Bots))])
	for idx in range(len(Bots)):
		print "{0:>3}:".format(idx+1),
		for otheridx in range(len(Bots)):
			if idx == otheridx:
				print "---",
			else:
				results = Bots[botmap[idx]].getmatches(botmap[otheridx])
				s = sum(results)
				total += s
				draws += results[1]
				if s < 5:
					print "???",
				else:
					print "{0:>2}%".format(int(float(results[2]) / float(s) * 100.0)),
		print
	print

	print "Draw percentage:", 100.0 * draws / total

	print "Rank", ' '*maxnamelen, "TS    ", "mu    ", "sigma "
	bots = [i for i in range(len(Bots))]
	bots.sort(key=lambda bot: Bots[bot].skill[0] - 3*Bots[bot].skill[1], reverse=True)
	for rank in range(len(bots)):
		bot = Bots[bots[rank]]
		print "{0:>3}: {name:{maxnamelen}} {trueskill:<6.5} {mu:<6.5} {sigma:<6.5}".format(rank+1,
			name=bot.name, maxnamelen=maxnamelen,
			trueskill=bot.skill[0] - 3*bot.skill[1],
			mu=bot.skill[0], sigma=bot.skill[1])

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Collect statistics form replays.')
	parser.add_argument('--bot', metavar='BOT', type=str, help='Show details for the given bot')

	args = parser.parse_args()

	main()
