#!/usr/bin/env python

import argparse
import json
import os

class Bot:
	def __init__(self, idx, name):
		self.idx = idx
		self.name = name
		self.games = []
		self.status = {}
		self.results = []

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
	bots = [Bots[idx] for idx in botidx]

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
	loadreplays()

	maxnamelen = max([len(bot.name) for bot in Bots])
	status = set()
	for bot in Bots:
		status.update(bot.status.iterkeys())

	print "{0} games played".format(NrGames)
	print

	if args.bot:
		for bot in Bots:
			if bot.name == args.bot:
				print "Details for {0}".format(bot.name)
				details = {}
				for gameidx,botid,status in bot.games:
					if not status in details:
						details[status] = []
					details[status].append((gameidx, botid))
				for status in details.iterkeys():
					print "  {0}:".format(status),
					print ', '.join(["{0}:{1}".format(gameidx,botid) for gameidx,botid in details[status]])
				print

	statusfields = [(s, max(10, len(s))) for s in status]
	print "Bots", ' ' * maxnamelen, ' '.join(["{0:>{1}}".format(s, l) for s, l in statusfields])
	for idx in range(len(Bots)):
		bot = Bots[idx]
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

	print "    ", ' '.join(["{0:>3}".format(i+1) for i in range(len(Bots))])
	for idx in range(len(Bots)):
		print "{0:>3}:".format(idx+1),
		for otheridx in range(len(Bots)):
			if idx == otheridx:
				print "---",
			else:
				results = Bots[idx].getmatches(otheridx)
				s = sum(results)
				if s < 5:
					print "???",
				else:
					print "{0:>2}%".format(int(float(results[2]) / float(s) * 100.0)),
		print

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Collect statistics form replays.')
	parser.add_argument('--bot', metavar='BOT', type=str, help='Show details for the given bot')

	args = parser.parse_args()

	main()
