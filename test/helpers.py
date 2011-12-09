# Some common helpers

import os
import stat

def listmaps(dirname):
	maps = []
	for name in os.listdir(dirname):
		path = dirname + '/' + name
		st = os.stat(path)
		if stat.S_ISDIR(st.st_mode):
			maps += listmaps(path)
		elif name.endswith(".map"):
			rec = {"path": path}
			with open(path, 'r') as filp:
				for line in filp:
					parts = line.split()
					if not parts:
						continue
					if parts[0] == "players":
						rec["players"] = int(parts[1])
			if not "players" in rec:
				print "Map {0} has no players info".format(path)
				continue

			print "Found {rec[players]} player map {rec[path]}".format(**locals())
			maps.append(rec)
	return maps

