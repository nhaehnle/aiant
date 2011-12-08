# Extract timing information from a debug.txt

import re
import sys

turn=None
v={}
v["moves"]=None
v["evals"]=None
v["time"]=None

turnre = re.compile("turn (\d+):")
timere = re.compile("time taken in turn (\d+): ([0-9.]+)ms")
moveevalre = re.compile("Total number of generated moves: (\d+), evaluated pairs: (\d+)")

print ' '*3, ' '.join(['{0:10}'.format(key) for key in v.iterkeys()])

for line in sys.stdin:
	m = turnre.match(line)
	if m:
		turn=int(m.groups()[0])
		for key in v.iterkeys():
			v[key] = None
	else:
		m = timere.match(line)
		if m:
			v["time"] = m.groups()[1]
		else:
			m = moveevalre.match(line)
			if m:
				v["moves"] = m.groups()[0]
				v["evals"] = m.groups()[1]

	if turn != None:
		for key in v.iterkeys():
			if v[key] == None:
				break
		else:
			print '{0:03}'.format(turn), ' '.join('{0:>10}'.format(v[key]) for key in v.iterkeys())
			for key in v.iterkeys():
				v[key] = None
				turn = None
