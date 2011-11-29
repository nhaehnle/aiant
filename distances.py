
radius = 8
map = [ [c*c + r*r  for c in range(-radius, radius+1)] for r in range(-radius, radius+1)]

for row in map:
  print ''.join(['%4d' % (v) for v in row])
