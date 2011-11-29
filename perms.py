
def makeperms(lst):
  if not lst:
    return [[]]
  out = []
  for idx in range(len(lst)):
    sub = makeperms(lst[:idx] + lst[idx+1:])
    out += [ [lst[idx]] + s for s in sub ]
  return out

def prettyprint(lst):
  for perm in lst:
    print "\t{", ', '.join([str(elt) for elt in perm]), "},"

prettyprint(makeperms([-1,0,1,2,3]))
