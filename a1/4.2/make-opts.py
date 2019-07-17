import itertools

f = open('out.txt', 'w')

#opts = ["-march=native", "-fomit-frame-pointer", "-floop-block", "-floop-interchange", "-floop-strip-mine", "-funroll-loops", "-flto"]
opts = ["-march=native", "-xHost", "-unroll", "-ipo"]
a=range(0,len(opts))
f.write('\n'.join('\n'.join((''.join(str(i) for i in c) + ' ' + ' '.join(opts[i] for i in c)) for c in itertools.combinations(a, i)) for i in range(1, len(a)+1)))

f.close()
