import sys
import numpy as np

crit = 2.0

# species need to be added
sp = []
with open('slct_species.out','r') as f:
  for row in f:
    data = row.split()
    if float(data[0])>crit:
      sp.append(data[1])
      print sp[-1]

print len(sp)
f = open('SpeciesNames.txt','w')
f1= open('SpeciesNames4.txt','r')

# save format
for i in range(0,28):
  line = f1.readline()
  f.write(line)

# check add species
while(line != ''):
  line = f1.readline()
  data = line.split()
  if (line != ''):
    print data[0]
    for i in range(0,len(sp)-1):
      print data[0],sp[i]
      if data[0] == sp[i]:
        f.write(line)
        break



