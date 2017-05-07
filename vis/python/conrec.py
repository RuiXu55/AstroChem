import sys
import numpy as np

# reaction length is six
rec = []
with open('reaction1.out') as hf:
  for line in hf:
    rec = np.append(rec,line.split())

with open('reaction10.out') as hf:
  for line in hf:
    rec = np.append(rec,line.split())
with open('reaction100.out') as hf:
  for line in hf:
    rec = np.append(rec,line.split())

# construct a new Speciesname output
f0 = open('ChemReactions6.txt','r')
f1 = open('ChemReactions.txt','w')

for i in range(0,13):
 line = f0.readline()
 data = line.split()
 f1.write(line)

# select species
for line in f0:
 data = line.split()
 if data[0] =='#':
   f1.write(line)
   break
 # full name
 for i in range (0,len(rec)/6):
   ind = i*6
   if data[1] == rec[ind] and data[2]==rec[ind+1] \
     and data[3] ==rec[ind+2] and data[4]==rec[ind+3] \
     and data[5] == rec[ind+4] and data[6]== rec[ind+5]:
     f1.write(line)
     break

for line in f0:
   f1.write(line)

f1.close()
f0.close()
