import sys
import numpy as np

# combine all the necessary species
dir = '../matlab/Data/'
threshold = 0.00001
name     = []
for i in range(0,3):
  r = [1,10,100]
  filename = 'r'+str(r[i])+'_z0.out'
  with open(dir+filename) as hf:
    for line in hf:
      data   = line.split()
      weight = float(data[0])
      if weight > threshold:
        name   = np.append(name,data[1])
  
print (name,len(name))

# construct a new Speciesname output
f0 = open('SpeciesNames4.txt','r')
f1 = open('SpeciesNames.txt','w')
while(1):
 line = f0.readline()
 data = line.split()
 f1.write(line)
 if data[0] =='#' and data[1] =='Species':
   break

# select species
for line in f0:
 data = line.split()
 # full name
 in_list = False
 for i in range (0,len(name)):
   if data[0] == name[i]:
     f1.write(line)
     in_list = True
     break
 # ion's neutral counterpart
 if(in_list==False):
   # add the neutral species
   # if not already in the list
   if(in_list==False):
     for i in range (0,len(name)):
       ion = str(name[i])
       ion = ion.replace('+','')
       if data[0] == ion:
         f1.write(line)
         break

f1.close()
f0.close()
