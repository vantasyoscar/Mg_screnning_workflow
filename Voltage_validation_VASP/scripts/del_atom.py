#!/home/cwj/opt/anaconda3/bin/python
from ase.io import read,write
from ase.build import sort
import sys
import random
fname=sys.argv[1]
ele=sys.argv[2]
ndel=int(sys.argv[3])
savefname=sys.argv[4]

a=sort(read(fname))
for ind in range(len(a)):
    if(a[ind].symbol==ele):
        begin=ind
        break
for ind in range(begin,len(a)+1):
    try:
        if(a[ind].symbol!=ele):
            end=ind
            break
    except:
        end=len(a)
        break
if(end-begin<ndel):
    print("Error!!!")
    print("The atom number (%d) of element %s is less than the number (%d) you want delete!"
            %(end-begin,ele,ndel))
    exit()

list=random.sample(range(begin,end),ndel)

del a[list]

write(savefname,a)