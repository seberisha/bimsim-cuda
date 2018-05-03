#!/usr/bin/python3


import subprocess, sys

nuStart = 800
nuStop = 4000
nuStep = 10

if len(sys.argv) > 1:
	nuStart = int(sys.argv[1])
if len(sys.argv) > 2:
	nuStop = int(sys.argv[2])
if len(sys.argv) > 3:
	nuStep = int(sys.argv[3])

command = "./bimsim"

#images
intImage = "out_i.bmp"
absImage = "out_a"
#detector specs
dsize = 208
dsample = 1
res = int(dsize / dsample)
sampling = 1
padding = 1
#incident field
source = ""
order = 100
mc = 400
#sphere
x = 0
y = 0
z = 0
a = 5
n = 1.4
k = 0.0
#spectral samples
iters = int((nuStop - nuStart) / nuStep) + 1
#optics
NAin = 0.2
NAout = 0.5


#set the position of the image plane
command += " -u " + str(-dsize/2)
command += " -v " + str(-dsize/2)
command += " -w " + str(a)
command += " -U " + str(dsize/2)
command += " -V " + str(dsize/2)
command += " -W " + str(a)
command += " --plane-norm-x " + str(0)
command += " --plane-norm-y " + str(0)
command += " --plane-norm-z " + str(1)
command += " -R " + str(res)
command += " --supersample " + str(sampling)

if source != "":
	command += " -X " + source
command += " --field-order " + str(order)
command += " -d " + str(padding)
command += " -s " + str(mc)
command += " -I " + intImage
command += " -A " + absImage
command += " --append"
#sphere
command += " -x " + str(x)
command += " -y " + str(y)
command += " -z " + str(z)
command += " -r " + str(a)
command += " -n " + str(n)
command += " -k " + str(k)
#set the optics
command += " -c " + str(NAin)
command += " -C " + str(NAout)
command += " -o " + str(NAin)
command += " -O " + str(NAout)

for inu in range(0, iters):
    print("Iteration # " + str(inu + 1) + "/" + str(iters))
    nu = nuStart + inu * nuStep
    
    runcommand = command + " --nu " + str(nu)
    print(runcommand)
    subprocess.call(runcommand, shell=True)

#print("Hello world!")
