import time

print "hello world"

time.sleep(1)

print "Dr. Grace Hopper, 1906-1992, creator of the first compiler"

time.sleep(2)

f = open("hopper.txt")
for line in f:
  print line.rstrip()
  time.sleep(0.2)

f.close()

