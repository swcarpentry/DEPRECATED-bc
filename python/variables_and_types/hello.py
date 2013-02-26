import time

print "hello world"

time.sleep(1)

print "is it me you're looking for?"

time.sleep(2)

f = open("richie.txt")
for line in f:
  print line.rstrip()
  time.sleep(0.2)

f.close()

