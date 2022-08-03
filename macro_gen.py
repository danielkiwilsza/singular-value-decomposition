import math

a = 0
while(a <= 2*math.pi):
    print("{:.5f}f,".format(math.cos(a)))
    a += math.pi / 64