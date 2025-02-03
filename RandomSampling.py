from sage.all import *
import random
from time import time

def algorithm1(b,bb,U,B):
    n = b.nrows() - 1
    for i in range(1,n+1):
        bb[i] = b[i]
        for j in range(1,i + 1):
            if i == j:
                U[i,j] = 1
                continue
            U[i,j] = b[i].dot_product(bb[j]) / (bb[j].norm() ** 2)
            bb[i] = bb[i] - U[i,j] * bb[j]
    for i in range(1,n+1):
        B[i] = bb[i].norm() ** 2

n = 40
b0v = zero_matrix(QQ,1,n+1)
b0r = zero_matrix(QQ,n,1)
b1 = zero_matrix(QQ,n)
for i in range(n):
    j = 0
    while j <= i:
        if i >= j:
            b1[i,j] = randint(-128,128)
            if b1[i,j] != 0:
                j = j + 1
b = b0r.augment(b1)
b = b0v.stack(b)
b[1,1] = 128


bb = zero_matrix(QQ,n+1)
U = zero_matrix(QQ,n+1)
B = zero_vector(QQ,n+1)
v = zero_vector(QQ,n+1)
p = zero_vector(QQ,n+1)
s = zero_vector(QQ,n+1)
u = 4
g = 0.6

start_time = time()

algorithm1(b,bb,U,B)

num = 0
while num < 2**24 :
    v[n] = 1
    p[n] = B[n]
    for i in range(n-1,0,-1):
        sum = 0
        for h in range(i+1,n):
            sum = sum + v[h] * U[h,i]
        s[i] = sum + U[n,i]
        if i < n-u:
            v[i] = -round(s[i])
        else:
            j = math.ceil(-1-s[i])
            candidates = [j,j+1]
            v[i] = random.choice(candidates)
            candidates.remove(v[i])
            tt = candidates[0]
        s[i] = s[i] + v[i]
        p[i] = p[i+1] + (s[i]**2)*B[i]
    if p[1] < g * B[1]:
        print(v)
        break
    num = num + 1

end_time = time()

if num == 1000:
    print("未発見")
execution_time = end_time - start_time
print(v.norm()**2)
print(execution_time)
