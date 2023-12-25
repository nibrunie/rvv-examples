def BfromA(x, y):
    return y, (2 * x + 3 * y) % 5

for j in range(5):
    for i in range(5):
        y = i
        xDual = (j - 3 * y)
        # x = {0: 0, 2: 1, 4: 2, 3: 4, 1: 3}[xDual]
        x = (xDual * 3) % 5 # 3 is 1/2 in the field
        #print(f"B[{i},{j}] = A[{x},{y}]", "")
        print(f"{(x + 5 * y) * 8}, ", end="")
    print("")

# python script to expand rotation amounts in B coordinates    
x = 1; y = 0;
Arots = {}
Brots = {}
# Iterate over ((0 1)(2 3))^t * (1 0) for 0 ≤ t ≤ 23 
for t in range(24):
    # Compute the rotation constant r = (t+1)(t+2)/2
    r = ((t+1)*(t+2)//2)%64
    Arots[(x,y)] = r
    print(f"A[{x, y}] <<< {r} -> ", end="")
    # Compute ((0 1)(2 3)) * (x y)
    tmpY = (2*x+3*y)%5; x = y; y = tmpY
    print(f"B[{x, y}]")
    Brots[(x,y)] = r

Brots[(0, 0)] = 0
print("B rotations")
for j in range(5):
    for i in range(5):
        print(f"{Brots[(i, j)]}, ", end="")
    print("")