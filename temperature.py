import numpy as np
import matplotlib.pyplot as plt

## user input
check = True
while(check):
    str = input("How many grid lines per side? ")
    try:
        grid = int(str)
        check = False
    except ValueError:
        print("Please enter an integer")
        check = True
## grid is now validated to be an integer

## this will be the number of mesh points on the grid
sq = grid*grid

## create matrices for maths process
m = np.zeros([sq, sq])
t = np.zeros([sq])
b = np.zeros([sq])

## create matrix for data visualization
f = np.zeros([grid, grid])

## fill m, the transformation matrix
## each row gets 0.25 at each column that represents a bordering grid point
for i in range(sq):
    for j in range(sq):
        if i == j:
            continue
        elif int(i/grid) == int(j/grid):
            if j - i == 1 or i - j == 1:
                m[i][j] = 0.25
        elif i%grid == j%grid:
            if int(i/grid) - int(j/grid) == 1 or int(j/grid) - int(i/grid) == 1:
                m[i][j] = 0.25

## fill b, the constant matrix
## each point gets 0.25 for each of the bordering grid points that
## are part of the 1° part of border plate's border
for i in range(sq):
    if(i%grid > (grid-1)/2):
        if int(i/grid) == 0 or int(i/grid) == grid-1:
            b[i] += 0.25
        if i%grid == 0 or i%grid == grid-1:
            b[i] += 0.25
## fix midpoint values for oddly sided grids
if grid%2 == 1:
    b[int((grid-1)/2)] = 0.125
    b[sq-1 - int((grid-1)/2)] = 0.125

## iterate Jacobi's Method
for i in range(int(grid**3)):
    t = np.matmul(m,t) + b

## fill the 2d temp array to prep for image display
for i in range(sq):
    f[int(i/grid)][i%grid] = t[i]

## display cool-warm plot of the temperature distribution with 0°C minimum value and 1°C maximum value
img = plt.imshow(f, "coolwarm", vmin=0, vmax=1)
plt.show()
