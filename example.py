from CGMoms import *

print("Example test case.")

print("")

MU =  [3, 4]
#MU =  [0,0]
SIGMA = [[5,6],[6,7]]
#SIGMA = [[5,0],[0,5]]

alpha = [1,2]
alpha_vals = [[0,0],
              [1,0],
              [0,1],
              [2,0],
              [1,1],
              [0,2],
              [3,0],
              [2,1],
              [1,2],
              [0,3]]

print("alpha:")
print(alpha)
print("MU:")
print(MU)
print("SIGMA:")
print(SIGMA)

print("")

print("Moment of order alpha:")
print(CGMoms_Kan(alpha, MU, SIGMA))

print("")

print("alpha_vals:")
print(alpha_vals)

moms=CGMoms(alpha_vals, MU, SIGMA)
alpha_vals = np.array(alpha_vals) # ensure alpha_vals is a np.array
MU = np.array(MU) # ensure MU is a np.array
MU = np.reshape(MU, MU.size) # ensure MU is a vector
SIGMA = np.array(SIGMA)

print("Moments of order alpha_vals:")
print(np.hstack((alpha_vals,np.reshape(moms,(moms.size, 1)))))
