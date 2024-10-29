# Extrapolation if f_1 values

from constants import *
from custom_classes import *
from recurrence import *
import pandas as pd #For data in file functions
import matplotlib.pyplot as plt

# Lets take n values of f1:
# h_i is for 1/v_i and a_i are the approximations calculated with h_i

aa = []
h = []

# Open the file and read it line by line
with open('dades.txt', 'r') as file:
    for line in file:
        # Split each line into two parts and convert them to Decimal
        values = line.split()
        if len(values) == 2:  # Ensure there are two values on the line
            aa.append(real(values[1],prec))
            h.append(1./real(values[0],prec))
n=len(h)

a = np.zeros((n,n), dtype = real)
for i in range(n):
    a[0,i] = aa[i]

print(f"n={n}")

# EXTRAPOLATION:

# A[j,i] = [A_{i,0}^j,...,A_{i,n}^j]
A = np.empty((n,n),dtype = object) 
for j in range(n):
    for i in range(n):
        A[j, i] = np.zeros(n,dtype=real)

# A_{i,n}^0 = h_i^n-1:
for i in range(n):
    for k in range(n):
        A[0,i][k] = h[i]**(k)

for j in range(1,n):
    for i in range(n-j):
        for k in range(n):
            A[j,i][k] = (A[j-1,i+1][j]*A[j-1,i][k] - A[j-1,i][j]*A[j-1,i+1][k]) / (A[j-1,i+1][j] - A[j-1,i][j])

for j in range(1,n):
    for i in range(n-j):
        a[j,i] = (A[j-1,i+1][j] * a[j-1,i] - A[j-1,i][j] * a[j-1,i+1]) / (A[j-1,i+1][j] - A[j-1,i][j])


print(f'\n Results:\n')
print('Right diagonal:')
for j in range(n):
    print(a[j,n-1-j])

print('1st step:')
for j in range(n):
    print(a[0,j])

print('2nd step:')
for j in range(n):
    print(a[1,j])

print('3rd step:')
for j in range(n):
    print(a[2,j])


k_values = np.arange(n)

# Plot the coefficients
plt.figure(figsize=(8, 6))

# Function to get ordinal suffix
def ordinal(n):
    if 10 <= n % 100 <= 20:
        suffix = 'th'
    else:
        suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')
    return str(n) + suffix


for j in range(4):
    plt.plot(k_values[j+1:], np.abs(A[j,3])[j+1:], marker='o', linestyle='-', label=f'{ordinal(j+1)} Step')

plt.yscale('log')  
plt.title("Coefficients A[j,3][k]")
plt.xlabel("k")
plt.ylabel("log(|A[j,3][k]|)")
plt.grid(True)
plt.legend()
plt.show()