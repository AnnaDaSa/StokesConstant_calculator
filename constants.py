# Here I store my constants all together 

# LIBRARIES
import numpy as np 
import mpmath as mp # PI
import heyoka as hy # Arbitrary precision arithetics and integrator
import psutil # For memory usage
import time # For time consumption

# CONSTANTS

prec=200 # Precision (in bits)

real = hy.real # Funtion real from heyoka
mp.pretty = True # Option for the type of output from mpmath
mp.dps = int(prec * np.log(2)/np.log(10)) # Digits of PI

PI = real(str(mp.pi),prec)

a = real("3.6",prec)
alpha_ctt = real("1.05",prec)
# r1 = real("0.",prec)
# r2 = real("1.",prec)
# r3 = real("1.",prec)
r1 = real("0.06",prec)
r2 = real("0.0008",prec)
r3 = real("0.0",prec)

nu = ((4*PI)/(a*alpha_ctt))**2

K = 60 

print(f"Precision in bits: {prec}")
print(f"Precision in digits: {mp.dps}")

#v_0=-R-iM
R = real("30.",prec)
M = real("40.",prec)
#print(f"e^(-M)={np.exp(-M)}")

#theta_0=theta0r + theta0i
theta0r = -R
theta0i = real("-0.000005",prec)


# FUNCTIONS

# Function to get current memory usage
def get_memory_usage():
    process = psutil.Process()
    mem_info = process.memory_info()
    return mem_info.rss / 1024**2  # Return memory in MB

# Function to measure execution time
def measure_time(func, *args, **kwargs):
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    duration = end_time - start_time
    return result, duration

# Hamiltonian 
def real_hamiltonian(x):
    Ah = 1+r1*np.cos(x[4])*np.cosh(x[5])+r2*np.cos(2*x[4])*np.cosh(2*x[5])
    Bh = -r1*np.sin(x[4])*np.sinh(x[5])-r2*np.sin(2*x[4])*np.sinh(2*x[5])
    E = x[6]+nu/2*(x[6]**2-x[7]**2)+2*(x[0]**2-x[1]**2)*(x[2]**2-x[3]**2)\
        -8*x[0]*x[1]*x[2]*x[3]-(Ah*(x[0]**2-x[1]**2)+2*Bh*x[0]*x[1])/(8*(x[0]**2+x[1]**2)**2)
    return E

def im_hamiltonian(x):
    Ah = 1+r1*np.cos(x[4])*np.cosh(x[5])+r2*np.cos(2*x[4])*np.cosh(2*x[5])
    Bh = -r1*np.sin(x[4])*np.sinh(x[5])-r2*np.sin(2*x[4])*np.sinh(2*x[5])
    E = x[7]+nu*x[6]*x[7]+4*((x[0]**2-x[1]**2)*x[2]*x[3]+(x[2]**2-x[3]**2)*x[0]*x[1])\
        -(-2*Ah*x[0]*x[1]+Bh*(x[0]**2-x[1]**2))/(8*(x[0]**2+x[1]**2)**2)
    return E