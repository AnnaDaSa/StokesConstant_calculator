# Calculation algorithm of f_1

from constants import *
from custom_classes import *
from recurrence import *

# Record initial memory usage
initial_memory = get_memory_usage()

# VARIABLES :
# mv = v1 + v2 i
# w  = w1 + w2 i
# theta = theta1 + theta2 i
# J = J1 + J2 i
v1, v2, w1, w2, theta1, theta2, j1, j2 = hy.make_vars("v1","v2","w1","w2","theta1","theta2","j1","j2")

A = 1 + r1*hy.cos(theta1)*hy.cosh(theta2) + r2*hy.cos(2*theta1)*hy.cosh(2*theta2) + r3*hy.cos(3*theta1)*hy.cosh(3*theta2)
B = - r1*hy.sin(theta1)*hy.sinh(theta2) - r2*hy.sin(2*theta1)*hy.sinh(2*theta2) - r3*hy.sin(3*theta1)*hy.sinh(3*theta2)
C = - r1*hy.sin(theta1)*hy.cosh(theta2) - 2*r2*hy.sin(2*theta1)*hy.cosh(2*theta2) - 3*r3*hy.sin(3*theta1)*hy.cosh(3*theta2)
D = - r1*hy.cos(theta1)*hy.sinh(theta2) - 2*r2*hy.cos(2*theta1)*hy.sinh(2*theta2) - 3*r3*hy.cos(3*theta1)*hy.sinh(3*theta2)

# EQUATIONS OF MOTION :
dv1 = 4*(v1**2-v2**2)*w1 - 8*v1*v2*w2
dv2 = 4*(v1**2-v2**2)*w2 + 8*v1*v2*w1
dw1 = -4*v1*(w1**2 - w2**2) + 8*v2*w1*w2 - (A*(v1**3 - 3*v1*v2**2)+B*(3*v2*v1**2 - v2**3))/(4*(v1**6 + v2**6 + 3*(v1**4*v2**2 + v1**2*v2**4)))
dw2 = -4*v2*(w1**2-w2**2)-8*v1*w1*w2-(B*(v1**3-3*v1*v2**2)-A*(3*v2*v1**2-v2**3))/(4*(v1**6+v2**6+3*(v1**4*v2**2+v1**2*v2**4)))
dtheta1 = 1+nu*j1
dtheta2 = nu*j2
dj1 = (C*(v1**2-v2**2)+2*D*v1*v2)/(8*(v1**2+v2**2)**2)
dj2 = (D*(v1**2-v2**2)-2*C*v1*v2)/(8*(v1**2+v2**2)**2)


# RECURRENCE :

V = TrigPol([0,r1,r2],[0,0,0])
v0 = CustomComplex(-R,-M)
theta0 = CustomComplex(theta0r,theta0i)

result, duration = measure_time(recurrencia, K, V)
alpha, dalpha, _ = result
print(f"Done calculating T up to order K={K}")
print(f"Time taken for recurrencia: {duration:.2f} seconds")

#Look for the initial condition:
start_time = time.time()

#w0 = dv T^k(v0,theta0)
#J0 = dtheta T^k(v0,theta0)

j0 = CustomComplex(real("0",prec),real("0",prec))
w0 = CustomComplex(real("0",prec),real("0",prec))

for i in range(0,K+1):
    j0 += dalpha[i].eval(theta0)/v0**(i+1)
    w0 += -(i+1)*alpha[i].eval(theta0)/v0**(i+2)

# We create a terminal event that triggers when Re(v1)=0:

t_ev = hy.t_event(
        # The left-hand side of the event equation
        v1,
        # Specify this is an event
        # for arbitrary-precision
        # integration.
        fp_type = real)

def F(theta_r, theta_i):
    j0 = CustomComplex(real("0",prec),real("0",prec))
    w0 = CustomComplex(real("0",prec),real("0",prec))

    for i in range(0,K+1):
        j0 += dalpha[i].eval(CustomComplex(theta_r,theta_i))/v0**(i+1)
        w0 += -(i+1)*alpha[i].eval(CustomComplex(theta_r,theta_i))/v0**(i+2)

    ta = hy.taylor_adaptive(
                    # The ODEs.
                    [(v1,dv1),(v2,dv2),
                        (w1,dw1),(w2,dw2),
                        (theta1,dtheta1),(theta2,dtheta2),
                        (j1,dj1),(j2,dj2)],
                    # Initial conditions.
                    [v0.real,
                        v0.imag,
                        w0.real,
                        w0.imag,
                        theta_r,
                        theta_i,
                        j0.real,
                        j0.imag],
                    # Terminal events.
                    t_events = [t_ev],
                    # Specify that the integrator must operate in arbitrary precision.
                    fp_type = real,
                    # Explicitly specify the precision.
                    prec = prec)
    time_grid = np.linspace(0,R+0.1,20000,dtype=hy.real)
    # We are integrating until Re(v1)=0 by the t_ev.
    ta.propagate_grid(time_grid)    
    return np.array([ta.state[4],ta.state[5]],dtype=real)

# Central differences: 

J_0 = np.zeros((2,2),dtype=real)

h = hy.real('1e-10',prec)
J_0[0][0] = (F(theta0r+h,theta0i)[0]-F(theta0r-h,theta0i)[0])/(2*h)
J_0[0][1] = (F(theta0r,theta0i+h)[0]-F(theta0r,theta0i-h)[0])/(2*h)
J_0[1][0] = (F(theta0r+h,theta0i)[1]-F(theta0r-h,theta0i)[1])/(2*h)
J_0[1][1] = (F(theta0r,theta0i+h)[1]-F(theta0r,theta0i-h)[1])/(2*h)


def BroydenStep(J,f0,f1,x0,x1):
    df = f1-f0
    dx = x1-x0
    norm2 = np.dot(dx,dx)
    return J + np.outer((df-np.dot(J,dx))/norm2,dx)

def NewtonStep(xn,fn,Jn):
    return xn - np.dot(np.linalg.inv(Jn),fn)

x0 = np.array([theta0r,theta0i],dtype=real)
f0 = F(x0[0],x0[1])

x1 = NewtonStep(x0,f0,J_0)

tol = mp.power(10,-mp.dps+5)
max_iter = 100

J_1 = np.zeros((2,2),dtype=real)

# Broyden's method iteration:
for i in range(max_iter):
    f1 = F(x1[0],x1[1])
    if np.linalg.norm(f1) < tol:
        print(f"Converged in {i+1} iterations")
        break
    J_1 = BroydenStep(J_0,f0,f1,x0,x1)
    x0, f0, J_0 = x1, f1, J_1
    x1 = NewtonStep(x0, f0, J_0)
    
if(i==max_iter-1):
    print("ERROR BROYDEN'S METHOD")
else:
    thetaf = CustomComplex(x1[0],x1[1])
    #print(f"Solution found: Theta = {thetaf}")

# Recalculate w0 and J0:
j0 = CustomComplex(real("0",prec),real("0",prec))
w0 = CustomComplex(real("0",prec),real("0",prec))

for i in range(0,K+1):
    j0 += dalpha[i].eval(thetaf)/v0**(i+1)
    w0 += -(i+1)*alpha[i].eval(thetaf)/v0**(i+2)

# Lets check what we found:
ic = np.array([v0.real, v0.imag, w0.real, w0.imag,thetaf.real, thetaf.imag, j0.real,j0.imag],dtype=real)
print(f"Initial condition: {ic}")

ta_check = hy.taylor_adaptive(
                    # The ODEs.
                    [(v1,dv1),(v2,dv2),
                        (w1,dw1),(w2,dw2),
                        (theta1,dtheta1),(theta2,dtheta2),
                        (j1,dj1),(j2,dj2)],
                    # Initial conditions.
                    [v0.real,
                        v0.imag,
                        w0.real,
                        w0.imag,
                        thetaf.real,
                        thetaf.imag,
                        j0.real,
                        j0.imag],
                    # Terminal events.
                    t_events = [t_ev],
                    # Specify that the integrator must operate in arbitrary precision.
                    fp_type = real,
                    # Explicitly specify the precision.
                    prec = prec)

time_grid = np.linspace(0,R+0.1,20000,dtype=hy.real)
ta_check.propagate_grid(time_grid)

end_time = time.time()
duration = end_time - start_time

print("Initial condition found")
print(f"Time taken: {duration}")

print(f"Final time of integration: {ta_check.time}")
print(f"Final condition: {ta_check.state}")

# Aproximation of f1:
v1 = CustomComplex(ta_check.state[0],ta_check.state[1])
w1 = CustomComplex(ta_check.state[2],ta_check.state[3])
f1 = - 2 * w1.imag * np.exp(- v1.imag)
v1_norm = np.sqrt(v1.real*v1.real + v1.imag*v1.imag)

print("\n")
print(f"Result : {v1_norm}\t{f1}")
print("\n")

# Record memory usage after calculation
final_memory = get_memory_usage()

# Print the memory usage difference
print(f"Memory usage difference: {final_memory - initial_memory:.2f} MB")
