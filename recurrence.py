# In this code we apply the recurrence from Theorem 5.1  to aproximate the unstable manifold T

#LIBRARIES
from constants import *
from custom_classes import *

# FORMULES :
def formula_18(alpha,dalpha,epsilon,i,k):
    ###print(f"Executing formula 18 for i+1={i+1}, k={k}")
    if(epsilon[i+1][k] == TrigPol([0],[0])):
        s1 = nu * dalpha[i] * dalpha[k]
        s2 = (real("4",prec)*(i+1)*(k+1)) * (alpha[i] * alpha[k])
        s12 = s1 + s2
        epsilon[i+1][k] = epsilon[i+2][k-1] + s12
        ###print(f"Inside formula_18: epsilon[{i+1}][{k}]")
    return 1


def formula_15_mean(alpha_mean, epsilon, k):
    ###print(f"Executing formula 15_mean for k={k}")
    if len(alpha_mean) <= k:
        alpha_mean.extend([real("0",prec)] * (k - len(alpha_mean) + 1))  # Extend the list with zeros
    alpha_mean[k] = real("1.0",prec) / (2 * PI * (k + 1)) * epsilon[2][k - 1].average()
    
    return 1

def formula_15_osc(alpha_osc, dalpha, epsilon, k):
    ##print(f"Executing formula 15_osc for k={k}")
    if len(dalpha) <= k:
        dalpha.extend([TrigPol([0], [0])] * (k - len(dalpha) + 1))  # Extend the list with zero-PolTrig
    if len(alpha_osc) <= k:
        alpha_osc.extend([TrigPol([0], [0])] * (k - len(alpha_osc) + 1))  # Extend the list with zero-PolTrig

    dalpha[k] = real("-1.",prec)* epsilon[1][k - 1]
    alpha_osc[k] = dalpha[k].integral()
    return 1

def formula_19(alpha, dalpha, epsilon, k):
    ##print(f"Executing formula 19 for k={k}")
    if(epsilon[k][k] == TrigPol([0],[0])):
        s1 = nu * dalpha[k - 1] * dalpha[k]
        s2 = (real("4.",prec) * (k + 1) * k) * alpha[k - 1] * alpha[k]
        epsilon[k][k] = s1 + s2
        #print(f"Inside formula_19: epsilon[{k}][{k}]")
    
    return 1

def formula_20(alpha,dalpha,epsilon,k):
    #print(f"Executing formula 20 for k={k}")
    if(epsilon[k+1][k] == TrigPol([0],[0])):
        s1 = nu/2 * dalpha[k] *dalpha[k]
        s2 = (real("2",prec)*(k+1)*(k+1)) * alpha[k] * alpha[k]
        epsilon[k+1][k] = s1 + s2
        #print(f"Inside formula_20: epsilon[{k+1}][{k}] ")
    
    return 1

def formula_17(alpha,epsilon,k):
    #print(f"Executing formula 17 for k={k}")
    if(epsilon[1][k] == TrigPol([0],[0])):
        s1 = epsilon[2][k-1]
        s2 = (-k-real("1.",prec)) * alpha[k]
        epsilon[1][k] = s1 + s2
        #print(f"Inside formula_17: epsilon[1][{k}]")
    
    return 1


def recurrencia(K,V):
    #print("He entrat a la funcio")
    alpha_mean = [real("-0.25",prec),0] 

    V_int = V.integral()
    alpha_osc_1 = real("0.125",prec) * V_int
    alpha_osc_0 = TrigPol([real("0.",prec)],[real("0.",prec)])

    alpha_osc = [alpha_osc_0 , alpha_osc_1] 

    alpha = [alpha_mean[0] + alpha_osc[0] , alpha_osc[1]] 

    dalpha = [TrigPol([real("0.",prec)], [real("0.",prec)]), real("0.125",prec) * V] 

    # MATRIX EPSILON :
    # epsilon_i^k = epsilon[i][k]
    epsilon = np.full((K+2, K+3), TrigPol([real("0.",prec)], [real("0.",prec)]), dtype=TrigPol)

    epsilon[1][1] = real("8.",prec) * alpha_mean[0] * alpha[1]
    epsilon[2][1] = nu/2. * dalpha[1] * dalpha[1] + real("8.",prec) * alpha[1] * alpha[1]


    # RECURRENCE :

    for k in range(2,K+1):
        #print(f"Calculating k={k}")
        if (k==2):
            formula_15_mean(alpha_mean,epsilon,2)
            formula_15_osc(alpha_osc,dalpha,epsilon,2)
            if len(alpha) <= k:
                alpha.extend([TrigPol([real("0.",prec)], [real("0.",prec)])] * (k - len(alpha) + 1))  # Extend the list with zero-PolTrig
            alpha[k] = alpha_mean[k] + alpha_osc[k]
            #print(alpha[k])
            formula_17(alpha,epsilon,k)
            formula_19(alpha,dalpha,epsilon,k)
            formula_20(alpha,dalpha,epsilon,k)
        elif (k==3):
            formula_15_mean(alpha_mean,epsilon,3)
            formula_15_osc(alpha_osc,dalpha,epsilon,3)
            if len(alpha) <= k:
                alpha.extend([TrigPol([real("0.",prec)], [real("0.",prec)])] * (k - len(alpha) + 1))  # Extend the list with zero-PolTrig
            alpha[k] = alpha_mean[k] + alpha_osc[k]
            ##print(alpha[k])
            formula_17(alpha,epsilon,k)
            formula_18(alpha,dalpha,epsilon,1,k)
            formula_19(alpha,dalpha,epsilon,k)
            formula_20(alpha,dalpha,epsilon,k)
        elif (k>3):
            formula_15_mean(alpha_mean,epsilon,k)
            formula_15_osc(alpha_osc,dalpha,epsilon,k)
            if len(alpha) <= k:
                alpha.extend([TrigPol([real("0.",prec)], [real("0.",prec)])] * (k - len(alpha) + 1))  # Extend the list with zero-PolTrig
            alpha[k] = alpha_mean[k] + alpha_osc[k]
            ##print(alpha[k])
            formula_17(alpha,epsilon,k)
            for i in range(1,k-2 +1):
                formula_18(alpha,dalpha,epsilon,i,k)
            formula_19(alpha,dalpha,epsilon,k)
            formula_20(alpha,dalpha,epsilon,k)
    return alpha, dalpha, epsilon 

            
        
        

