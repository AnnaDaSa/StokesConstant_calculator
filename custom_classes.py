# Here we create the classes CustomComplex, TrigPol and ExpPol

#LIBRARIES
from constants import *

class CustomComplex:
    def __init__(self, real, imag):
        self.real = real
        self.imag = imag

    def __add__(self, other):
        if isinstance(other, CustomComplex):
            return CustomComplex(self.real + other.real, self.imag + other.imag)
        elif isinstance(other, (int, float,hy.real)):
            return CustomComplex(self.real + other, self.imag)
        else:
            raise TypeError("Unsupported operand type for addition")

    def __sub__(self, other):
        if isinstance(other, CustomComplex):
            return CustomComplex(self.real - other.real, self.imag - other.imag)
        elif isinstance(other, (int, float,hy.real)):
            return CustomComplex(self.real - other, self.imag)
        else:
            raise TypeError("Unsupported operand type for addition")

    def __radd__(self, other):
        return self.__add__(other)
    
    def __mul__(self, other):
        if isinstance(other, CustomComplex):
            real_part = self.real * other.real - self.imag * other.imag
            imag_part = self.real * other.imag + self.imag * other.real
            return CustomComplex(real_part, imag_part)
        elif isinstance(other, (int, float,hy.real)):
            return CustomComplex(self.real * other, self.imag * other)
        else:
            raise TypeError("Unsupported operand type for multiplication")

    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __truediv__(self, other):
        if isinstance(other, (int, float, hy.real)):
            if other == 0:
                raise ZeroDivisionError("Division by zero")
            # division by a scalar
            return CustomComplex(self.real / other, self.imag / other)
        elif isinstance(other, CustomComplex):
            if other.real == 0 and other.imag == 0:
                raise ZeroDivisionError("Division by zero")
            # division by another complex
            conj = CustomComplex(other.real, -other.imag)
            denominator = (other.real * other.real + other.imag * other.imag)
            numerator = self * conj
            return CustomComplex(numerator.real / denominator, numerator.imag / denominator)
        else:
            raise TypeError("Unsupported operand type for division")

    def __rtruediv__(self, other):
        if isinstance(other, (int, float, hy.real)):
            if self.real == 0 and self.imag == 0:
                raise ZeroDivisionError("Division by zero")
            # reciplocal of a comples (1/x) and right division by a complex
            conj = CustomComplex(self.real, -self.imag)
            denominator = (self.real * self.real + self.imag * self.imag)
            numerator = CustomComplex(other, 0) * conj
            return CustomComplex(numerator.real / denominator, numerator.imag / denominator)
        else:
            raise TypeError("Unsupported operand type for division")


    def __pow__(self, exponent):
        if isinstance(exponent, (int, float, hy.real)):
            r = np.sqrt(self.real**2 + self.imag**2)
            theta = np.arctan2(self.imag, self.real)
            r_pow = r ** exponent
            theta_pow = theta * exponent
            return CustomComplex(r_pow * np.cos(theta_pow), r_pow * np.sin(theta_pow))
        else:
            raise TypeError("Unsupported operand type for power operation")

    def real_part(self):
        return self.real

    def imag_part(self):
        return self.imag
    
    def exp(self):
        r = np.sqrt(self.real**2+self.imag**2)
        theta = np.arctan2(self.imag,self.real)
        return r*CustomComplex(np.cos(theta),np.sin(theta))

    def __repr__(self):
        return f"{self.real} + {self.imag} i"
    
    def __pow__(self, exponent):
        if isinstance(exponent, int):
            # Exponentiation by an integer using repeated multiplication
            result = CustomComplex(1, 0)  
            base = self
            if exponent == 0:
                return CustomComplex(1, 0)  
            elif exponent < 0:
                base = CustomComplex(1, 0) / self  # For negative exponents, use reciprocal
                exponent = -exponent

            for _ in range(exponent):
                result *= base  # Multiply base repeatedly

            return result
        else:
            raise TypeError("Unsupported operand type for power operation")
    
    
class TrigPol:
    def __init__(self, coef_cos, coef_sin):
        self.s = coef_sin
        self.c = coef_cos
        #a_k = a[k] for k=0..n
        #b_k = b[k]
        # they must have the same lenght:
        # append zeros if necessary
        self.n = len(coef_cos) - 1

    def __add__(self, other):
        if isinstance(other, TrigPol):
            n1 = self.n
            cos_1 = self.c
            sin_1 = self.s
            n2 = other.n
            cos_2 = other.c
            sin_2 = other.s

            cos_3 = np.zeros(max(n1, n2) + 1, dtype=hy.real)
            sin_3 = np.zeros(max(n1, n2) + 1, dtype=hy.real)
            
            if n1 <= n2:
                for k in range(n1 + 1):
                    cos_3[k] = cos_1[k] + cos_2[k]
                    sin_3[k] = sin_1[k] + sin_2[k]
                for k in range(n1 + 1, n2 + 1):
                    cos_3[k] = cos_2[k]
                    sin_3[k] = sin_2[k]
            else:
                for k in range(n2 + 1):
                    cos_3[k] = cos_1[k] + cos_2[k]
                    sin_3[k] = sin_1[k] + sin_2[k]
                for k in range(n2 + 1, n1 + 1):
                    cos_3[k] = cos_1[k]
                    sin_3[k] = sin_1[k]


            result = TrigPol(cos_3, sin_3)
            return result
        elif isinstance(other, (int, float,hy.real)):
            a = self.c
            b = self.s
            a[0]+= other
            result = TrigPol(a,b)
            return result 
        else:
            raise TypeError("Unsupported operand type for addition")

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        if isinstance(other,TrigPol):
            self_exp = Trig_to_Exp(self)
            other_exp = Trig_to_Exp(other)
            result_exp = self_exp * other_exp
            result_trig = Exp_to_Trig(result_exp)
            return result_trig
        elif isinstance(other,(int,float,hy.real)):
            a = self.c
            b = self.s
            trig_pol = TrigPol([other*x for x in a], [other*x for x in b])
            return trig_pol
        else:
            raise TypeError("Unsupported operand type for multiplication")
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def average(self):
        a = self.c
        return 2*PI*a[0]
    
    def integral(self):
        a = self.c
        b = self.s
        n = self.n
        c = np.zeros(n+1,dtype=hy.real) 
        d = np.zeros(n+1,dtype=hy.real)

        for k in range(1,n+1):
            c[k] = -b[k]/k
            d[k] = a[k]/k
            
        result = TrigPol(c,d)
        return result 

    def __eq__(self, other):
        if isinstance(other, TrigPol):
            return np.array_equal(self.c, other.c) and np.array_equal(self.s, other.s)
        return False

    def __str__(self):
        return f"TrigPol(c={self.c}, s={self.s})"

    def __repr__(self):
        return f"TrigPol(c={self.c}, s={self.s})"

    def eval(self, t):
        if isinstance(t,CustomComplex):
            # t = a + bi
            a = t.real
            b = t.imag
            sum_val = 0.
            for k in range(self.n + 1):
                # cos(a+bi)=cos(a)cosh(b)-i sin(a)sinh(b)
                # sin(a+bi)=sin(a)cosh(b)+i cos(a)sinh(b)
                cskt = CustomComplex(np.cos(k*a)*np.cosh(k*b),-np.sin(k*a)*np.sinh(k*b))
                snkt = CustomComplex(np.sin(k*a)*np.cosh(k*b), np.cos(k*a)*np.sinh(k*b))
                sum_val += self.c[k] * cskt + self.s[k] * snkt
            return sum_val
        else:
            sum_val = 0.
            for k in range(self.n + 1):
                sum_val += self.c[k] * np.cos(k * t) + self.s[k] * np.sin(k * t)
            return sum_val

class ExpPol:
    def __init__(self,coef):
        # Coefficients go from -n to n
        # so ours will go from 0 to 2n
        self.c = coef 
        #c_k = c[n+k] for k=-n..n
        self.n = (len(coef)-1) // 2 # len(coef) = 2n+1
    
    def __mul__(self,other):
        if isinstance(other,ExpPol):
            n1 = self.n
            n2 = other.n
            m = n1+n2 
            c = self.c
            d = other.c
            e = np.zeros(2*m+1, dtype=CustomComplex)
            for i in range(2*n1+1):
                for j in range(2*n2+1):
                    e[i+j] += c[i] * d[j]
            exp_pol = ExpPol(e)
            return exp_pol
        else:
            raise TypeError("Unsupported operand type for multiplication")
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __str__(self):
        return f"ExpPol({self.c})"

    def __repr__(self):
        return f"ExpPol({self.c})"

def Trig_to_Exp(trig_pol):
    a = trig_pol.c
    b = trig_pol.s 
    n = trig_pol.n
    c = np.zeros(2*n+1, dtype=object)
    for k in range(0, n):
        c[k] = CustomComplex(a[n-k], b[n-k]) / 2
    c[n] = a[0]
    for k in range(1, n+1):
        c[n+k] = CustomComplex(a[k], -b[k]) / 2
    exp_pol =ExpPol(c)
    return exp_pol

def Exp_to_Trig(exp_pol):
    c = exp_pol.c
    n = exp_pol.n
    a = [c[n].real_part() if isinstance(c[n], CustomComplex) else c[n]]
    b = [0]
    for k in range(1,n+1):
        a.append(hy.real("2.",prec) * c[n+k].real_part())
        b.append(-hy.real("2.",prec) * c[n+k].imag_part())
   
    trig_pol =TrigPol(a,b)
    return trig_pol
