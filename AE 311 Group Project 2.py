#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, pi, sin, acos

m, p, theta = symbols('m, p, theta')

alpha_lo = -m/(pi * p**2) * (-2 * p * theta + 3/2 * theta + 2*p*sin(theta) + 1/4*sin(2*theta)                               - 2 * sin(theta))            - m/(pi * (1-p)**2) * (2 * p * theta - 3/2 * theta - 2*p*sin(theta) - 1/4*sin(2*theta) +                                    2 * sin(theta) - 2*pi*p + 3/2*pi)

def calculate_theta(p):
    if p == 0:
        return np.pi
    return np.arccos(1 - 2 * p)

theta_1 = calculate_theta(0.4)
theta_2 = calculate_theta(0.3)
theta_3 = calculate_theta(0)

alo1 = alpha_lo.subs({p: 0.4, m: 0.03, theta: theta_1})
alo2 = alpha_lo.subs({p: 0.3, m: 0.03, theta: theta_2})
alo3 = alpha_lo.subs({p: 0, m: 0, theta: theta_3})
alo4 = 0.0089858 - 0.162990/2

#cl = 2pi(alpha - alphalo)
def cl(alpha, alo):
    return 2*pi*(alpha - alo)

#cl= 2pi(Ao + A1/2 )
def cl4412(alpha, alo):
    return 2*pi*((alpha - 0.0089858) + 0.162990/2)

alpha_start = np.deg2rad(-5)
alpha_stop = np.deg2rad(10)


alpha = np.linspace(alpha_start, alpha_stop, 100)

cl1 = cl(alpha, alo1)
cl2 = cl(alpha, alo2)
cl3 = cl(alpha, alo3)
cl4 = cl4412(alpha, alo4)

plt.figure(figsize=(14, 8))
plt.plot(np.rad2deg(alpha), cl1, label='Airfoil 1 (p=0.4, m=0.03)')
plt.plot(np.rad2deg(alpha), cl2, label='Airfoil 2 (p=0.3, m=0.03)')
plt.plot(np.rad2deg(alpha), cl3, label='Airfoil 3 NACA 0012')
plt.plot(np.rad2deg(alpha), cl4, label='Airfoil 4 NACA 4412')
plt.xlabel('Angle of Attack (deg)')
plt.ylabel('Lift Coefficient ($C_L$)')
plt.title('Lift Coefficient vs. Angle of Attack')
plt.grid(True)
plt.legend()
plt.show()

print(alo4)


# In[6]:


# Cmac
A1 = 0.162990
A2 = 0.0277226
Cmac = pi/4 * (A2 - A1)
Cmac.evalf()


# In[7]:


#A values for airfoil 1
from sympy import *
import numpy as np
x_prime, t, n = symbols('x_prime, t, n')

#m = 0.03
#p = 0.4
eq1 = (0.03/ 0.4**2) * (2 * 0.4 * x_prime - x_prime**2)

eq2 = (0.03/ (1- 0.4)**2) * ( 1 - 2*0.4+ 2*0.4*x_prime - x_prime**2)

EQ1 = diff(eq1, x_prime)
EQ2 = diff(eq2, x_prime)


integral1 = integrate(EQ1.subs(x_prime, (0.5 * (1 - cos(t))))*(cos(n*t)), (t, 0, acos(1 - 2 *0.4)))

integral2 = integrate(EQ2.subs(x_prime, (0.5 * (1 - cos(t))))*(cos(n*t)), (t, acos(1- 2*0.4), pi))

A1 = (2/pi) * (integral1.subs(n,1) + integral2.subs(n,1))
A2 = (2/pi) * (integral1.subs(n,2) + integral2.subs(n,2))

Cmac1 = pi/4 * (A2 - A1)
Cmac1.evalf()


# In[8]:


#A values for airfoil 2
from sympy import *
import numpy as np
x_prime, t, n = symbols('x_prime, t, n')

#m = 0.03
#p = 0.3
eq1 = (0.03/ 0.3**2) * (2 * 0.3 * x_prime - x_prime**2)

eq2 = (0.03/ (1- 0.3)**2) * ( 1 - 2*0.3+ 2*0.3*x_prime - x_prime**2)

EQ1 = diff(eq1, x_prime)
EQ2 = diff(eq2, x_prime)


integral1 = integrate(EQ1.subs(x_prime, (0.5 * (1 - cos(t))))*(cos(n*t)), (t, 0, acos(1 - 2 *0.3)))

integral2 = integrate(EQ2.subs(x_prime, (0.5 * (1 - cos(t))))*(cos(n*t)), (t, acos(1- 2*0.3), pi))

A1 = (2/pi) * (integral1.subs(n,1) + integral2.subs(n,1))
A2 = (2/pi) * (integral1.subs(n,2) + integral2.subs(n,2))

Cmac2 = pi/4 * (A2 - A1)
Cmac2.evalf()


# In[9]:


import numpy as np
import matplotlib.pyplot as plt

# Experimental data for Reynolds number 3.1e6
data_reynolds_3e6 = {
    'alpha': [-14, -12, -10, -8, 0, 8, 10, 12, 14, 16, 18, 20],
    'cl': [-0.7, -0.8, -0.7, -0.5, 0.4, 1.2, 1.3, 1.4, 1.5, 1.4, 1.3, 1.25]
}

# Experimental data for Reynolds number 5.7e6
data_reynolds_5e6 = {
    'alpha': [-14, -12, -10, -8, 0, 8, 10, 12, 14, 16, 18, 20],
    'cl': [-0.7, -0.8, -0.7, -0.5, 0.4, 1.2, 1.4, 1.6, 1.7, 1.6, 1.5, 1.3]
}

# Experimental data for Reynolds number 8.9e6
data_reynolds_9e6 = {
    'alpha': [-14, -12, -10, -8, 0, 8, 10, 12, 14, 16, 18, 20],
    'cl': [-0.7, -0.8, -0.7, -0.5, 0.4, 1.2, 1.4, 1.5, 1.6, 1.7, 1.7, 1.5]
}

# Experimental data for Reynolds number 6e9 Standard Roughness
data_reynolds_6e9 = {
    'alpha': [-14, -12, -10, -8, 0, 8, 10, 12, 14, 16],
    'cl': [-0.7, -0.8, -0.7, -0.5, 0.4, 1.15, 1.3, 1.4, 1.35, 1.25]
}


plt.figure(figsize=(14, 8))

plt.plot(data_reynolds_3e6['alpha'], data_reynolds_3e6['cl'], label='Reynolds Number 3.1e6', color='blue', marker='o', linestyle='dotted')
plt.plot(data_reynolds_5e6['alpha'], data_reynolds_5e6['cl'], label='Reynolds Number 5.7e6', color='red', marker='s', linestyle='dotted')
plt.plot(data_reynolds_9e6['alpha'], data_reynolds_9e6['cl'], label='Reynolds Number 8.9e6', color='green', marker='d', linestyle='dotted')
plt.plot(data_reynolds_6e9['alpha'], data_reynolds_6e9['cl'], label='Reynolds Number 5.7e9 Standard Roughness', color='orange', marker='^', linestyle='dotted')
plt.plot(np.rad2deg(alpha), cl4, label='TAT NACA 4412')
plt.xlabel('Angle of Attack (deg)')
plt.ylabel('Lift Coefficient ($C_L$)')
plt.title('Experimental Data for NACA 4412 Airfoil at Different Reynolds Numbers')
plt.grid(True)
plt.legend()
plt.show()


# In[10]:


#Reynolds numbers 

#assuming standard sea level conditions

rho = 1.225 #kg/m^3
v = 10 #m/s
l = 0.7 #m
mu = 1.789e-5 #kg/m/s

Re = rho * v * l / mu
Re


# In[11]:


#Reynolds numbers 

#assuming standard sea level conditions

rho = 1.225 #kg/m^3
v = 15 #m/s
l = 0.7 #m
mu = 1.789e-5 #kg/m/s

Re = rho * v * l / mu
Re


# In[12]:


#iii
import matplotlib.pyplot as plt

# Data for each Reynolds number
re_1e6 = {
    -5: -0.08045488,
    -4.5: -0.02460377,
    -4.279136: 1.279155E-16,
    -4: 0.03135513,
    -3.5: 0.0871898,
    -3: 0.1428386,
    -2.5: 0.1982197,
    -1.5: 0.3092319,
    -1: 0.3645405,
    -0.5: 0.4198483,
    0: 0.4745808,
    0.5: 0.5273312,
    1: 0.5742739,
    1.5: 0.6445718,
    2: 0.6985486,
    2.5: 0.7521058,
    3: 0.8061427,
    3.5: 0.8603995,
    4: 0.9143056,
    4.5: 0.9683901,
    5: 1.020651,
    5.5: 1.07386,
    6: 1.125563,
    6.5: 1.175025,
    7: 1.222905,
    7.5: 1.265063,
    8: 1.304998,
    8.5: 1.340643,
    9: 1.372703,
    9.5: 1.403653,
    10: 1.432877
}

re_3e6 = {
    -5: -0.08752828,
    -4.5: -0.03054719,
    -4: 0.02604101,
    -3.5: 0.08279049,
    -3: 0.1394263,
    -2: 0.2528879,
    -1.5: 0.3096269,
    -1: 0.3664806,
    -0.5: 0.4230609,
    0: 0.4794562,
    0.5: 0.5364795,
    1: 0.5928533,
    2: 0.7041269,
    2.5: 0.7572531,
    3: 0.8177858,
    3.5: 0.8739894,
    4: 0.9285115,
    4.5: 0.9823137,
    5: 1.035165,
    5.5: 1.086014,
    6: 1.136507,
    6.5: 1.18322,
    7: 1.229648,
    7.5: 1.273923,
    8: 1.318661,
    8.5: 1.364223,
    9: 1.408925,
    9.5: 1.452796,
    10: 1.493554
}

re_6e6 = {
    -5: -0.09057773,
    -4.5: -0.03342218,
    -4: 0.02378626,
    -3.5: 0.08099343,
    -3: 0.1384873,
    -2.5: 0.1958374,
    -2: 0.253154,
    -1.5: 0.3101549,
    -1: 0.367865,
    -0.5: 0.4251601,
    0: 0.4826249,
    0.5: 0.5399974,
    1: 0.5975494,
    1.5: 0.6544252,
    2: 0.7119825,
    2.5: 0.7678224,
    3: 0.8240888,
    3.5: 0.8782773,
    4: 0.926019,
    4.5: 0.9847641,
    5: 1.036561,
    5.5: 1.087346,
    6: 1.136982,
    6.5: 1.186577,
    7: 1.23443,
    7.5: 1.2836,
    8: 1.332944,
    8.5: 1.381447,
    9: 1.428556,
    9.5: 1.475993,
    10: 1.521032
}

re_9e6 = {
    -5: -0.09233561,
    -4.5: -0.03470476,
    -4: 0.0228367,
    -3.5: 0.08050562,
    -3: 0.1378788,
    -2.5: 0.195326,
    -2: 0.2532258,
    -1.5: 0.3108144,
    -1: 0.3686301,
    -0.5: 0.4263497,
    0: 0.4843076,
    0.5: 0.541671,
    1: 0.5997366,
    1.5: 0.6570536,
    2: 0.7142495,
    2.5: 0.7704882,
    3: 0.8258907,
    3.5: 0.8805864,
    4: 0.9340116,
    4.5: 0.9862219,
    5: 1.036842,
    5.5: 1.087341,
    6.5: 1.189432,
    7: 1.239572,
    7.5: 1.290603,
    8: 1.340924,
    8.5: 1.390184,
    9: 1.43957,
    9.5: 1.487618,
    10: 1.535667
}

plt.figure(figsize=(20, 15))


plt.plot(list(re_1e6.keys()), list(re_1e6.values()), marker='o', label='Re = 1e6')
plt.plot(list(re_3e6.keys()), list(re_3e6.values()), marker='s', label='Re = 3e6')
plt.plot(list(re_6e6.keys()), list(re_6e6.values()), marker='^', label='Re = 6e6')
plt.plot(list(re_9e6.keys()), list(re_9e6.values()), marker='x', label='Re = 9e6')
plt.plot(np.rad2deg(alpha), cl4, label='TAT NACA 4412')

plt.xlabel('Alpha (degrees)')
plt.ylabel('Cl')
plt.title('Cl vs Alpha for Different Reynolds Numbers')
plt.grid(True)
plt.legend()


plt.show()


# In[ ]:





# In[ ]:




