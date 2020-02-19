import numpy as np
from scipy.integrate import odeint
import scipy as sp
import scipy.optimize
import matplotlib.pyplot as plt

####################################################################################################################################
######################################################## ADYUVANTES ############################
###################################################################################################

# Evaluando el comportamiento de las poblaciones

print('Variando la masa total de adyuvante m')

m = 9.72     # Masa total inicial de glisfosato. Para 3 L de herbicida: 5.112 g
phi = 0.00548  # se fumiga dos veces al año: 2/365
pi = 0.7      # proporción del herbicida que se va al aire. 1-pi es lo que va a la tierra
ga = 1/5      # 1/ga = lo que dura el adyuvante en el aire: 5 días
gt = 1/20     # 1/gt = lo que dura el glifosato en la tierra: 20-60 días
k = 0.8       # Proporción del glifosato en agua que se va al aire. 1-k se va directamente a la persona.
G = m * (1-pi) * (1-k) * 2 / 160   # Cantidad en masa de glifosato que consume una persona al día. 2 L/día * 1058x10^-3 g/L
mw = 1/5      # 1/mw = Lo que dura el glifosato en el agua: 5 días
mp = 1/3      # 1/mp = Lo que dura el glifosato en el cuerpo de una persona: 3 días.
n1 = 0.8      # Proporción del glifosato que va del aire a la persona. (1-n1) es la proporción que se sale del domo.


def modelAdyuv(y,t,m,phi,pi,ga,gt,G,k,mw,mp,n1):
    A, T, W, P=y
    Adot = m*phi*pi + k * mw * W - ga * n1 * A - ga * (1-n1) * A
    Tdot = m * phi * (1-pi) - gt * T
    Wdot = gt * T - mw * k * W - G * (1-k)
    Pdot = G * (1-k) + ga * n1 * A - mp * P
    dydt=[Adot,Tdot,Wdot,Pdot]
    return dydt

y01 = [m*(pi),m*(1-pi),0,0]
t = np.linspace(0,100,100)

sol=odeint(modelAdyuv,y01,t,args=(m,phi,pi,ga,gt,G,k,mw,mp,n1))

plt.plot(t,sol[:,0],'g',label='A(t)')
plt.plot(t,sol[:,1],'b',label='T(t)')
plt.plot(t,sol[:,2],'r',label='P(t)')
plt.plot(t,sol[:,3],'p',label='W(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylabel('Gramos de adyuvante')
plt.title('Comportamiento del POE-15: 0.27% w/w')
plt.show()

m11 = 398.52  # Gramos de POE-15 en 3L de herbicida al 11.07% w/w
# Esa masa la saque con la densidad del POE-15 que saqué de:
# https://www.stepan.com/products/Surfactants/TOXIMUL%C2%AE/TOXIMUL%C2%AE-TA-15.aspx

t11 = np.linspace(0,100,100)
y11 = [m*(pi),m*(1-pi),0,0]

sol11=odeint(modelAdyuv,y01,t11,args=(m11,phi,pi,ga,gt,G,k,mw,mp,n1))

plt.plot(t11,sol11[:,0],'g',label='A(t)')
plt.plot(t11,sol11[:,1],'b',label='T(t)')
plt.plot(t11,sol11[:,2],'r',label='P(t)')
plt.plot(t11,sol11[:,3],'p',label='W(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylabel('Gramos de adyuvante')
plt.title('Comportamiento del POE-15: 11.07% w/w, tiempo extendido')
plt.show()