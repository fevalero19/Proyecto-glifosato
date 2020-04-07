import numpy as np
from scipy.integrate import odeint
import scipy as sp
import scipy.optimize
import matplotlib.pyplot as plt

####################################################################################################################################
######################################################## Glifosato ############################
###################################################################################################

# proyecto final 

import numpy as np
from scipy.integrate import odeint
import scipy as sp
import scipy.optimize
import matplotlib.pyplot as plt
import sympy as sym 

# Las unidades temporales son días

# def modelroots(y):
#     A, T, W, P=y
#     M = 5.112     # Masa total inicial de glisfosato
#     phi = 0.00548 # se fumiga dos veces al año: 2/365
#     pi = 0.7      # proporción del herbicida que se va al aire. 1-pi es lo que va a la tierra
#     ga = 1/5      # 1/ga = lo que dura el glifosato en el aire: 5 días
#     gt = 1/20     # 1/gt = lo que dura el glifosato en la tierra: 20-60 días
#     G = 0.002116  # Cantidad en masa de glifosato que consume una persona al día. 2 L/día * 1058x10^-3 g/L
#     k = 0.8       # Proporción del glifosato en agua que se va al aire. 1-k se va directamente a la persona.
#     mw = 1/5      # 1/mw = Lo que dura el glifosato en el agua: 5 días
#     mp = 1/30      # 1/mp = Lo que dura el glifosato en el cuerpo de una persona: 3 días. 
#     n1 = 0.8 
#     Adot = M*phi*pi + k * mw * W - ga * n1 * A - ga * (1-n1) * A
#     Tdot = M * phi * (1-pi) * gt * T 
#     Wdot = gt * T - mw * k * W - G * (1-k)
#     Pdot = G * (1-k) + ga * n1 * A - mp * P
#     dydt=[Adot,Tdot,Wdot,Pdot]
#     return dydt

# # Condiciones iniciales
# y0 = [0.01,0.01,0.01,0.01]

# # Calculando equilibrios


# eqs = sp.optimize.root(modelroots,y0,method='hybr')

# print('Con los mismos parámetros del modelo y método hybr ',eqs.x)

# # Evaluando el comportamiento de las poblaciones

# M = 5.112     # Masa total inicial de glisfosato. Para 3 L de herbicida: 5.112 g
# phi = 0.00548 # se fumiga dos veces al año: 2/365
# pi = 0.7      # proporción del herbicida que se va al aire. 1-pi es lo que va a la tierra
# ga = 1/5      # 1/ga = lo que dura el glifosato en el aire: 5 días
# gt = 1/20     # 1/gt = lo que dura el glifosato en la tierra: 20-60 días
# G = 0.002116  # Cantidad en masa de glifosato que consume una persona al día. 2 L/día * 1058x10^-3 g/L
# k = 0.8       # Proporción del glifosato en agua que se va al aire. 1-k se va directamente a la persona.
# mw = 1/5      # 1/mw = Lo que dura el glifosato en el agua: 5 días
# mp = 1/3      # 1/mp = Lo que dura el glifosato en el cuerpo de una persona: 3 días. 
# n1 = 0.8      # Proporción del glifosato que va del aire a la persona. (1-n1) es la proporción que se sale del domo.


# def modelGlif(y,t,M,phi,pi,ga,gt,G,k,mw,mp,n1):
#     A, T, W, P=y
#     Adot = M*phi*pi + k * mw * W - ga * n1 * A - ga * (1-n1) * A
#     Tdot = M * phi * (1-pi) * gt * T 
#     Wdot = gt * T - mw * k * W - G * (1-k)
#     Pdot = G * (1-k) + ga * n1 * A - mp * P
#     dydt=[Adot,Tdot,Wdot,Pdot]
#     return dydt

# y01 = [M*(pi),M*(1-pi),0,0]
# t = np.linspace(0,100,100)

# sol=odeint(modelGlif,y01,t,args=(M,phi,pi,ga,gt,G,k,mw,mp,n1))

# plt.plot(t,sol[:,0],'g',label='A(t)')
# plt.plot(t,sol[:,1],'b',label='T(t)')
# plt.plot(t,sol[:,2],'r',label='P(t)')
# plt.plot(t,sol[:,3],'p',label='W(t)')
# plt.legend(loc='best')
# plt.xlabel('t')
# plt.ylabel('Gramos de glifosato')
# plt.title('Comportamiento del modelo')
# plt.show()


####################################################################################################################################
######################################################## ADYUVANTES ############################
###################################################################################################

# Evaluando el comportamiento de las poblaciones

print('Variando la masa total de adyuvante m')

m = 9.72     # Masa total inicial de adyuvantes. En 3 L de herbicida * (1,02 Kg/L) = 3,6 Kg y 2,7% w/w de herbiicida --> 
             # m = 3,6 Kg * (0,0027) = 9,72 kg
phi = 0.00548  # se fumiga dos veces al año: 2/365
pi = 0.7      # proporción del herbicida que se va al aire. 1-pi es lo que va a la tierra
ga = 1/5      # 1/ga = lo que dura el adyuvante en el aire: 5 días
gt = 1/20     # 1/gt = lo que dura el glifosato en la tierra: 20-60 días
k = 0.8       # Proporción del glifosato en agua que se va al aire. 1-k se va directamente a la persona.
G = m * (1-pi) * (1-k) * 2 / 160   # Cantidad en masa de glifosato que consume una persona al día. 2 L/día * 1058x10^-3 g/L
# Se va a asumir que todo el adyuvante queda en la capa superior
# del lago --> atraviesa un área de 1 m^2
# D= 9m^2/s * 1/1 m^2

mw = 1/5    # 1/mw = Lo que dura el glifosato en el agua: 5 días --> difusividad de agua

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

################################################################################
################################################################################
################################################################################
print('Variando la difusividad aire-agua del surfactante')

# Se va a asumir que todo el adyuvante queda en la capa superior
# del lago --> atraviesa un área de 1 m^2
# D= 9m^2/s * 1/1 m^2
mw = 0.000104     

# Primero con m --> el herbicida con adyuvante al 0.27% w/w
sol=odeint(modelAdyuv,y01,t,args=(m,phi,pi,ga,gt,G,k,mw,mp,n1))

plt.plot(t,sol[:,0],'g',label='A(t)')
plt.plot(t,sol[:,1],'b',label='T(t)')
plt.plot(t,sol[:,2],'r',label='P(t)')
plt.plot(t,sol[:,3],'p',label='W(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylabel('Gramos de adyuvante')
plt.title('Comportamiento del POE-15: 0.27% w/w, Dw = 0.000104')
plt.show()

################################################################################
# Luego con mw = 0.000104 y m11 --> Herbicida al 11/07% w/w de adyuvante

sol11=odeint(modelAdyuv,y01,t11,args=(m11,phi,pi,ga,gt,G,k,mw,mp,n1))

plt.plot(t11,sol11[:,0],'g',label='A(t)')
plt.plot(t11,sol11[:,1],'b',label='T(t)')
plt.plot(t11,sol11[:,2],'r',label='P(t)')
plt.plot(t11,sol11[:,3],'p',label='W(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylabel('Gramos de adyuvante')
plt.title('Comportamiento del POE-15: 11.07% w/w, tiempo extendido, Daw = 0.000104 m2/s')
plt.show()



################################################################################

# Ahora variamos por primera vez la mw y probamos para ambos casos de masa


# Se va a asumir que todo el adyuvante queda en la capa superior
# del lago --> atraviesa un área de 1 m^2
# D= 9m^2/s * 1/1 m^2 --> y (60*60=3600) seg / 1 día --> unidades finales 1/día
mw = 0.0000324     

sol=odeint(modelAdyuv,y01,t,args=(m,phi,pi,ga,gt,G,k,mw,mp,n1))

plt.plot(t,sol[:,0],'g',label='A(t)')
plt.plot(t,sol[:,1],'b',label='T(t)')
plt.plot(t,sol[:,2],'r',label='P(t)')
plt.plot(t,sol[:,3],'p',label='W(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylabel('Gramos de adyuvante')
plt.title('Comportamiento del POE-15: 0.27% w/w , Daw = 0.0000234 m2/s')
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
plt.title('Comportamiento del POE-15: 11.07% w/w, tiempo extendido, Daw = 2.8 m2/s')
plt.show()

################################################################################
# Gráficas 3d
################################################################################

difusividades = np.linspace(0.0000324,0.00548,100)
y01 = [m*(pi),m*(1-pi),0,0]
t11 = np.linspace(0,1000,100)
m = 9.72
phi = 0.00548  # se fumiga dos veces al año: 2/365
pi = 0.7      # proporción del herbicida que se va al aire. 1-pi es lo que va a la tierra
ga = 1/5      # 1/ga = lo que dura el adyuvante en el aire: 5 días
gt = 1/20     # 1/gt = lo que dura el glifosato en la tierra: 20-60 días
k = 0.8       # Proporción del glifosato en agua que se va al aire. 1-k se va directamente a la persona.
G = m * (1-pi) * (1-k) * 2 / 160   # Cantidad en masa de glifosato que consume una persona al día. 2 L/día * 1058x10^-3 g/L

# Se va a asumir que todo el adyuvante queda en la capa superior
# del lago --> atraviesa un área de 1 m^2
# D= 9m^2/s * 1/1 m^2 --> y 1 s / (60*60=3600) día --> unidades finales 1/día
contador =0
mp = 1/3      # 1/mp = Lo que dura el glifosato en el cuerpo de una persona: 3 días.
n1 = 0.8  
solis = np.zeros((100,100))
# pgt = np.zeros((50,50))
for i  in range(len(difusividades)):
    mw = difusividades[i]
    lasol = odeint(modelAdyuv,y01,t11,args=(m,phi,pi,ga,gt,G,k,mw,mp,n1))
    tempo = lasol[:,2] # personas
    # pgt[:,contador]=lasol[:,2]
    # tempo = np.expand_dims(tempo,axis=1)
    solis[:,contador]=tempo
    contador = contador + 1



tv,gtv=np.meshgrid(t11,difusividades)

fig = plt.figure()
ax = fig.gca(projection='3d')
z=np.linspace(0,40)
surf = ax.plot_surface(tv,gtv, solis, cmap=cm.coolwarm, linewidth=0)

fig.colorbar(surf, shrink=0.5, aspect=5)

ax.set_xlabel('t')
ax.set_ylabel('1/mw: Difusividad del adyuvante en agua')
ax.set_zlabel('Gramos de adyuvante en persona')

ax.view_init(40,40)

plt.show()