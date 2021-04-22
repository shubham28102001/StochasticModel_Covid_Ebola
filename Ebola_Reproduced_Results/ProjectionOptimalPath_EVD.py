from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = plt.axes(projection="3d")


mu=0.00005
Bi=0.5
Bd=0.6
Bh=0.00016
N=1000
sigma=0.1
gir=0.07
mue=0.12
tau=0.2
delta=0.33
ghr=0.1
k=0

S = []
E = []
I = []
R = []
D = []
H = []

S.append(6200000)
E.append(5000000)
I.append(2500)
R.append(0)
D.append(0)
H.append(0)

i = 0
while (i < 10000):
    D.append(D[i] + ((mue*I[i]-((delta+mu)*D[i]))/N))
    H.append(H[i] + ((tau*I[i]-((ghr+mue+mu)*H[i]))/N))
    R.append(R[i] + ((gir*I[i]+ghr*H[i]-mu*R[i])/N))
    I.append(I[i] + ((sigma*(E[i]/N))-((gir+mue+tau+mu)*(I[i]/N))))
    E.append(E[i] + ((((Bi*I[i]*S[i])+(Bd*D[i]*S[i])+(Bh*H[i]*S[i]))/N**2)-((mu+sigma)*(E[i]/N))+k*(S[i]/N)))
    S.append(S[i] + (mu-(((Bi*I[i]*S[i])+(Bd*D[i]*S[i])+(Bh*H[i]*S[i]))/N**2)-mu*(S[i]/N)-k*(S[i]/N)))
    i += 1
    
plt.title("Optimal Path to Extinction")
ax.plot3D(S, E, I, 'blue')
ax.set_ylabel('Exposed')
ax.set_xlabel('Susceptible ')
ax.set_zlabel('Infected')
plt.show()
