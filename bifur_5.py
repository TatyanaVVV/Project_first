import matplotlib.pyplot as plt
import numpy as np

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

k2=0.95
k3=0.0032
k3m=0.002
#k1=0.12
#k3m=0.004
#alk = k3m/(k3m + k3)
num=1000
x=np.linspace(0, 0.999, num)
y=np.zeros(num)
z=np.zeros(num)
k1=np.zeros(num)
k1m=np.zeros(num)
K1=np.zeros(num)
K1m=np.zeros(num)
sp=np.zeros(num)
DI=np.zeros(num)
delta_a=np.zeros(num)
# СЂР°СЃС‡РµС‚ Р»РёРЅРёРё РЅРµР№С‚СЂР°Р»СЊРЅРѕСЃС‚Рё
for i in range(num):
    y[i]=(1-x[i])*x[i]*k3/(k3m+2*k3*x[i])
    z[i] = (1-x[i]-2*y[i])
    k1m[i] = (k2*z[i]*z[i]*(x[i]-z[i])-k3*z[i]*(1+x[i]-2*y[i])+k3m*(3*y[i]+x[i]-1))/(1-2*y[i])
    k1[i]= (k1m[i]*x[i]+k3*x[i]*z[i]-k3m*y[i]+k2*x[i]*z[i]*z[i])/z[i]
    a11 = -k1[i]-k1m[i]-k2*z[i]*z[i]+2*k2*x[i]*z[i]-k3*(1-2*y[i]-2*x[i])
    a12 = -2*k1[i]+2*k3*x[i]+k3m+4*k2*x[i]*z[i]
    a21 = k3*(1-2*x[i]-2*y[i])
    a22 = -2*k3*x[i]-k3m
    sp[i]= a11+a22
    delta_a[i]=a11*a22-a12*a21
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(k1, k1m, label='Линия нейтральности', color='r')
# СЂР°СЃС‡РµС‚ Р»РёРЅРёРё РєСЂР°С‚РЅРѕСЃС‚Рё
for i in range(num):
    y[i]=(1-x[i])*x[i]*k3/(k3m+2*k3*x[i])
    z[i] = (1-x[i]-2*y[i])
    #K1m[i] = (k3m*k3m*y[i]+k3m*k3*z[i]*(2*y[i]-x[i])-k2*k3m*z[i]*z[i]*(z[i]-x[i]))/(4*k3*x[i]*z[i]+k3m*(x[i]+z[i]))
    K1m[i] = (k3m*k3m*y[i]+k3m*k3*z[i]*(2*y[i]-x[i])-k2*k3m*z[i]*z[i]*(z[i]-x[i])-2*k3*k3*x[i]*z[i]*z[i])/(4*k3*x[i]*z[i]+k3m*(x[i]+z[i]))
    K1[i]= (K1m[i]*x[i]+k3*x[i]*z[i]-k3m*y[i]+k2*x[i]*z[i]*z[i])/z[i]
    a11 = -K1[i]-K1m[i]-k2*z[i]*z[i]+2*k2*x[i]*z[i]-k3*(1-2*y[i]-2*x[i])
    a12 = -2*K1[i]+2*k3*x[i]+k3m+4*k2*x[i]*z[i]
    a21 = k3*(1-2*x[i]-2*y[i])
    a22 = -2*k3*x[i]-k3m
    sp[i]= a11+a22
    delta_a[i]=a11*a22-a12*a21
    if (abs(K1m[i]-k1m[i])<1e-6):
        ax1.plot(K1[i],K1m[i], 'r*', marker='o', label='точка бифуркации коразмерности 2 Такенса Богданова')
for i in range(1,num-1):
    if (K1m[i]>K1m[i-1]) and (K1m[i]>K1m[i+1]):
        ax1.plot(K1[i],K1m[i], 'g*', marker='o',label='точка бифуркации коразмерности 2 трехкратное равновесие')
ax1 = fig.add_subplot(111)
ax1.plot(K1, K1m, label='Линия кратности', color='b')
#ax1.plot(x, y, label='xy', color='b')
plt.xlabel(u'K1', fontdict=font)
plt.ylabel(u'K1m', fontdict=font)
plt.legend()
plt.xlim([0.0, 0.15])
plt.ylim([0.0, 0.02])

k2=0.95
k3=0.0032
k1m=0.01
#k1=0.12
k3m=0.004 
num=1000
x=np.linspace(0, 0.999, num)
yh=[]
xh=[]
k1h=[]
ysn=[]
xsn=[]
k1sn=[]
ydi=[]
xdi=[]
k1di=[]
j1 = 0
j2 = 0
j3 = 0
i=1
y[i]=(1-x[i])*x[i]*k3/(k3m+2*k3*x[i])
z[i] = (1-x[i]-2*y[i])
K1[i]= (k1m*x[i]+k3*x[i]*z[i]-k3m*y[i]+k2*x[i]*z[i]*z[i])/z[i]
a11 = -K1[i]-k1m-k2*z[i]*z[i]+2*k2*x[i]*z[i]-k3*(1-2*y[i]-2*x[i])
a12 = -2*K1[i]+2*k3*x[i]+k3m+4*k2*x[i]*z[i]
a21 = k3*(1-2*x[i]-2*y[i])
a22 = -2*k3*x[i]-k3m
sp[i]= a11+a22
delta_a[i]=a11*a22-a12*a21
DI[i] = sp[i]*sp[i]-4*delta_a[i]
fig2 = plt.figure()
for i in range(2,num):
    y[i]=((1-x[i])*x[i]*k3)/(k3m+2*k3*x[i])
    z[i] = 1-x[i]-2*y[i]
    K1[i]= (k1m*x[i]+k3*x[i]*z[i]-k3m*y[i]+k2*x[i]*z[i]*z[i])/z[i]
    a11 = -K1[i]-k1m-k2*z[i]*z[i]+2*k2*x[i]*z[i]-k3*(1-2*y[i]-2*x[i])
    a12 = -2*K1[i]+2*k3*x[i]+k3m+4*k2*x[i]*z[i]
    a21 = k3*(1-2*x[i]-2*y[i])
    a22 = -2*k3*x[i]-k3m
    sp[i]= a11+a22
    delta_a[i]=a11*a22-a12*a21
    DI[i] = sp[i]*sp[i]-4*delta_a[i]
    if (sp[i]*sp[i-1]<=0):
        yh.append(y[i])
        xh.append(x[i])
        k1h.append(K1[i])
        plt.plot(K1[i], x[i], 'go', label='hopf')
        plt.plot(K1[i], y[i], 'go', label='hopf')
    if (delta_a[i]*delta_a[i-1]<=0):
        ysn.append(y[i])
        xsn.append(x[i])
        k1sn.append(K1[i])
        plt.plot(K1[i], x[i], 'r*', label='saddle-node')
        plt.plot(K1[i], y[i], 'r*', label='saddle-node')
    if (DI[i]*DI[i-1]<=0):
        ydi.append(y[i])
        xdi.append(x[i])
        k1di.append(K1[i])
        tst=plt.plot(K1[i], x[i], 'b^', label='node')
        plt.plot(K1[i], y[i], 'b^', label='node')

k1x_plot = fig2.add_subplot(111)
k1x_plot.plot(K1, x, label='k1(x)')
k1y_plot = fig2.add_subplot(111)
k1y_plot.plot(K1, y, label='k1(y)')
plt.title('Однопараметрический анализ k3m=0.004', fontdict=font)
plt.xlabel('XY', fontdict=font)
plt.ylabel('k1', fontdict=font)
plt.xlim([0.0, 0.2])
plt.ylim([0.0, 1])
plt.legend()

plt.show()