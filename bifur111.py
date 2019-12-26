from sympy import Symbol, solve, lambdify, Matrix
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sympy
from matplotlib.legend_handler import HandlerLine2D

def TwoParam(k1val, k1mval, k2val, k3val, k3mval):
    solution = solve([eq1, eq2], k1, y)
    k1_sol = solution[0][0]
    y_sol = solution[0][1]
    A = Matrix([eq1, eq2])
    var_vector = Matrix([x, y])
    jacA = A.jacobian(var_vector)
    detA = jacA.det()
    traceA = jacA.trace()

    k1_trace_sol = solve(traceA.subs(y,y_sol),k1)[0]
    k1m_trace_sol = solve(k1_trace_sol - k1_sol,k1m)[0]
    true_k1_trace_sol = k1_sol.subs(k1m,k1m_trace_sol)
    k1m_trace_y = lambdify((x, k1, k1m, k2, k3, k3m),k1m_trace_sol)
    k1_trace_y = lambdify((x, k1, k1m, k2, k3, k3m),true_k1_trace_sol)

    k1_det_sol = solve(detA.subs(y,y_sol),k1)[0]
    k1m_det_sol = solve(k1_det_sol - k1_sol,k1m)[0]
    true_k1_det_sol = k1_sol.subs(k1m,k1m_det_sol)
    k1m_det_y = lambdify((x, k1, k1m, k2, k3, k3m),k1m_det_sol)
    k1_det_y = lambdify((x, k1, k1m, k2, k3, k3m),true_k1_det_sol)

    plt.plot(k1_trace_y(X, k1val, k1mval, k2val, k3val, k3mval), k1m_trace_y(X, k1val, k1mval, k2val, k3val, k3mval),label='Линия нейтральности', color='r')
    plt.plot(k1_det_y(X, k1val, k1mval, k2val, k3val, k3mval), k1m_det_y(X, k1val, k1mval, k2val, k3val, k3mval),label='Линия кратности', color='b')
    plt.xlabel('k1')
    plt.ylabel('k1m')
    plt.xlim([0.0, 0.15])
    plt.ylim([0.0, 0.02])
    plt.legend()
    plt.show()
    return


def solveSystem(init, k1val, k1mval, k2val, k3val, k3mval, dt, iterations):

    f1 = lambdify((x, y, k1, k1m, k2, k3, k3m), eq1)
    f2 = lambdify((x, y, k3, k3m), eq2)
    def rhs(xy, times):
        return [f1(xy[0], xy[1], k1val, k1mval, k2val, k3val, k3mval), f2(xy[0], xy[1], k3val, k3mval)]
    times = np.arange(iterations) * dt
    return odeint(rhs, init, times), times


def autocol():
    res, times = solveSystem([0.5, 0.25], 0.12, 0.01, 0.95, 0.0032, 0.002, 0.01, 1e6)
    ax = plt.subplot(211)
    plt.plot(times, res[:, 0], 'b')
    plt.title('Релаксационные колебания')
    plt.ylabel('x')
    plt.xlabel('t')
    plt.grid()
    ax1 = plt.subplot(212)
    plt.plot(times, res[:, 1], 'r')
    plt.xlabel('t')
    plt.ylabel('y')
    plt.grid()
    plt.show()
    return
def streamplot(k1val, k1mval, k2val, k3val, k3mval):
    f1 = lambdify((x, y, k1, k1m, k2, k3, k3m), eq1)
    f2 = lambdify((x, y, k3, k3m), eq2)
    Y, X = np.mgrid[0:0.5:5000j, 0:1:5000j]
    U = f1(X, Y, k1val, k1mval, k2val, k3val, k3mval)
    V = f2(X, Y, k3val, k3mval)
    velocity = np.sqrt(U*U+V*V)
    plt.streamplot(X, Y, U, V, density = [2, 2], color=velocity)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Фазовый портрет системы')
    plt.show()
    
    return

def analysis_k1(elem, cur_k, k1val, k1mval, k2val, k3val, k3mval):
    if cur_k == 0:
        k1mval = elem
    if cur_k == 1:
        k3mval = elem
    fl = 0
    fl_delta = 0
    fl_di = 0
    bifurcation_point = []
    S_A = lambdify((x,y,k1,k1m,k2,k3,k3m),traceA)
    delta_A = lambdify((x,y,k1,k1m,k2,k3,k3m),detA)
    for i in range(N):
        cur_y = y_func(X[i],k1val,k1mval,k3val,k3mval)
        cur_k1 = k1_func(X[i],cur_y,k1mval,k2val,k3mval,k3val)
        k1_f[i] = cur_k1
        y_f[i] = cur_y
        
        sa = S_A(X[i], cur_y, k1val,k1mval,k2val,k3mval,k3val)
        deltaa = delta_A(X[i],cur_y,k1val,k1mval,k2val,k3mval,k3val)
        di = sa*sa - 4*deltaa
        if sa < 0:
            if fl != 0 and fl == 1:
                #print(X[i], cur_y)
                bifurcation_point.append([X[i],cur_y])
                plt.plot(k1_f[i],X[i],'r', marker='o', label="hopf_node")
                plt.plot(k1_f[i],y_f[i],'r', marker='o')
            fl = -1
            #print ('меньше нуля',i, sa, cur_y, cur_k2, fl)
        else:
            if fl != 0 and fl == -1:
                #print(X[i], cur_y)
                bifurcation_point.append([X[i],cur_y])
                plt.plot(k1_f[i],X[i],'r', marker='o', label="hopf_node")
                plt.plot(k1_f[i],y_f[i],'r', marker='o')
            fl = 1
            #print('больше нуля',i, sa, cur_y, cur_k2, fl)

        if deltaa < 0:
            if fl_delta != 0 and fl_delta == 1:
                #print ('определитель меньше нуля',i, deltaa, cur_y, cur_k2, fl)
                plt.plot(k1_f[i],X[i],'k*', marker='o', label="saddle_node")
                plt.plot(k1_f[i],y_f[i],'k*', marker='o')
            fl_delta = -1
        else:
            if fl_delta != 0 and fl_delta == -1:
                #print('определитель больше нуля',i, deltaa, cur_y, cur_k2, fl)
                plt.plot(k1_f[i],X[i],'k*', marker='o', label="saddle_node")
                plt.plot(k1_f[i],y_f[i],'k*', marker='o')
            fl_delta = 1
    line1, = plt.plot(k1_f, X, 'b--', label="x")
    line2, = plt.plot(k1_f, y_f, 'k', label="y")

    plt.legend(handler_map={line1: HandlerLine2D(numpoints=4)}) 
    if cur_k == 0:
        plt.title('Однопараметрический анализ. Зависимость стационарных решений от параметра k1, k_1 = ' + str(k1mval))
    else:
        plt.title('Однопараметрический анализ. Зависимость стационарных решений от параметра k1, k_3 = ' + str(k3mval))
    plt.plot(k1_f,X, 'b--')
    plt.plot(k1_f,y_f,'k')
    plt.xlabel('k1')
    plt.ylabel('x,y')
    plt.show()

k1 = Symbol('k1', Positive=True)
k1m = Symbol('k1m', Positive=True)
k2 = Symbol('k2', Positive=True)
k3 = Symbol('k3', Positive=True)
k3m = Symbol('k3m', Positive=True)
x = Symbol("x", Positive=True)
y = Symbol("y", Positive=True)
eq1 = k1 * (1 - x - 2*y)- k1m*x - k3*x*(1-x-2*y)+k3m*y-k2*x*(1-x-2*y)*(1-x-2*y)
eq2 = k3*x*(1-x-2*y)-k3m*y

X = np.linspace(0.001,0.987,1000)
N = np.size(X)


#Однопараметрический анализ 2 - значения
k_1_range = [0.001, 0.005, 0.01, 0.015, 0.02]
k_3_range = [0.0005, 0.001, 0.002, 0.003, 0.004]
k1_f = np.zeros(N)
y_f = np.zeros(N)
solution = solve([eq1, eq2], k1, y)
k1_sol = solution[0][0]
y_sol = solution[0][1]
y_func = lambdify((x,k1m,k2,k3,k3m),y_sol)
k1_func = lambdify((x,y,k1m,k2,k3,k3m),k1_sol) 
A = Matrix([eq1, eq2])
var_vector = Matrix([x, y])
jacA = A.jacobian(var_vector)
detA = jacA.det()
traceA = jacA.trace()
#----------------------ОДНОПАРАМЕТРИЧЕСКИЙ АНАЛИЗ ПО K1--------------------------------
for elem in k_3_range:
   analysis_k1(elem,0,0.12, 0.01, 0.95, 0.0032, 0.002)

for elem in k_1_range:
  analysis_k1(elem,1,0.12, 0.01, 0.95, 0.0032, 0.002)

streamplot(0.12, 0.01, 0.95, 0.0032, 0.002)
#----------------------ДВУХПАРАМЕТРИЧЕСКИЙ АНАЛИЗ ПО K1,K2--------------------------------
TwoParam(0.12, 0.01, 0.95, 0.0032, 0.002)
#колебания
autocol()