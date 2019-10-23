import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from tkinter import Button, Tk, Entry, Label, LabelFrame, Scale, HORIZONTAL, ttk, DISABLED
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.integrate import odeint
import math
import time
from matplotlib.animation import FuncAnimation

class Point:
    def __init__(self, x = 0, y = 0, u = 0, v = 0, m = 1, color = 'red', time = 1):
        self.x = x
        self.y = y
        self.u = u
        self.v = v
        self.m = m
        self.color = color
        self.time = time

class Emitter:
    def __init__(self, x = 0, y = 0):
        self.x = x
        self.y = y

    def GeneratePoint(self, u = 1, v = 1, m = 1, color = 'red', time = 1 ):
        point = Point(x = self.x, y = self.y, m = m, color = color, time = time)
        return point


point_list = []
emitter = Emitter()


def PointsScatter(point_list):
    # fig, ax = plt.subplots()
    N = len(point_list)
    x = np.zeros([N])
    y = np.zeros([N])
    s = np.zeros([N])
    c = []
    for i, point in enumerate(point_list):
        x[i] = point.x
        y[i] = point.y
        s[i] = point.m
        c.append(point.color)
    scat = plt.scatter(x, y, s, c)
    # plt.show()
    return scat


def PointsScatterSolarSystem(point_list):
    # fig, ax = plt.subplots()
    N = len(point_list)
    x = np.zeros([N])
    y = np.zeros([N])
    s = np.zeros([N])
    c = []
    for i, point in enumerate(point_list):
        x[i] = point.x
        y[i] = point.y
        s[i] = 1
        c.append(point.color)
    scat = plt.scatter(x, y, s, c)
    # plt.show()
    return scat

def UpdatePointList(point_list, y_array, delta_t):
    del_indexes = []
    for i, point in enumerate(point_list):
        point.x = y_array[4*i]
        point.y = y_array[4*i + 1]
        point.u = y_array[4*i + 2]
        point.v = y_array[4*i + 3]
        point.t -= delta_t
        if point.t <= 0:
            del_indexes.append[i]

    for index in del_indexes:
        point_list.pop(index)
    return point_list

def GenerateInitialData(point_list):
    N = len(point_list)
    y0 = np.zeros([4*N])
    m_array = np.zeros([N])
    for i, point in enumerate(point_list):
        y0[4*i] = point.x
        y0[4 * i + 1] = point.y
        y0[4 * i + 2] = point.u
        y0[4 * i + 3] = point.v
        m_array[i] = point.m
    return y0, m_array

def RightPartFunc(y, t, m_array):
    N = y.size//4
    f = np.zeros_like(y)
    G = 6.67e-11
    for i in range(0,N):
        f[4*i] = y[4*i+2]
        f[4*i + 1] = y[4*i + 3]
        for j in range(0, N):
            if not j==i:
                z = math.pow((y[4*i] - y[4*j])**2 + (y[4*i+1] - y[4*j+1])**2, 1.5)
                f[4*i + 2] += (y[4*j] - y[4*i]) * G * m_array[j] / z / m_array[i]
                f[4*i + 3] += (y[4*j+1] - y[4*i+1]) * G * m_array[j] / z / m_array[i]
    return f

def Jacobian (y,t,m_array):
    N = y.size//4
    G = 6.67e-11
    f = np.zeros([4*N, 4*N])
    for i in range(0, N):
        f[4*i,4*i+2] = 1
        f[4*i+1,4*i+3] = 1
        for j in range (0,N):
            if not j==i:
                z = math.pow((y[4 * i] - y[4 * j]) ** 2 + (y[4 * i + 1] - y[4 * j + 1]) ** 2, 1.5)
                f[4 * i + 2, 4 * j] = G * m_array[j]/z
                f[4 * i + 3, 4 * j + 1] = G * m_array[j] / z
                f[4 * i + 2, 4 * i] -= G * m_array[j] / z
                f[4 * i + 3, 4 * i + 1] -= G * m_array[j] / z
    return f

def CalculateAccelerations (y, m_array):
    N = y.size//4
    a = np.zeros([2*N])
    G = 6.67e-11
    for i in range(0, N):
        for j in range(0, N):
            if not j==i:
                z = math.pow((y[4*i] - y[4*j])**2 + (y[4*i+1] - y[4*j+1])**2, 1.5)
                a[2*i] += (y[4*j] - y[4*i]) * G * m_array[j] / z / m_array[i]
                a[2*i + 1] += (y[4*j+1] - y[4*i+1]) * G * m_array[j] / z / m_array[i]
    return a

def VerletSolve (y0, t_array, m_array):
    N = y0.size//4
    T = len(t_array)
    y = np.zeros([T, 4*N]) #вывод такого же формата как в odeint
    y[0] = y0
    for index, t in enumerate(t_array):
        if t > 0:
            y_pred = y[index - 1]
            t_pred = t_array[index - 1]
            a_pred = CalculateAccelerations(y_pred, m_array)
            print(a_pred[2 * 0], a_pred[0 + 1])
            delta_t = t - t_pred
            y_cur = np.zeros([4*N])
            for i in range (0,N):
                y_cur[4*i] = y_pred[4*i] + y_pred[4*i + 2] * delta_t + 0.5 * a_pred[2*i] * delta_t**2
                y_cur[4*i + 1] = y_pred[4*i + 1] + y_pred[4*i + 3] * delta_t + 0.5 * a_pred[2*i + 1] * delta_t**2
            a_cur = CalculateAccelerations(y_cur, m_array)
            for i in range(0,N):
                y_cur[4*i + 2] = y_pred[4*i + 2] + 0.5 * (a_pred[2*i] + a_cur[2*i]) * delta_t
                y_cur[4*i + 3] = y_pred[4*i + 3] + 0.5 * (a_pred[2*i + 1] + a_cur[2*i + 1]) * delta_t
            y[index] = y_cur
    return y


def VerletNextLayer(point_list, delta_t):
    y0, m_array = GenerateInitialData(point_list)
    N = y0.size//4
    y = np.zeros_like(y0)

    a_pred = CalculateAccelerations(y0, m_array)
    for i in range (0,N):
        y[4*i] = y0[4*i] + y0[4*i + 2] * delta_t + 0.5 * a_pred[2*i] * delta_t**2
        y[4*i + 1] = y0[4*i + 1] + y0[4*i + 3] * delta_t + 0.5 * a_pred[2*i + 1] * delta_t**2

    a_cur = CalculateAccelerations(y, m_array)
    for i in range(0,N):
        y[4*i + 2] = y0[4*i + 2] + 0.5 * (a_pred[2*i] + a_cur[2*i]) * delta_t
        y[4*i + 3] = y0[4*i + 3] + 0.5 * (a_pred[2*i + 1] + a_cur[2*i + 1]) * delta_t
    point_list = UpdatePointList(point_list, y, delta_t)
    print("next layer")
    return point_list


def OdeintNextLayer(point_list, delta_t):
    y0, m_array = GenerateInitialData(point_list)
    t_array = np.ndarray([0, delta_t])
    y = odeint(RightPartFunc, y0, t_array, (m_array,), Jacobian)[1]
    point_list = UpdatePointList(point_list, y, delta_t)
    return point_list


def GeneratePointsFromFile(path):
    point_list = []
    with open(path, 'r') as fileread:
        lines = fileread.readlines()
    for line in lines:
        point_info = line.split()
        x = float(point_info[0])
        y = float(point_info[1])
        u = float(point_info[2])
        v = float(point_info[3])
        m = float(point_info[4])
        c = point_info[5]
        t = float(point_info[6])
        point = Point(x, y, u, v, m, c, t)
        point_list.append(point)
    return point_list



def LoadButtonClicked():
    point_list = GeneratePointsFromFile("C:/Users/user/Desktop/SolarSystemInfo.txt")
    fig = plt.figure()
    graph_canvas = FigureCanvasTkAgg(fig, master=root)
    graph_canvas.get_tk_widget().grid(column = 2)

    ax = fig.add_axes()
    scat = PointsScatterSolarSystem(point_list)
    # plt.show()
    if combo.get() == "Verlet":
        nextlayer_func = VerletNextLayer
    else:
        nextlayer_func = OdeintNextLayer
    FuncAnimation(fig, lambda x: nextlayer_func(point_list, t_delta=1e10), interval=500)

def EmitterButtonClicked():
    x = float(x_entry.get())
    y = float(y_entry.get())
    emitter.x = x
    emitter.y = y
    load_button.configure(state=DISABLED)
    plt.scatter(x, y)
    print("new emitter created")
    return






root = Tk()
root.title("Гравитационная задача N тел")

emitter_info_label = Label(bg='white', fg='black', width=40, text='Emitter info:')
emitter_info_label.grid(column=0, row=0)

x_label = Label(bg='white', fg='black', width=40, text='Emitter x coordinate:')
x_label.grid(column=0, row=1)
x_entry = Entry(width=10, text = "0.0")
x_entry.grid(column=0, row=2)
y_label = Label(bg='white', fg='black', width=40, text='Emitter y coordinate:')
y_label.grid(column=0, row=3)
y_entry = Entry(width=10, text = "0.0")
y_entry.grid(column=0, row=4)

create_emitter_button = Button(text="Create/Update emitter", command = EmitterButtonClicked)
create_emitter_button.grid(column=0, row=5)


point_info_label = Label(root, bg='white', fg='black', width=40, text='Point info:')
point_info_label.grid(column=1, row=0)
u_label = Label(root, bg='white', fg='black', width=40, text='point u speed:')
u_label.grid(column=1, row=1)
u_entry = Entry(root, width=10, text = "1.0")
u_entry.grid(column=1, row=2)
v_label = Label(root, bg='white', fg='black', width=40, text='point v speed:')
v_label.grid(column=1, row=3)
v_entry = Entry(root, width=10, text = "1.0")
v_entry.grid(column=1, row=4)

m_label = Label(root, bg='white', fg='black', width=40, text='point mass:')
m_label.grid(column=1, row=5)
m_scale = Scale(root, orient=HORIZONTAL, from_=0, to=10)
m_scale.set(5)
m_scale.grid(column=1, row=6)
t_label = Label(root, bg='white', fg='black', width=40, text='point lifetime:')
t_label.grid(column=1, row=7)
t_entry = Entry(root, width=10, text = "100")
t_entry.grid(column=1, row=8)
color_label = Label(root, bg='white', fg='black', width=40, text='point color:')
color_label.grid(column=1, row=7)
color_entry = Entry(root, width=10, text = "green")
color_entry.grid(column=1, row=8)

create_point_button = Button(root, text='Create point')
create_point_button.grid(column=1, row=9)


method_label = Label(root, bg='white', fg='black', width=40, text='Choose method:')
method_label.grid(column=0, row=7)
combo = ttk.Combobox(root)
combo['values'] = ("odeint", "Verlet")
combo.current(1)
combo.grid(column=0, row=8)


load_button = Button(root, text="Load point_list from file", command = LoadButtonClicked)
load_button.grid(column = 0, row = 9)


root.geometry('1000x1200')
root.mainloop()
