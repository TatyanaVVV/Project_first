import numpy as np
import random
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import tkinter
from tkinter import Button
from scipy.integrate import odeint
import math
import time



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
    def __init__(self, x = 0, y = 0, u = 1, v = 1):
        self.x = x
        self.y = y
        self.u = u
        self.v = v

    def GeneratePoint(self, max_m = 1, max_distance = 1, max_v = 1, max_time = 1):
        color_list = ['red','blue','green','yellow','orange','black']
        color = color_list[random.randint(0,5)]
        m = max_m * random.random()
        d = max_distance * random.random()
        x = self.x + self.u*d
        y = self.y + self.v*d
        u = max_v * (2*random.random() - 1)
        v = max_v * (2*random.random() - 1)
        time = max_time * random.random()
#координаты новой точки
        point = Point(x, y, m = m, color = color, time = time)
        return point

    def GeneratePoints(self, max_m = 1, max_distance = 1, max_v = 1, max_time = 1, N = 10):
        point_list = []
        for i in range(0,N):
            point_list.append(self.GeneratePoint(max_m,max_distance,max_v,max_time))
        return point_list


def PointsVisualize(point_list):
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    ax.set_xlim(-5000, 5000)
    ax.set_ylim(-5000, 5000)
    for point in point_list:
        if (point.m > 1e7):
            width = math.pow(point.m, 0.33)
            height = math.pow(point.m, 0.33)
        else:
            width = math.pow(point.m * 1000, 0.33)
            height = math.pow(point.m * 1000, 0.33)

        ell = Ellipse((point.x, point.y), width= width, height=height, color=point.color)
        ax.add_artist(ell)
        ell.set_clip_box(ax.bbox)

    plt.show()
    return

'''
point1 = Point(1, 2, m=3, color='red')
point2 = Point(3,5,m=4,color='blue')
point3 = Point(-4,1,m=1,color='black')


# root = tkinter.Tk()
# root.mainloop()
emitter = Emitter(u=1,v=2)
button = Button(None, text="Press Me", command=lambda: PointsVisualize(emitter.GeneratePoints(max_m = 3, max_distance=10)))
button.pack()
button.mainloop()
'''
def GenerateSolarSystem ():
    point_list = []
    point_list.append(Point(x = 57.9, y = 0, m = 3.3, color = 'grey')) # mercury
    point_list.append(Point(x = 108.2*math.sin(2), y = 108.2*math.cos(2), m = 48.7, color = 'cyan')) #venus
    point_list.append(Point(x = 149.6*math.sin(1), y = 149.6*math.cos(1), m = 59.7, color = 'blue')) #earth
    point_list.append(Point(x = 227.9*math.sin(4), y = 227.9*math.cos(4), m = 6.42, color = 'red')) #mars
    point_list.append(Point(x = 778.3*math.sin(5.9), y = 778.3*math.cos(5.9), m = 19e3, color = 'orange')) #jupiter
    point_list.append(Point(x = 1427*math.sin(3), y = 1427*math.cos(3), m = 5.7e3, color = 'pink')) #saturn
    point_list.append(Point(x = 2869*math.sin(4.2), y = 2869*math.cos(4.2), m = 870, color = 'violet')) #uran
    point_list.append(Point(x = 4496*math.sin(5), y = 4496*math.cos(5), m = 1024, color = 'green')) #neptun
    point_list.append(Point(m = 2e7, color = 'yellow')) #sun
    return point_list

PointsVisualize(GenerateSolarSystem())

def GenerateInitialData (point_list):
    N = len(point_list)
    y0 = np.zeros([4*N])
    m_array = np.zeros([N])
    for i in range(0,N):
        point = point_list[i]
        y0[4*i] = point.x
        y0[4 * i + 1] = point.y
        y0[4 * i + 2] = point.u
        y0[4 * i + 3] = point.v
        m_array[i] = point.m
    return y0, m_array

def RightPartFunc(y,t,m_array):
    N = y.size//4
    f = np.zeros_like(y)
    G = 1
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
    G = 1
    f = np.zeros([4*N,4*N])
    for i in range (0,N):
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

def UpdatePointList (point_list, y_array):
    N = len(point_list)
    for i in range(0, N):
        point_list[i].x = y_array[4*i]
        point_list[i].y = y_array[4*i + 1]
        point_list[i].u = y_array[4*i + 2]
        point_list[i].v = y_array[4*i + 3]
    return point_list

t_array = range(0,100)
point_list = GenerateSolarSystem()
N = len(point_list)
y0, m_array = GenerateInitialData(point_list)
m = (m_array,)
y = odeint(RightPartFunc, y0, t_array, (m_array, ), Jacobian)
for t in t_array:
    point_list = UpdatePointList(point_list, y[t])
    PointsVisualize(point_list)
    time.sleep(1)
