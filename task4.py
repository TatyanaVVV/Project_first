from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

#Create mesh and define function space
R = 1.0
circle_r = Circle(Point(0, 0), R) #задача в круге
mesh = generate_mesh(circle_r, 64)
V = FunctionSpace(mesh, 'P', 1)
bounds = MeshFunction("size_t", mesh, 1)

#-------------------ГРАНИЧНЫЕ УСЛОВИЯ----------------------------------
class boundary1(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary or (x[1] >= 0)
b1 = boundary1()
b1.mark(bounds, 0)
class boundary2(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary or (x[1] < 0)
b2 = boundary2()
b2.mark(bounds, 1)

alpha = 1
_D = Expression('-2*sin(x[0])*cos(x[1]) - 2*sin(x[1])*cos(x[0])', degree = 2) 
f = Expression('sin(x[1])', degree = 2, alpha = alpha)
g =Expression('cos(x[1])', degree = 2, R = R)


#ВАРИАЦИОННАЯ ЗАДАЧА
u = TrialFunction(V)
v = TestFunction(V)
bc = DirichletBC(V, u_D, bounds, 1)
a = dot(grad(u), grad(v))*dx + alpha*dot(u,v)*dx
L = f*v*dx - g*v*ds(0, subdomain_data = bounds)
u = Function(V)
solve(a == L, u, bc)

#НОРМЫ L2 и максимум-норма
error_L2 = errornorm(u_D, u, 'L2')
vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
error_C = np.max(np.abs(vertex_values_u - vertex_values_u_D))
print('L2-error = ', error_L2)
print('C-error = ', error_C)

#Визуализация с помощью matplotlib
n = mesh.num_vertices()
d = mesh.geometry().dim()
mesh_coordinates = mesh.coordinates().reshape((n, d))
triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
triangulation = tri.Triangulation(mesh_coordinates[:, 0], mesh_coordinates[:, 1], triangles)
plt.figure()
zfaces = np.asarray([u_D(cell.midpoint()) for cell in cells(mesh)])
plt.tripcolor(triangulation, facecolors=zfaces, edgecolors='k')
plt.colorbar()
plt.title("Аналитическое решение")
plt.savefig('аналитическое.png')

plt.figure()
zfaces = np.asarray([u(cell.midpoint()) for cell in cells(mesh)])
plt.tripcolor(triangulation, facecolors=zfaces, edgecolors='k')
plt.colorbar()
plt.title("Численное решение")
plt.savefig('численное.png')