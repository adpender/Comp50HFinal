from fenics import * 
import time

#Subdivisions of x,y, and t sections 
totalTime = 10.0 
interval = 50
dt = totalTime/interval 

nx = 1000 
ny = 25

#Point definitions of mesh
a = Point(0,0)
b = Point(2.5, 2.5)

#Mesh definition and function space definition 
mesh = RectangleMesh(a,b,nx, ny)
V = FunctionSpace(mesh, 'P', 1)

#Boundary definition of function as 0 on boundary case
def boundary(x, onBoundary): 
    return onBoundary 

bc = DirichletBC(V, Constant(0), boundary)

#Gaussian diffusion equation with diffusion coefficient 3, initial value
u_O = Expression('exp(-a*pow(x[0], 2) - a*pow(x[1],2))', degree = 2, a = 3)

#Interpolation of function u_n as initial value function 
u_n = interpolate(u_O, V)
u_n.rename('u', 'initial value')

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)

#Calculation of general function of grad u*v 
F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
 
#Solution for a and l 
a, L = lhs(F), rhs(F)



#Exports to file for visualization
vtkfile = File('gaussian_diffusion.pvd')
vtkfile << (u_n, 0.0)

#Generating specific solution for all time intervals 
u = Function(V)
t = 0 
for n in range(steps): 
    t += dt
    solve (a == L, u, bc)
    
    vtkfile << (u, (float(t)))
    plot(u)
    time.sleep(0.1)
    
    u_n.assign(u)