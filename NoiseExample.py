from fenics import * 
#Subdivisions of x and t sections 
xSection = 100
tSection = 100

#Point definitions for boundary of domain
a = Point(0,0)
b = Point(100,100)

#Rectangular region of defined mesh 
mesh = RectangleMesh(a,b,xSection, tSection)
V = FunctionSpace(mesh, 'P', 1)

#Functional representation of B term in Van Deemter Equation
u_D = Expression('2*x[0]*x[0]/x[1]', degree = 2)

#Defines boundary function and conditions of boundary 
def boundary(x, onBoundary): 
    return onBoundary 
bc = DirichletBC(V, Constant(0), boundary)

#Vector field intepretations of the field V with constant f(0) as a solution
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)

#Interpolation of vector field by function 
u_n = interpolate(u_D, V)
u_n.rename('u', 'initial value')


#Use of vector fields and known solution allows function F as expansion 
F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx

#Left hand side and right hand sides of vector field 
a, L = lhs(F), rhs(F)

u = Function(V)
u.rename('u', 'solution')

#Function to find conditions such that a==L & u & bc 
#Specific solution for PDE system 
solve(a == L, u, bc)

plot (u)
plot(mesh)

#Exporting to vtk format for visualization
vtkfile = File('NoiseSolution.pvd')
vtkfile << u 


interactive()
