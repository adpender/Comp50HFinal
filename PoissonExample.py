from fenics import *

# Create mesh and define function space
mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'P', 1)

# Boundary Conditions of Function
u_D = Expression('1 + x[0]*x[0]*x[0] + 2*x[1]*x[1]', degree=3)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Plot solution and mesh
plot(u)
plot(mesh)

# VTK format of Solution 
vtkfile = File('poisson/solution2.pvd')
vtkfile << u

# Hold plot
interactive()