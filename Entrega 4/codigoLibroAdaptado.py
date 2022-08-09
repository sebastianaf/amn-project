""" From "COMPUTATIONAL PHYSICS" & "COMPUTER PROBLEMS in PHYSICS"
    by RH Landau, MJ Paez, and CC Bordeianu (deceased)
    Copyright R Landau, Oregon State Unv, MJ Paez, Univ Antioquia,
    C Bordeianu, Univ Bucharest, 2017.
    Please respect copyright & acknowledge our work."""

# Beam.py: solves Navier-Stokes equation for the flow around a beam

from numpy import *  # Needed for range
import pylab as p
from mpl_toolkits.mplot3d import Axes3D;
from matplotlib import image;

Nxmax = 25
Nymax = 25;  # Grid parameters
u = zeros((Nxmax + 1, Nymax + 1), float)  # Stream
w = zeros((Nxmax + 1, Nymax + 1), float)  # Vorticity
V0 = 1.0  # Initial v
omega = 0.1  # Relaxation param
IL = 7  # Geometry  Columna
H = 7   #Fila
T = 7
IL2 = 18
H2 = 7
T2 = 7
h = 1.
nu = 1.  # Viscosity
iter = 0  # Number iterations
R = V0 * h / nu  # Reynold number, normal units
print("Working, wait for the figure, count to 30")


def Mostrar(m):
    print("La matriz es la siguiente:")
    for fila in m:
        for valor in fila:
            print("\t", valor, end=" ")
        print()

def rellenar():
    # Viga 1
    for i in range(IL+1, IL + T):
        for j in range(0, H+1):
            u[i, j] = 0
            w[i, j] = 0

    for i in range(IL2, IL2 + T2):
        for j in range(Nymax - H2, Nymax):
            u[i, j] = 0
            w[i, j] = 0

def borders():  # Initialize stream,vorticity, sets BC
    for i in range(0, Nxmax + 1):  # Initialize stream function
        for j in range(0, Nymax + 1):  # Init vorticity
            w[i, j] = 0.
            u[i, j] = j * V0
    for i in range(0, Nxmax + 1):  # Fluid surface
        u[i, Nymax] = u[i, Nymax - 1] + V0 * h
        w[i, Nymax - 1] = 0.
    for j in range(0, Nymax + 1):
        u[1, j] = u[0, j]
        w[0, j] = 0.  # Inlet
    for i in range(0, Nxmax + 1):  # Centerline
        if i <= IL and i >= IL + T:
            u[i, 0] = 0.
            w[i, 0] = 0.
    for j in range(1, Nymax):  # Outlet
        w[Nxmax, j] = w[Nxmax - 1, j]
        u[Nxmax, j] = u[Nxmax - 1, j]  # Borders


def beam():  # BC for the beam
    for j in range(0, H + 1):  # Beam sides
        w[IL, j] = -2 * u[IL - 1, j] / (h * h)  # Front side
        w[IL + T, j] = -2 * u[IL + T + 1, j] / (h * h)  # Back side
    for i in range(IL, IL + T + 1):
        w[i, H - 1] = -2 * u[i, H] / (h * h);  # Top
    for i in range(IL, IL + T + 1):
        for j in range(0, H + 1):
            u[IL, j] = 0.  # Front
            u[IL + T, j] = 0.  # Back
            u[i, H] = 0;  # Top

def beam2():  # BC for the beam
    for j in range(Nymax - H2 - 1, Nymax):  # Beam sides
        w[IL2 - 1, j] = -2 * u[IL2 - 2, j] / (h * h)  # Front side
    for i in range(IL2, IL2 + T2 + 1):
        w[i, Nymax - H2 - 1] = -2 * u[i, Nymax - H2] / (h * h);  # Bottom
    for i in range(IL2, IL2 + T2 + 1):
        for j in range(Nymax - H2 - 1, Nymax):
            x=2
            u[IL2 - 1, j] = 0  # Front
            u[i, Nymax - H2 - 1] = 0  # Top

def relax():  # Method to relax stream
     # Reset conditions at beam
    beam()
    beam2()
    for i in range(1, Nxmax):  # Relax stream function
        for j in range(1, Nymax):
            r1 = omega * ((u[i + 1, j] + u[i - 1, j] + u[i, j + 1] + u[i, j - 1] + h * h * w[i, j]) / 4 - u[i, j])
            u[i, j] += r1
    for i in range(1, Nxmax):  # Relax vorticity
        for j in range(1, Nymax):
            a1 = w[i + 1, j] + w[i - 1, j] + w[i, j + 1] + w[i, j - 1]
            a2 = (u[i, j + 1] - u[i, j - 1]) * (w[i + 1, j] - w[i - 1, j])
            a3 = (u[i + 1, j] - u[i - 1, j]) * (w[i, j + 1] - w[i, j - 1])
            r2 = omega * ((a1 - (R / 4.) * (a2 - a3)) / 4.0 - w[i, j])
            w[i, j] += r2


m = 0
borders()
while (iter <= 2100):
    iter += 1
    if iter % 10 == 0:
        print(m)
        m += 1
    relax()
for i in range(0, Nxmax + 1):
    #v = 0.2
    for j in range(0, Nymax + 1):
        u[i, j] = u[i, j] / (V0 * h) #round(v-0.001,3)  # stream in V0h units
        #v=v-0.001
# u.resize((70,70));
# w.resize((70,70));
x = list(range(0, Nxmax - 1))  # to plot lines in x axis
y = list(range(0, Nymax - 1))
# x=range(0,69)                   #to plot lines in x axis
# y=range(0,69)
X, Y = p.meshgrid(x, y)  # grid for position and time

##############################################################
# Normalizamos las matrices halladas
def normalizar(u, w):
    for j in range(Nymax - 1):
        for i in range(Nxmax - 1):
            m = math.sqrt(pow(u[i][j], 2) + pow(w[i][j], 2))
            if (m != 0):
                u[i][j] = u[i][j] / m
                w[i][j] = w[i][j] / m

normalizar(u,w)
rellenar()

def functz(u):  # returns stream flow to plot
    z = u[X, Y]  # for several iterations
    return z


def functz1(w):  # returns stream flow to plot
    z1 = w[X, Y]  # for several iterations
    return z1


Z = functz(u)
Z1 = functz1(w)
Mostrar(Z)
Mostrar(Z1)
fig1 = p.figure()
p.title('Stream function - 2D Flow over a beam')
p.imshow(Z, origin='lower');
p.quiver(X,Y,Z,Z1)
p.colorbar();
fig2 = p.figure()
p.title('Vorticity - 2D Flow over a beam')
p.imshow(Z1, origin='lower');
p.quiver(X,Y,Z,Z1)
p.colorbar();
p.show()  # Shows the figure, close Python shell to
# Finish watching the figure
print("finished")