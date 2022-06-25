import matplotlib.pyplot as plt
import numpy as np
import math

# Valores temporalmente por defecto
# Distancia entre puntos
h = 1
# Parametro de relajacion
omega = 0.1
Nxmax = 20
Nymax = 20
# Viga 1
inicio = 5
alto1 = 6
ancho1 = 6
# Viga 2
inicio2 = 15
alto2 = 4
ancho2 = 6
V0 = 25  # Velocidad inicial

u = np.zeros((Nxmax, Nymax), float)
w = np.zeros((Nxmax, Nymax), float)

u2 = np.zeros((Nxmax, Nymax), float)
w2 = np.zeros((Nxmax, Nymax), float)

def inicializar():
    for i in range(0, Nxmax):  # Inicializamos los valores de la matriz
        for j in range(0, Nymax):
            w[i, j] = 0
            u[i, j] = V0 - j

inicializar()

def viga1():
    for i in range(Nymax - alto1 + 1, Nymax):
        w[i, inicio] = - 2 * (u[i, inicio - 1] - u[i, inicio]) / (h * h)  # Front
        w[i, inicio + ancho1 - 1] = - 2 * (u[i, inicio + ancho1] - u[i, inicio + ancho1 - 1]) / (h * h)  # Back
    for j in range(inicio, inicio + ancho1):
        w[Nymax - alto1, j] = - 2 * (u[Nymax - alto1 - 1, j] - u[Nymax - alto1, j]) / (h * h);  # top
    for j in range(inicio, inicio + ancho1):
        for i in range(Nymax - alto1 + 1, Nymax):
            u[i, inicio] = 0  # Front
            u[i, inicio + ancho1 - 1] = 0  # Back
            u[Nymax - alto1, j] = 0  # Top


def viga2():
    for i in range(0, alto2):
        w[i, inicio2] = 2 * (u[i, inicio2 - 1] - u[i, inicio2]) / (h * h)  # Front
    for j in range(inicio2, inicio2 + ancho2 - 1):
        w[alto2 - 1, j] = 2 * (u[alto2, j] - u[alto2 - 1, j]) / (h * h)  # Bottom
    for j in range(inicio2, inicio2 + ancho2 - 1):
        for i in range(0, alto2):
            u[i, inicio2] = 0.  # Back
            u[alto2 - 1, j] = 0.  # Bottom

def rellenar():
    for i in range(0, alto2 - 1):
        for j in range(inicio2 + 1, inicio2 + ancho2 - 1):
            u[i, j] = 0
            w[i, j] = 0
    for i in range(Nxmax - alto1 + 1, Nxmax):
        for j in range(inicio + 1, inicio + ancho1 - 1):
            u[i, j] = 0
            w[i, j] = 0

def Mostrar(m):
    print("La matriz es la siguiente:")
    for fila in m:
        for valor in fila:
            print("\t", valor, end=" ")
        print()

def run():
    viga1()
    viga2()
    for i in range(1, Nxmax - 1):  # Relax stream
        for j in range(1, Nymax - 1):
            r1 = omega * ((u[i + 1, j] + u[i - 1, j] + u[i, j + 1] + u[i, j - 1]
                           + h * h * w[i, j]) / 4 - u[i, j])
            if i >= alto2 and (j < inicio or j > inicio + ancho1 - 1):
                u[i, j] += r1
            if alto2 <= i < (Nymax - alto1) and (inicio <= j <= inicio + ancho1 - 1):
                u[i, j] += r1
            if i < alto2 and j <= Nxmax - ancho2:
                u[i, j] = r1
    for i in range(1, Nxmax - 1):  # Relax vorticity
        for j in range(1, Nymax - 1):
            r2 = omega * ((w[i + 1, j] + w[i - 1, j] + w[i, j + 1] + w[i, j - 1]
                           + h * h * u[i, j]) / 4 - w[i, j])
            if i >= alto2 and (j < inicio or j > inicio + ancho1 - 1):
                w[i, j] += r2
            if alto2 <= i < (Nymax - alto1) and (inicio <= j <= inicio + ancho1 - 1):
                w[i, j] += r2
            if i < alto2 and j <= Nxmax - ancho2:
                w[i, j] += r2
    rellenar()

for i in range(10):
    run()

Mostrar(u)
Mostrar(w)

#############################################################
# Matriz de magnitudes

magn = np.zeros((Nxmax, Nxmax))

for j in range(Nxmax):
    for i in range(Nxmax):
        a = math.sqrt(pow(u[i][j], 2) + pow(w[i][j], 2))
        magn[i][j] = a

# Mostrar(magn)

##############################################################
# Normalizamos las matrices halladas
def normalizar():
    for j in range(Nymax):
        for i in range(Nxmax):
            m = math.sqrt(pow(u[i][j], 2) + pow(w[i][j], 2))
            if (m != 0):
                u[i][j] = u[i][j] / m
                w[i][j] = w[i][j] / m

# normalizar()

##############################################################
# Vectores

x = np.linspace(0, 19, 20)
y = np.linspace(0, 19, 20)

xmesh, ymesh = np.meshgrid(x, y)

umesh = u
vmesh = w

#####################################################
# Graficar
plt.imshow(magn)
plt.quiver(xmesh, ymesh, umesh, vmesh)
plt.colorbar()
plt.show()