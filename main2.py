import matplotlib.pyplot as plt
import numpy as np
import math
from pprint import pprint
from numpy import array, zeros, diag, diagflat, dot

# Valores temporalmente por defecto
# Distancia entre puntos
h = 1
omega = 0.00001
# Dimensiones de la matriz
Nxmax = 25
Nymax = 25
# Viga 1
inicio = 7
alto1 = 7
ancho1 = 7
# Viga 2
inicio2 = 19
alto2 = 5
ancho2 = 7
V0 = 25  # Velocidad inicial

u = np.zeros((Nxmax, Nymax), float)
w = np.zeros((Nxmax, Nymax), float)

bu = np.zeros(Nxmax * Nymax)
bw = np.zeros(Nxmax * Nymax)

for i in range(Nxmax * Nymax):
  bu[i] = h/8
  bw[i] = h/8

xu0 = np.zeros(Nxmax*Nymax)
xw0 = np.zeros(Nxmax*Nymax)

def Mostrar(m):
    print("La matriz es la siguiente:")
    for fila in m:
        for valor in fila:
            print("\t", valor, end=" ")
        print()

def inicializar():
    for i in range(Nxmax):
        V = 0.3
        for j in range(Nymax):
            u[i, j] = round(V - 0.001, 3)
            V = V - 0.001

inicializar()

def gen_matriz_sis_lineal(n,u,w,h,tipo):
    a = 0
    b = 0
    c = 0
    d = 0
    matriz = []
    for i in range(1, (n ** 2) + 1):
        fila = []
        for j in range(1, (n ** 2) + 1):
            if (tipo==1):
                a = (h / 8) * u[math.ceil(i/n) - 1, math.ceil(j/n) - 1] - 1 / 4
                b = -(h / 8) * u[math.ceil(i/n) - 1, math.ceil(j/n) - 1] - 1 / 4
                c = (h / 8) * w[math.ceil(i/n) - 1, math.ceil(j/n) - 1] - 1 / 4
                d = -(h / 8) * w[math.ceil(i/n) - 1, math.ceil(j/n) - 1] - 1 / 4
            else:
                a = (h / 8) * w[math.ceil(i/n) - 1, math.ceil(j/n) - 1] - 1 / 4
                b = -(h / 8) * w[math.ceil(i/n) - 1, math.ceil(j/n) - 1] - 1 / 4
                c = (h / 8) * u[math.ceil(i/n) - 1, math.ceil(j/n) - 1] - 1 / 4
                d = -(h / 8) * u[math.ceil(i/n) - 1, math.ceil(j/n) - 1] - 1 / 4
            # Primera diagonal
            if (i == j):
                fila.append(1)

            # Segunda diagonal superior
            elif (
                    (i % n != 0 and i + 1 == j)):
                fila.append(a)

            # Segunda diagonal inferior
            elif (
                    ((i % n) - 1 != 0 and i == j + 1)):
                fila.append(b)

            # Tercera diagonal superior
            elif (i + n == j):
                fila.append(c)

            # Tercera diagonal inferior
            elif (i == j + n):
                fila.append(d)

            else:
                fila.append(0)

        matriz.append(fila)
    return matriz


def surfaceG(m, b, flag):
    n = len(u)

    for i in range(0, n):
        if flag == "u":
            b[i] = 2 * h * V0
        else:
            b[i] = 0
        for j in range(0, n * n):
            # print(f'{i},{j}')
            m[i][j] = 0
            if i == j:
                m[i][j] = 1
    if flag == "u":
        for i in range(0, n):
            for j in range(0, n * n):
                if (i == j + 1 and j < n - 1):
                    m[i][j] = -1
                if (j == i + 1 and j < n):
                    m[i][j] = 1
    return m, b


def InletF(m, b, flag):
    n = len(u)
    for i in range(0, n * n, n):
        if flag == "u":
            b[i] = 2 * h * V0
        else:
            b[i] = 0
        for j in range(0, n * n):
            m[i][j] = 0
            if i == j:
                m[i][j] = 1
    if flag == "u":
        for i in range(0, n * n, n):
            for j in range(0, n * n):
                if (i + 1 == j):
                    m[i][j] = -1

    else:
        for i in range(0, n * n, n):
            for j in range(0, n * n):
                if (i + n == j):
                    m[i][j] = -1
                if (i == j + n):
                    m[i][j] = 1
    return m, b


def outlet(m, b, flag):
    n = len(u)
    for i in range(n - 1, n * n, n):
        b[i] = 0
        for j in range(0, n * n):
            m[i][j] = 0
            if i == j:
                m[i][j] = 1

    for i in range(n - 1, n * n, n):
        for j in range(0, n * n):
            if (i + n == j):
                m[i][j] = 1
            if (i == j + n):
                m[i][j] = -1

    return m, b

def centerLine(m,b,flag):
  n=len(u)
  for i in range(n*n-n,n*n):
    b[i] = 0
    for j in range(0,n*n):
      m[i][j] = 0
      if i == j:
        m[i][j] = 1
  return m,b

def richardson(A,x,b,N):

  for i in range(N):
    r=b-np.dot(A,x)
    x = x + r

  return x

uJac = gen_matriz_sis_lineal(Nxmax, u, w, 1, 1)
wJac = gen_matriz_sis_lineal(Nxmax, u, w, 1, 2)
surfaceG(uJac, bu, "u")
surfaceG(wJac, bw, "w")
InletF(uJac, bu, "u")
InletF(wJac, bw, "w")
outlet(uJac, bu, "u")
outlet(wJac, bw, "w")
centerLine(uJac, bu, "u")
centerLine(wJac, bw, "w")
# Mostrar(uJac)
# Mostrar(wJac)
for i in range(200):
    xu0 = richardson(uJac, xu0, bu, 1)
    xw0 = richardson(wJac, xw0, bw, 1)
    it = 0
    for i in range(0, Nxmax):
        for j in range(0, Nymax):
            u[i, j] = u[i, j] + omega * xu0[it]
            w[i, j] = w[i, j] + omega * xw0[it]
            it += 1
    uJac = gen_matriz_sis_lineal(Nxmax, u, w, 1, 1)
    wJac = gen_matriz_sis_lineal(Nxmax, u, w, 1, 2)

# Mostrar(u)
# Mostrar(w)

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

#############################################################
# Matriz de magnitudes

magn = np.zeros((Nxmax, Nxmax))

for j in range(Nxmax):
    for i in range(Nxmax):
        a = math.sqrt(pow(u[i][j], 2) + pow(w[i][j], 2))
        magn[i][j] = a

##############################################################
# Vectores

x = np.linspace(0, Nxmax - 1, Nxmax)
y = np.linspace(0, Nymax - 1, Nymax)

xmesh, ymesh = np.meshgrid(x, y)

umesh = u
vmesh = w

#####################################################
# Graficar
plt.imshow(magn)
plt.quiver(xmesh, ymesh, umesh, vmesh)
plt.colorbar()
plt.show()

