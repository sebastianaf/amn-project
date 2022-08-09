import matplotlib.pyplot as plt
import numpy as np
import math
from pprint import pprint
from numpy import array, zeros, diag, diagflat, dot

# Valores temporalmente por defecto
# Distancia entre puntos
h = 1
omegaX = 0.1
omegaY = 0.1
# Dimensiones de la matriz
Nxmax = 25
Nymax = 25
# Viga 1
inicio = 7
alto1 = 7
ancho1 = 7
# Viga 2
inicio2 = 18
alto2 = 7
ancho2 = 7
V0 = 0.2  # Velocidad inicial

u = np.zeros((Nxmax, Nymax), float)
w = np.zeros((Nxmax, Nymax), float)

bu = np.zeros(Nxmax * Nymax)
bw = np.zeros(Nxmax * Nymax)

for i in range(Nxmax * Nymax):
  bu[i] = h/8
  bw[i] = 0

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
        V = 0.2
        for j in range(Nymax):
            u[i, j] = round(V - 0.001, 3)
            V = V - 0.001
    '''
    for i in range(Nxmax):
        for j in range(Nymax):
            u[i, j] = 0.01
    '''

inicializar()

def rellenar(m1, m2):
    for i in range(Nxmax - alto1, Nxmax):
        for j in range(inicio, inicio + ancho1):
            m1[i, j] = 0
            m2[i, j] = 0

    for i in range(0, alto2):
        for j in range(inicio2, inicio2 + ancho2):
            m1[i, j] = 0
            m2[i, j] = 0

rellenar(u, w)

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
            # b[i] = V0
            b[i] = u[1][i]
        else:
            b[i] = 0
        for j in range(0, n * n):
            m[i][j] = 0
            if i == j:
                m[i][j] = 1
    return m, b


def InletF(m, b, flag):
    n = len(u)
    for i in range(0, n * n, n):
        if flag == "u":
            b[i] = V0
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
        i2 = math.ceil(i / Nymax) - 1
        if flag == "u":
            b[i] = u[i2][Nxmax - 2]
        else:
            b[i] = 0
        for j in range(0, n * n):
            m[i][j] = 0
            if i == j:
                m[i][j] = 1
    '''
    for i in range(n - 1, n * n, n):
        for j in range(0, n * n):
            if (i + n == j):
                m[i][j] = 1
            if (i == j + n):
                m[i][j] = -1
    '''

    return m, b

def centerLine(m, b, flag):
    n = len(u)
    for i in range(n * n - n, n * n):
        i2 = math.ceil(i / Nymax) - 1
        if flag == "u":
            b[i] = u[Nymax - 2][i2]
        else:
            b[i] = 0
        for j in range(0, n * n):
            m[i][j] = 0
            if i == j:
                m[i][j] = 1
    return m, b

def transformPairToOne(i, j):
    return (i+1) * Nxmax - (Nxmax - j)

def transform(j):
    return j%Nxmax

def llenar(m, b, flag):
    n = len(u)
    # Rellenar
    # Viga 1
    for i in range(Nxmax - alto1, Nxmax):
        for k in range(inicio, inicio + ancho1):
            for j in range(0, n * n):
                j2 = math.ceil(j / Nymax) - 1
                m[transformPairToOne(i, k)][j] = 0
                if transformPairToOne(i, k) == j:
                    m[transformPairToOne(i, k)][j] = 1

    for i in range(Nxmax - alto1, Nxmax):
        for k in range(inicio, inicio + ancho1):
            b[transformPairToOne(i, k)] = 0

    # Viga 2
    for i in range(0, alto2):
        for k in range(inicio2, inicio2 + ancho2):
            for j in range(0, n * n):
                j2 = math.ceil(j / Nymax) - 1
                m[transformPairToOne(i, k)][j] = 0
                if transformPairToOne(i, k) == j:
                    m[transformPairToOne(i, k)][j] = 1

    for i in range(0, alto2):
        for k in range(inicio2, inicio2 + ancho2):
            b[transformPairToOne(i, k)] = 0

def viga1(m, b, flag):
    n = len(u)

    # Pared izquierda viga 1
    for i in range(n - alto1, n):
        for j in range(0, n * n):
            if flag == "u":
                b[i * n + inicio - 1] = 0
            else:
                j2 = math.ceil(j / Nymax) - 1
                b[i * n + inicio - 1] = -2 * (u[i][j2 - 1] - u[i][j2]) / h * h
    # Pared superior viga 1
    f = n - alto1 - 1
    for j in range(inicio, inicio + ancho1 + 1):
        if flag == "u":
            b[transformPairToOne(f, j)] = 0
        else:
            j2 = math.ceil(j / Nymax) - 1
            i = n - alto1 - 1
            b[transformPairToOne(f, j)] = -2 * (u[i - 1][j] - u[i][j]) / h * h
    # Pared derecha viga 1
    for i in range(n - alto1 - 1, n):
        for j in range(0, n * n):
            if flag == "u":
                b[i * n + inicio + ancho1] = 0
            else:
                j2 = math.ceil(j / Nymax) - 1
                if j2 == inicio + ancho1:
                    b[i * n + inicio + ancho1] = -2 * (u[i][j2 + 1] - u[i][j2]) / h * h

    if flag == "u":
        # Pared izquierda viga 1
        for i in range(n - alto1, n):
            for j in range(0, n * n):
                m[i * n + inicio - 1][j] = 0
                if i * n + inicio - 1 == j:
                    m[i * n + inicio - 1][j] = 1
        # Pared superior viga 1
        f = n - alto1 - 1
        for j in range(inicio, inicio + ancho1):
            for k in range(0, n * n):
                m[transformPairToOne(f, j)][k] = 0
                if transformPairToOne(f, j) == k:
                    m[transformPairToOne(f, j)][k] = 1
        # Pared derecha viga 1
        for i in range(n - alto1 - 1, n):
            for j in range(0, n * n):
                m[i * n + inicio + ancho1][j] = 0
                if i * n + inicio + ancho1 == j:
                    m[i * n + inicio + ancho1][j] = 1
    else:
        # Pared izquierda viga 1
        for i in range(n - alto1, n):
            for j in range(0, n * n):
                m[i * n + inicio - 1][j] = 0
                if i * n + inicio - 1 == j:
                    m[i * n + inicio - 1][j] = 1

        # Pared superior viga 1
        f = n - alto1 - 1
        for j in range(inicio, inicio + ancho1):
            for k in range(0, n * n):
                m[transformPairToOne(f, j)][k] = 0
                if transformPairToOne(f, j) == k:
                    m[transformPairToOne(f, j)][k] = 1
        # Pared derecha viga 1
        for i in range(n - alto1 - 1, n):
            for j in range(0, n * n):
                m[i * n + inicio + ancho1][j] = 0
                if i * n + inicio + ancho1 == j:
                    m[i * n + inicio + ancho1][j] = 1

def viga2(m, b, flag):
    n = len(u)

    # Pared izquierda viga 2
    for i in range(0, alto2 + 1):
        for j in range(0, n * n):
            if flag == "u":
                b[i * n + inicio2 - 1] = 0
            else:
                j2 = Nymax - ancho2 - 1
                # print(i, j2)
                b[i * n + inicio2 - 1] = -2 * (u[i][j2 - 1] - u[i][j2]) / h * h
                break

    # Pared inferior viga 2
    f = alto2
    for j in range(inicio2, inicio2 + ancho2):
        if flag == "u":
            b[transformPairToOne(f, j)] = 0
        else:
            j2 = math.ceil(j / Nymax) - 1
            i = alto2
            b[transformPairToOne(f, j)] = -2 * (u[i + 1][j] - u[i][j]) / h * h

    if flag == "u":
        # Pared izquierda viga 2
        for i in range(0, alto2 + 1):
            for j in range(0, n * n):
                m[i * n + inicio2 - 1][j] = 0
                if i * n + inicio2 - 1 == j:
                    m[i * n + inicio2 - 1][j] = 1
        # Pared inferior viga 2
        f = alto2
        for j in range(inicio2, inicio2 + ancho2):
            for k in range(0, n * n):
                m[transformPairToOne(f, j)][k] = 0
                if transformPairToOne(f, j) == k:
                    m[transformPairToOne(f, j)][k] = 1
    else:
        # Pared izquierda viga 2
        for i in range(0, alto2 + 1):
            for j in range(0, n * n):
                m[i * n + inicio2 - 1][j] = 0
                if i * n + inicio2 - 1 == j:
                    m[i * n + inicio2 - 1][j] = 1

        # Pared inferior viga 2
        f = alto2
        for j in range(inicio2, inicio2 + ancho2):
            for k in range(0, n * n):
                m[transformPairToOne(f, j)][k] = 0
                if transformPairToOne(f, j) == k:
                    m[transformPairToOne(f, j)][k] = 1

def Jacobi(A,x,b,N):
  D = np.diag(A).reshape(x.shape)
  R = A - np.diagflat(D)
  for i in range(N):
    x = (b - np.dot(R,x))/D
  return x

uJac = gen_matriz_sis_lineal(Nxmax, u, w, 1, 1)
wJac = gen_matriz_sis_lineal(Nxmax, u, w, 1, 2)
def condiciones(mu, mw, b1, b2):
    surfaceG(mu, b1, "u")
    surfaceG(mw, b2, "w")
    InletF(mu, b1, "u")
    InletF(mw, b2, "w")
    outlet(mu, b1, "u")
    outlet(mw, b2, "w")
    centerLine(mu, b1, "u")
    centerLine(mw, b2, "w")
    llenar(mu, b1, "u")
    llenar(mw, b2, "w")
    viga1(mu, b1, "u")
    viga1(mw, b2, "w")
    viga2(mu, b1, "u")
    viga2(mw, b2, "w")

condiciones(uJac, wJac, bu, bw)
# Mostrar(uJac)
# Mostrar(wJac)
for i in range(1):
    xu0 = Jacobi(uJac, xu0, bu, 10)
    xw0 = Jacobi(wJac, xw0, bw, 10)
    it = 0
    for i in range(0, Nxmax):
        for j in range(0, Nymax):
            u[i, j] = u[i, j] + omegaX * xu0[it]
            w[i, j] = w[i, j] + omegaY * xw0[it]
            it += 1
    newUJac = gen_matriz_sis_lineal(Nxmax, u, w, 1, 1)
    newWJac = gen_matriz_sis_lineal(Nxmax, u, w, 1, 2)
    condiciones(newUJac, newWJac, bu, bw)


Mostrar(u)
Mostrar(w)

#############################################################
# Matriz de magnitudes

magn = np.zeros((Nxmax, Nxmax))

for j in range(Nxmax):
    for i in range(Nxmax):
        a = math.sqrt(pow(u[i][j], 2) + pow(w[i][j], 2))
        magn[i][j] = a
        # magn[i][j] = abs(u[i][j]) + abs(w[i][j])
        # magn[i][j] = abs(u[i][j]+w[i][j])/2
        # magn[i][j] = [i][j]
Mostrar(magn)
##############################################################
# Normalizamos las matrices halladas
def normalizar():
    for j in range(Nymax):
        for i in range(Nxmax):
            m = math.sqrt(pow(u[i][j], 2) + pow(w[i][j], 2))
            if (m != 0):
                u[i][j] = u[i][j] / m
                w[i][j] = w[i][j] / m

normalizar()

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
plt.colorbar()
plt.quiver(xmesh, ymesh, umesh, vmesh)
plt.show()