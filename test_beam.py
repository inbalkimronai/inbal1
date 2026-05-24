import math

E = 200e9
b = 0.002
A = b * b
I = b**4 / 12
L = 100.0

k = [
    [E * A / L, 0.0, 0.0, -E * A / L, 0.0, 0.0],
    [0.0, 12 * E * I / L**3, 6 * E * I / L**2, 0.0, -12 * E * I / L**3, 6 * E * I / L**2],
    [0.0, 6 * E * I / L**2, 4 * E * I / L, 0.0, -6 * E * I / L**2, 2 * E * I / L],
    [-E * A / L, 0.0, 0.0, E * A / L, 0.0, 0.0],
    [0.0, -12 * E * I / L**3, -6 * E * I / L**2, 0.0, 12 * E * I / L**3, -6 * E * I / L**2],
    [0.0, 6 * E * I / L**2, 2 * E * I / L, 0.0, -6 * E * I / L**2, 4 * E * I / L],
]
F = [0.0, -1000.0, 0.0, 0.0, 0.0, 0.0]
fixed = [0, 1, 2]
free = [3, 4, 5]
Kff = [[k[i][j] for j in free] for i in free]
Ff = [F[i] for i in free]

n = len(free)
for i in range(n):
    pivot = max(range(i, n), key=lambda r: abs(Kff[r][i]))
    Kff[i], Kff[pivot] = Kff[pivot], Kff[i]
    Ff[i], Ff[pivot] = Ff[pivot], Ff[i]
    piv = Kff[i][i]
    for j in range(i, n):
        Kff[i][j] /= piv
    Ff[i] /= piv
    for r in range(n):
        if r != i:
            factor = Kff[r][i]
            for c in range(i, n):
                Kff[r][c] -= factor * Kff[i][c]
            Ff[r] -= factor * Ff[i]

u = [0.0] * 6
for i, idx in enumerate(free):
    u[idx] = Ff[i]

R = [sum(k[i][j] * u[j] for j in range(6)) - F[i] for i in range(6)]
print(R)
