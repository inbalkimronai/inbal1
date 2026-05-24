import math

E = 200e9
b = 0.002
A = b * b
I = b**4 / 12
Lside = 100.0
Ldiag = math.sqrt(2) * Lside

nodes = {
    'A': (0.0, 0.0),
    'B': (0.0, 100.0),
    'C': (100.0, 100.0),
    'D': (100.0, 0.0),
}

members = [
    ('AB', 'A', 'B', Lside),
    ('BC', 'B', 'C', Lside),
    ('CD', 'C', 'D', Lside),
    ('DA', 'D', 'A', Lside),
    ('AC', 'A', 'C', Ldiag),
    ('BD', 'B', 'D', Ldiag),
]


def element_stiffness(E, A, I, L, dx, dy):
    c = dx / L
    s = dy / L

    k_local = [
        [E * A / L, 0.0, 0.0, -E * A / L, 0.0, 0.0],
        [0.0, 12 * E * I / L**3, 6 * E * I / L**2, 0.0, -12 * E * I / L**3, 6 * E * I / L**2],
        [0.0, 6 * E * I / L**2, 4 * E * I / L, 0.0, -6 * E * I / L**2, 2 * E * I / L],
        [-E * A / L, 0.0, 0.0, E * A / L, 0.0, 0.0],
        [0.0, -12 * E * I / L**3, -6 * E * I / L**2, 0.0, 12 * E * I / L**3, -6 * E * I / L**2],
        [0.0, 6 * E * I / L**2, 2 * E * I / L, 0.0, -6 * E * I / L**2, 4 * E * I / L],
    ]

    T = [
        [c, s, 0.0, 0.0, 0.0, 0.0],
        [-s, c, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, c, s, 0.0],
        [0.0, 0.0, 0.0, -s, c, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    ]

    # k = T^T * k_local * T
    M = [[sum(k_local[r][k] * T[k][c] for k in range(6)) for c in range(6)] for r in range(6)]
    K = [[sum(T[r][k] * M[k][c] for k in range(6)) for c in range(6)] for r in range(6)]
    return K


node_dof = {node: i * 3 for i, node in enumerate(nodes)}
K = [[0.0 for _ in range(12)] for _ in range(12)]

for _, n1, n2, L in members:
    x1, y1 = nodes[n1]
    x2, y2 = nodes[n2]
    dx = x2 - x1
    dy = y2 - y1
    k = element_stiffness(E, A, I, L, dx, dy)
    dofs = [node_dof[n1] + 0, node_dof[n1] + 1, node_dof[n1] + 2,
            node_dof[n2] + 0, node_dof[n2] + 1, node_dof[n2] + 2]
    for a in range(6):
        for b in range(6):
            K[dofs[a]][dofs[b]] += k[a][b]

F = [0.0] * 12
F[node_dof['D'] + 1] = -1000.0

fixed = [node_dof['A'] + 0, node_dof['A'] + 1, node_dof['A'] + 2,
         node_dof['B'] + 0, node_dof['B'] + 1, node_dof['B'] + 2]
free = [i for i in range(12) if i not in fixed]

Kff = [[K[i][j] for j in free] for i in free]
Ff = [F[i] for i in free]

# Gaussian elimination
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

u = [0.0] * 12
for i, idx in enumerate(free):
    u[idx] = Ff[i]

R = [0.0] * 12
for i in range(12):
    R[i] = sum(K[i][j] * u[j] for j in range(12)) - F[i]

print("A reactions:", R[node_dof['A'] + 0], R[node_dof['A'] + 1], R[node_dof['A'] + 2])
print("B reactions:", R[node_dof['B'] + 0], R[node_dof['B'] + 1], R[node_dof['B'] + 2])
