const E = 200e9;
const b = 0.002;
const A = b * b;
const I = Math.pow(b, 4) / 12;
const Lside = 100.0;
const Ldiag = Math.sqrt(2) * Lside;

const nodes = {
  A: [0.0, 0.0],
  B: [0.0, 100.0],
  C: [100.0, 100.0],
  D: [100.0, 0.0],
};

const members = [
  ['AB', 'A', 'B', Lside],
  ['BC', 'B', 'C', Lside],
  ['CD', 'C', 'D', Lside],
  ['DA', 'D', 'A', Lside],
  ['AC', 'A', 'C', Ldiag],
  ['BD', 'B', 'D', Ldiag],
];

function multiplyMatrices(a, b) {
  const rows = a.length;
  const cols = b[0].length;
  const shared = a[0].length;
  const result = Array.from({ length: rows }, () => Array(cols).fill(0.0));

  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < cols; j++) {
      let sum = 0.0;
      for (let k = 0; k < shared; k++) {
        sum += a[i][k] * b[k][j];
      }
      result[i][j] = sum;
    }
  }

  return result;
}

function elementStiffness(E, A, I, L, dx, dy) {
  const c = dx / L;
  const s = dy / L;

  const kLocal = [
    [E * A / L, 0.0, 0.0, -E * A / L, 0.0, 0.0],
    [0.0, 12 * E * I / Math.pow(L, 3), 6 * E * I / Math.pow(L, 2), 0.0, -12 * E * I / Math.pow(L, 3), 6 * E * I / Math.pow(L, 2)],
    [0.0, 6 * E * I / Math.pow(L, 2), 4 * E * I / L, 0.0, -6 * E * I / Math.pow(L, 2), 2 * E * I / L],
    [-E * A / L, 0.0, 0.0, E * A / L, 0.0, 0.0],
    [0.0, -12 * E * I / Math.pow(L, 3), -6 * E * I / Math.pow(L, 2), 0.0, 12 * E * I / Math.pow(L, 3), -6 * E * I / Math.pow(L, 2)],
    [0.0, 6 * E * I / Math.pow(L, 2), 2 * E * I / L, 0.0, -6 * E * I / Math.pow(L, 2), 4 * E * I / L],
  ];

  const T = [
    [c, s, 0.0, 0.0, 0.0, 0.0],
    [-s, c, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, c, s, 0.0],
    [0.0, 0.0, 0.0, -s, c, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
  ];

  return multiplyMatrices(multiplyMatrices(transpose(T), kLocal), T);
}

function transpose(matrix) {
  return matrix[0].map((_, col) => matrix.map(row => row[col]));
}

const nodeDof = Object.fromEntries(Object.keys(nodes).map((node, i) => [node, i * 3]));
const K = Array.from({ length: 12 }, () => Array(12).fill(0.0));

for (const [, n1, n2, L] of members) {
  const [x1, y1] = nodes[n1];
  const [x2, y2] = nodes[n2];
  const dx = x2 - x1;
  const dy = y2 - y1;
  const k = elementStiffness(E, A, I, L, dx, dy);
  const dofs = [nodeDof[n1], nodeDof[n1] + 1, nodeDof[n1] + 2, nodeDof[n2], nodeDof[n2] + 1, nodeDof[n2] + 2];

  for (let a = 0; a < 6; a++) {
    for (let b = 0; b < 6; b++) {
      K[dofs[a]][dofs[b]] += k[a][b];
    }
  }
}

const F = Array(12).fill(0.0);
F[nodeDof.D + 1] = -1000.0;

const fixed = [
  nodeDof.A, nodeDof.A + 1, nodeDof.A + 2,
  nodeDof.B, nodeDof.B + 1, nodeDof.B + 2,
];

const free = Array.from({ length: 12 }, (_, i) => i).filter(i => !fixed.includes(i));

const Kff = free.map(i => free.map(j => K[i][j]));
const Ff = free.map(i => F[i]);

for (let i = 0; i < free.length; i++) {
  let pivot = i;
  for (let r = i + 1; r < free.length; r++) {
    if (Math.abs(Kff[r][i]) > Math.abs(Kff[pivot][i])) {
      pivot = r;
    }
  }

  [Kff[i], Kff[pivot]] = [Kff[pivot], Kff[i]];
  [Ff[i], Ff[pivot]] = [Ff[pivot], Ff[i]];

  const piv = Kff[i][i];
  for (let j = i; j < free.length; j++) {
    Kff[i][j] /= piv;
  }
  Ff[i] /= piv;

  for (let r = 0; r < free.length; r++) {
    if (r === i) continue;
    const factor = Kff[r][i];
    for (let c = i; c < free.length; c++) {
      Kff[r][c] -= factor * Kff[i][c];
    }
    Ff[r] -= factor * Ff[i];
  }
}

const u = Array(12).fill(0.0);
for (let i = 0; i < free.length; i++) {
  u[free[i]] = Ff[i];
}

const R = Array.from({ length: 12 }, (_, i) => K[i].reduce((sum, val, j) => sum + val * u[j], 0) - F[i]);

console.log('A reactions:', R[nodeDof.A], R[nodeDof.A + 1], R[nodeDof.A + 2]);
console.log('B reactions:', R[nodeDof.B], R[nodeDof.B + 1], R[nodeDof.B + 2]);
