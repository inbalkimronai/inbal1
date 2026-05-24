const http = require('http');

function fetchProblem() {
  return new Promise((resolve, reject) => {
    http.get('http://192.168.60.104:8451/api/problem', (res) => {
      let data = '';
      res.on('data', (chunk) => {
        data += chunk;
      });
      res.on('end', () => {
        try {
          resolve(JSON.parse(data));
        } catch (err) {
          reject(err);
        }
      });
    }).on('error', reject);
  });
}

function transpose(matrix) {
  return matrix[0].map((_, col) => matrix.map((row) => row[col]));
}

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

function multiplyMatrixVector(matrix, vector) {
  return matrix.map((row) => row.reduce((sum, value, index) => sum + value * vector[index], 0.0));
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

function solveLinearSystem(K, F) {
  const n = K.length;
  const A = K.map((row, i) => row.slice().concat([F[i]]));

  for (let i = 0; i < n; i++) {
    let pivot = i;
    for (let r = i + 1; r < n; r++) {
      if (Math.abs(A[r][i]) > Math.abs(A[pivot][i])) {
        pivot = r;
      }
    }
    if (Math.abs(A[pivot][i]) < 1e-16) {
      throw new Error('Singular stiffness matrix');
    }

    [A[i], A[pivot]] = [A[pivot], A[i]];

    const piv = A[i][i];
    for (let j = i; j <= n; j++) {
      A[i][j] /= piv;
    }

    for (let r = 0; r < n; r++) {
      if (r === i) {
        continue;
      }
      const factor = A[r][i];
      for (let c = i; c <= n; c++) {
        A[r][c] -= factor * A[i][c];
      }
    }
  }

  return A.map((row) => row[n]);
}

async function main() {
  const problem = await fetchProblem();

  const E = 200e9;
  const rodSide = problem.rodSide;
  const A = rodSide * rodSide;
  const I = Math.pow(rodSide, 4) / 12;
  const yieldMPa = problem.yieldMPa * 1e6;

  const nodes = problem.nodes;
  const members = problem.members;
  const boundary = problem.boundary;

  const nodeCount = nodes.length;
  const dofCount = nodeCount * 3;

  const fixed = new Set();
  for (let i = 0; i < nodeCount; i++) {
    const [x] = nodes[i];
    if (x === boundary.x_min || x === boundary.x_max) {
      fixed.add(i * 3);
      fixed.add(i * 3 + 1);
      fixed.add(i * 3 + 2);
    }
  }

  const free = [];
  for (let i = 0; i < dofCount; i++) {
    if (!fixed.has(i)) {
      free.push(i);
    }
  }

  const K = Array.from({ length: dofCount }, () => Array(dofCount).fill(0.0));

  for (const [i, j] of members) {
    const [x1, y1] = nodes[i];
    const [x2, y2] = nodes[j];
    const dx = x2 - x1;
    const dy = y2 - y1;
    const L = Math.hypot(dx, dy);
    const k = elementStiffness(E, A, I, L, dx, dy);

    const dofs = [i * 3, i * 3 + 1, i * 3 + 2, j * 3, j * 3 + 1, j * 3 + 2];
    for (let a = 0; a < 6; a++) {
      for (let b = 0; b < 6; b++) {
        K[dofs[a]][dofs[b]] += k[a][b];
      }
    }
  }

  const F = Array(dofCount).fill(0.0);
  for (const load of problem.loads) {
    const [nodeIndex, fx, fy] = load;
    F[nodeIndex * 3] += fx;
    F[nodeIndex * 3 + 1] += fy;
  }

  const Kff = free.map((rowDof) => free.map((colDof) => K[rowDof][colDof]));
  const Ff = free.map((rowDof) => F[rowDof]);
  const u = Array(dofCount).fill(0.0);
  const freeDisp = solveLinearSystem(Kff, Ff);
  for (let i = 0; i < free.length; i++) {
    u[free[i]] = freeDisp[i];
  }

  const reactions = Array(dofCount).fill(0.0);
  for (let i = 0; i < dofCount; i++) {
    reactions[i] = K[i].reduce((sum, value, j) => sum + value * u[j], 0.0) - F[i];
  }

  const T = [
    [1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1],
  ];

  let maxStress = 0.0;
  const memberForces = [];

  for (const [i, j] of members) {
    const [x1, y1] = nodes[i];
    const [x2, y2] = nodes[j];
    const dx = x2 - x1;
    const dy = y2 - y1;
    const L = Math.hypot(dx, dy);
    const c = dx / L;
    const s = dy / L;
    const Tlocal = [
      [c, s, 0.0, 0.0, 0.0, 0.0],
      [-s, c, 0.0, 0.0, 0.0, 0.0],
      [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
      [0.0, 0.0, 0.0, c, s, 0.0],
      [0.0, 0.0, 0.0, -s, c, 0.0],
      [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    ];
    const kLocal = [
      [E * A / L, 0.0, 0.0, -E * A / L, 0.0, 0.0],
      [0.0, 12 * E * I / Math.pow(L, 3), 6 * E * I / Math.pow(L, 2), 0.0, -12 * E * I / Math.pow(L, 3), 6 * E * I / Math.pow(L, 2)],
      [0.0, 6 * E * I / Math.pow(L, 2), 4 * E * I / L, 0.0, -6 * E * I / Math.pow(L, 2), 2 * E * I / L],
      [-E * A / L, 0.0, 0.0, E * A / L, 0.0, 0.0],
      [0.0, -12 * E * I / Math.pow(L, 3), -6 * E * I / Math.pow(L, 2), 0.0, 12 * E * I / Math.pow(L, 3), -6 * E * I / Math.pow(L, 2)],
      [0.0, 6 * E * I / Math.pow(L, 2), 2 * E * I / L, 0.0, -6 * E * I / Math.pow(L, 2), 4 * E * I / L],
    ];

    const elementDisp = [u[i * 3], u[i * 3 + 1], u[i * 3 + 2], u[j * 3], u[j * 3 + 1], u[j * 3 + 2]];
    const localDisp = multiplyMatrixVector(Tlocal, elementDisp);
    const localForces = multiplyMatrixVector(kLocal, localDisp);
    const axial = Math.abs(localForces[0]);
    const stress = axial / A;
    memberForces.push({ i, j, axial, stress });
    maxStress = Math.max(maxStress, stress);
  }

  const safetyFactor = yieldMPa / maxStress;
  const forceBalance = problem.loads.reduce((sum, load) => sum + load[1] + load[2], 0.0);
  const reactionBalance = reactions.reduce((sum, value) => sum + value, 0.0);

  console.log(JSON.stringify({
    safetyFactor,
    maxStress,
    yieldMPa,
    loadCount: problem.loads.length,
    fixedCount: fixed.size / 3,
    memberCount: members.length,
    forceBalance,
    reactionBalance,
    maxAxial: Math.max(...memberForces.map((m) => m.axial)),
  }, null, 2));
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
