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
  return matrix.map((row) => row.reduce((sum, value, idx) => sum + value * vector[idx], 0.0));
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
    if (Math.abs(A[pivot][i]) < 1e-14) {
      throw new Error('Singular matrix');
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

function evaluateStructure(problem, members) {
  const E = 200e9;
  const rodSide = problem.rodSide;
  const A = rodSide * rodSide;
  const I = Math.pow(rodSide, 4) / 12;
  const yieldMPa = problem.yieldMPa * 1e6;
  const nodes = problem.nodes;
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

  let u;
  try {
    u = solveLinearSystem(Kff, Ff);
  } catch (err) {
    return { feasible: false, singular: true };
  }

  const disp = Array(dofCount).fill(0.0);
  for (let i = 0; i < free.length; i++) {
    disp[free[i]] = u[i];
  }

  const reactions = Array(dofCount).fill(0.0);
  for (let i = 0; i < dofCount; i++) {
    reactions[i] = K[i].reduce((sum, value, idx) => sum + value * disp[idx], 0.0) - F[i];
  }

  let maxStress = 0.0;
  let maxAxial = 0.0;
  const memberData = [];

  for (const [i, j] of members) {
    const [x1, y1] = nodes[i];
    const [x2, y2] = nodes[j];
    const dx = x2 - x1;
    const dy = y2 - y1;
    const L = Math.hypot(dx, dy);
    const c = dx / L;
    const s = dy / L;

    const localDisp = [
      c * disp[i * 3] + s * disp[i * 3 + 1],
      -s * disp[i * 3] + c * disp[i * 3 + 1],
      disp[i * 3 + 2],
      c * disp[j * 3] + s * disp[j * 3 + 1],
      -s * disp[j * 3] + c * disp[j * 3 + 1],
      disp[j * 3 + 2],
    ];

    const kLocal = [
      [E * A / L, 0.0, 0.0, -E * A / L, 0.0, 0.0],
      [0.0, 12 * E * I / Math.pow(L, 3), 6 * E * I / Math.pow(L, 2), 0.0, -12 * E * I / Math.pow(L, 3), 6 * E * I / Math.pow(L, 2)],
      [0.0, 6 * E * I / Math.pow(L, 2), 4 * E * I / L, 0.0, -6 * E * I / Math.pow(L, 2), 2 * E * I / L],
      [-E * A / L, 0.0, 0.0, E * A / L, 0.0, 0.0],
      [0.0, -12 * E * I / Math.pow(L, 3), -6 * E * I / Math.pow(L, 2), 0.0, 12 * E * I / Math.pow(L, 3), -6 * E * I / Math.pow(L, 2)],
      [0.0, 6 * E * I / Math.pow(L, 2), 2 * E * I / L, 0.0, -6 * E * I / Math.pow(L, 2), 4 * E * I / L],
    ];

    const localForces = multiplyMatrixVector(kLocal, localDisp);
    const axial = localForces[0];
    const stress = Math.abs(axial) / A;
    maxStress = Math.max(maxStress, stress);
    maxAxial = Math.max(maxAxial, Math.abs(axial));

    memberData.push({ i, j, axial, stress });
  }

  const safetyFactor = yieldMPa / maxStress;
  const totalFx = reactions.filter((_, idx) => idx % 3 === 0).reduce((sum, val) => sum + val, 0.0);
  const totalFy = reactions.filter((_, idx) => idx % 3 === 1).reduce((sum, val) => sum + val, 0.0);

  return {
    feasible: true,
    singular: false,
    safetyFactor,
    maxStress,
    maxAxial,
    reactions,
    memberData,
    totalFx,
    totalFy,
  };
}

async function main() {
  const problem = await fetchProblem();

  let currentMembers = problem.members.slice();
  let currentEval = evaluateStructure(problem, currentMembers);

  if (!currentEval.feasible) {
    throw new Error('Initial structure is infeasible');
  }

  const removals = [];

  while (true) {
    const sortedMembers = currentMembers
      .map((member, index) => ({ member, index }))
      .map(({ member, index }) => {
        const m = currentEval.memberData.find((d) => d.i === member[0] && d.j === member[1]);
        return {
          member,
          index,
          axial: m ? Math.abs(m.axial) : 0.0,
          stress: m ? m.stress : 0.0,
        };
      })
      .sort((a, b) => a.stress - b.stress || a.axial - b.axial);

    let chosen = null;
    let chosenEval = null;

    for (const candidate of sortedMembers) {
      const trialMembers = currentMembers.filter((_, idx) => idx !== candidate.index);
      const trialEval = evaluateStructure(problem, trialMembers);
      if (!trialEval.feasible) {
        continue;
      }
      if (trialEval.safetyFactor >= 1.0) {
        chosen = candidate;
        chosenEval = trialEval;
        break;
      }
    }

    if (!chosen) {
      break;
    }

    currentMembers = currentMembers.filter((_, idx) => idx !== chosen.index);
    currentEval = chosenEval;
    removals.push({ removed: chosen.member, safetyFactor: chosenEval.safetyFactor, maxStress: chosenEval.maxStress });
  }

  console.log(JSON.stringify({
    score: currentMembers.length,
    members: currentMembers.length,
    removed: removals.length,
    safetyFactor: currentEval.safetyFactor,
    maxStress: currentEval.maxStress,
    maxAxial: currentEval.maxAxial,
    totalFx: currentEval.totalFx,
    totalFy: currentEval.totalFy,
    removals,
  }, null, 2));
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
