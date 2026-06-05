// ============================================================================
// NMA Engine v5.3 - With Advanced Methods (Q-profile CI, Contributions, CINeMA, Threshold Analysis)
// ============================================================================

// Named constants for numerical tolerances
const TOLERANCE = {
    SINGULAR: 1e-14,      // Singularity detection
    SVD: 1e-10,           // SVD pseudoinverse cutoff
    CONVERGENCE: 1e-6,    // Iteration convergence
    SMALL_SE: 0.001,      // Warning for very small SE
    LARGE_ES: 5           // Warning for large effect size
};

const Engine = {
    // ========== Matrix Operations ==========
    zeros: (r, c) => Array(r).fill(0).map(() => Array(c).fill(0)),
    
    eye: (n) => {
        const I = Engine.zeros(n, n);
        for (let i = 0; i < n; i++) I[i][i] = 1;
        return I;
    },
    
    matVec: (A, v) => A.map(row => row.reduce((s, a, j) => s + a * v[j], 0)),
    
    matMul: (A, B) => {
        const m = A.length, n = B[0].length, k = B.length;
        const C = Engine.zeros(m, n);
        for (let i = 0; i < m; i++) {
            for (let j = 0; j < n; j++) {
                for (let p = 0; p < k; p++) {
                    C[i][j] += A[i][p] * B[p][j];
                }
            }
        }
        return C;
    },
    
    transpose: (A) => {
        const m = A.length, n = A[0].length;
        const T = Engine.zeros(n, m);
        for (let i = 0; i < m; i++) {
            for (let j = 0; j < n; j++) {
                T[j][i] = A[i][j];
            }
        }
        return T;
    },
    
    diag: (A) => A.map((r, i) => r[i]),
    
    trace: (A) => Engine.diag(A).reduce((s, v) => s + v, 0),
    
    diagMat: (v) => {
        const n = v.length;
        const D = Engine.zeros(n, n);
        for (let i = 0; i < n; i++) D[i][i] = v[i];
        return D;
    },
    
    // Matrix inversion with partial pivoting
    inv: (M) => {
        const n = M.length;
        const A = M.map(r => [...r]);
        const I = Engine.eye(n);
        
        for (let i = 0; i < n; i++) {
            let maxRow = i, maxVal = Math.abs(A[i][i]);
            for (let k = i + 1; k < n; k++) {
                if (Math.abs(A[k][i]) > maxVal) {
                    maxVal = Math.abs(A[k][i]);
                    maxRow = k;
                }
            }
            
            if (maxRow !== i) {
                [A[i], A[maxRow]] = [A[maxRow], A[i]];
                [I[i], I[maxRow]] = [I[maxRow], I[i]];
            }
            
            const pivot = A[i][i];
            if (Math.abs(pivot) < TOLERANCE.SINGULAR) {
                return null;
            }
            
            for (let j = 0; j < n; j++) {
                A[i][j] /= pivot;
                I[i][j] /= pivot;
            }
            
            for (let k = 0; k < n; k++) {
                if (k !== i) {
                    const f = A[k][i];
                    for (let j = 0; j < n; j++) {
                        A[k][j] -= f * A[i][j];
                        I[k][j] -= f * I[i][j];
                    }
                }
            }
        }
        return I;
    },
    
    // SVD via Jacobi rotations
    svd: (A, maxIter = 100, tol = 1e-12) => {
        const m = A.length, n = A[0].length;
        let U = A.map(r => [...r]);
        let V = Engine.eye(n);
        
        for (let iter = 0; iter < maxIter; iter++) {
            let converged = true;
            
            for (let i = 0; i < n - 1; i++) {
                for (let j = i + 1; j < n; j++) {
                    let a = 0, b = 0, c = 0;
                    for (let k = 0; k < m; k++) {
                        a += U[k][i] * U[k][i];
                        b += U[k][j] * U[k][j];
                        c += U[k][i] * U[k][j];
                    }
                    
                    if (Math.abs(c) < tol) continue;
                    converged = false;
                    
                    const zeta = (b - a) / (2 * c);
                    // sign(0) must be +1 here: when column norms are equal (zeta==0) the
                    // correct Jacobi rotation is 45deg (t=+1). Math.sign(0)===0 would give
                    // t=0 (no rotation), leaving anti-parallel columns un-orthogonalised and
                    // producing a 2x-too-large pseudoinverse on rank-deficient 2x2 Laplacians.
                    const sgnZeta = zeta < 0 ? -1 : 1;
                    const t = sgnZeta / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta));
                    const cos = 1 / Math.sqrt(1 + t * t);
                    const sin = cos * t;
                    
                    for (let k = 0; k < m; k++) {
                        const tmp = U[k][i];
                        U[k][i] = cos * tmp - sin * U[k][j];
                        U[k][j] = sin * tmp + cos * U[k][j];
                    }
                    
                    for (let k = 0; k < n; k++) {
                        const tmp = V[k][i];
                        V[k][i] = cos * tmp - sin * V[k][j];
                        V[k][j] = sin * tmp + cos * V[k][j];
                    }
                }
            }
            
            if (converged) break;
        }
        
        const S = [];
        for (let j = 0; j < n; j++) {
            let norm = 0;
            for (let i = 0; i < m; i++) norm += U[i][j] * U[i][j];
            norm = Math.sqrt(norm);
            S.push(norm);
            if (norm > TOLERANCE.SINGULAR) {
                for (let i = 0; i < m; i++) U[i][j] /= norm;
            }
        }
        
        return { U, S, V };
    },
    
    // Moore-Penrose pseudoinverse via SVD
    pinv: (A, tol = TOLERANCE.SVD) => {
        const { U, S, V } = Engine.svd(A);
        const Sinv = S.map(s => Math.abs(s) > tol ? 1 / s : 0);
        return Engine.matMul(Engine.matMul(V, Engine.diagMat(Sinv)), Engine.transpose(U));
    },
    
    // ========== Statistical Functions ==========
    
    pnorm: (z) => {
        const sign = z < 0 ? -1 : 1;
        const x = Math.abs(z) / Math.sqrt(2);
        const t = 1 / (1 + 0.3275911 * x);
        const poly = ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t;
        return 0.5 * (1 + sign * (1 - poly * Math.exp(-x * x)));
    },
    
    qt: (p, df) => {
        if (df <= 0) return 1.96;
        if (df >= 120) return Engine.qnorm(p);
        const z = Engine.qnorm(p);
        const g1 = (z * z * z + z) / 4;
        const g2 = (5 * z * z * z * z * z + 16 * z * z * z + 3 * z) / 96;
        return z + g1 / df + g2 / (df * df);
    },
    
    qnorm: (p) => {
        if (p <= 0) return -Infinity;
        if (p >= 1) return Infinity;
        if (p === 0.5) return 0;
        
        const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
        const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
        const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
        const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
        
        const pL = 0.02425, pH = 1 - pL;
        let q, r;
        
        if (p < pL) {
            q = Math.sqrt(-2 * Math.log(p));
            return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
        } else if (p <= pH) {
            q = p - 0.5;
            r = q * q;
            return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
        } else {
            q = Math.sqrt(-2 * Math.log(1 - p));
            return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
        }
    },
    
    gamma: (z) => {
        const g = 7;
        const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
        
        if (z < 0.5) return Math.PI / (Math.sin(Math.PI * z) * Engine.gamma(1 - z));
        
        z -= 1;
        let x = c[0];
        for (let i = 1; i < g + 2; i++) x += c[i] / (z + i);
        const t = z + g + 0.5;
        return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
    },
    
    pchisq: (x, df) => {
        if (x <= 0 || df <= 0) return 0;
        
        if (df > 100) {
            const z = Math.pow(x / df, 1 / 3) - (1 - 2 / (9 * df));
            return Engine.pnorm(z / Math.sqrt(2 / (9 * df)));
        }
        
        const a = df / 2;
        const zz = x / 2;
        let sum = 0, term = Math.exp(-zz) * Math.pow(zz, a) / Engine.gamma(a + 1);
        
        for (let k = 0; k < 200; k++) {
            sum += term;
            term *= zz / (a + k + 1);
            if (Math.abs(term) < 1e-14) break;
        }
        
        return sum;
    },
    
    // ========== NMA Core Functions ==========
    
    validateContrasts: (cs) => {
        const errors = [];
        const warnings = [];
        
        cs.forEach((c, i) => {
            const row = i + 1;
            if (isNaN(c.es) || !isFinite(c.es)) errors.push(`Row ${row}: Invalid effect size`);
            if (isNaN(c.se) || c.se <= 0) errors.push(`Row ${row}: SE must be positive`);
            if (!c.t1 || !c.t2) errors.push(`Row ${row}: Missing treatment names`);
            if (c.t1 === c.t2) errors.push(`Row ${row}: Same treatment in both arms`);
            if (c.se < TOLERANCE.SMALL_SE) warnings.push(`Row ${row}: Very small SE (${c.se.toFixed(4)})`);
            if (Math.abs(c.es) > TOLERANCE.LARGE_ES) warnings.push(`Row ${row}: Large ES (${c.es.toFixed(2)})`);
        });
        
        return { errors, warnings, valid: errors.length === 0 };
    },
    
    buildB: (cs, ts) => {
        const m = cs.length, n = ts.length;
        const idx = {};
        ts.forEach((t, i) => idx[t] = i);
        const B = Engine.zeros(m, n);
        cs.forEach((c, i) => {
            B[i][idx[c.t1]] = 1;
            B[i][idx[c.t2]] = -1;
        });
        return B;
    },
    
    buildL: (B, W) => {
        const n = B[0].length;
        const L = Engine.zeros(n, n);
        for (let i = 0; i < B.length; i++) {
            for (let j = 0; j < n; j++) {
                for (let k = 0; k < n; k++) {
                    L[j][k] += B[i][j] * W[i] * B[i][k];
                }
            }
        }
        return L;
    },
    
    estimate: (cs, ts, tau2 = 0) => {
        if (!cs.length || ts.length < 2) return null;
        
        const B = Engine.buildB(cs, ts);
        const W = cs.map(c => 1 / (c.se * c.se + tau2));
        const L = Engine.buildL(B, W);
        const Lp = Engine.pinv(L);
        
        if (!Lp) return null;
        
        const y = cs.map(c => c.es);
        const Bwy = ts.map((_, j) => cs.reduce((s, c, i) => s + B[i][j] * W[i] * y[i], 0));
        const theta = Engine.matVec(Lp, Bwy);
        const yhat = Engine.matVec(B, theta);
        const Q = y.reduce((s, yi, i) => s + W[i] * Math.pow(yi - yhat[i], 2), 0);
        
        const H = Engine.matMul(Engine.matMul(B, Lp), Engine.transpose(B));
        for (let i = 0; i < H.length; i++) {
            for (let j = 0; j < H[0].length; j++) {
                H[i][j] *= Math.sqrt(W[i] * W[j]);
            }
        }
        
        return { theta, V: Lp, Q, H, B, W, y, yhat };
    },
    
    tau2REML: (cs, ts, maxIter = 50) => {
        const m = cs.length, n = ts.length;
        const df = m - n + 1;
        
        if (df <= 0) return { tau2: 0, converged: true, iterations: 0, message: 'df ≤ 0' };
        
        let tau2 = 0.01;
        let converged = false;
        let atBoundary = false;
        let prevScore = null;
        let dampingFactor = 1.0;
        
        for (let iter = 0; iter < maxIter; iter++) {
            const est = Engine.estimate(cs, ts, tau2);
            if (!est) return { tau2: Math.max(0, tau2), converged: false, iterations: iter + 1, message: 'Estimation failed' };
            
            const W = cs.map(c => 1 / (c.se * c.se + tau2));
            const BLB = Engine.matMul(Engine.matMul(est.B, est.V), Engine.transpose(est.B));
            
            const P = Engine.zeros(m, m);
            for (let i = 0; i < m; i++) {
                for (let j = 0; j < m; j++) {
                    P[i][j] = (i === j ? W[i] : 0) - W[i] * BLB[i][j] * W[j];
                }
            }
            
            const r = cs.map((c, i) => c.es - est.yhat[i]);
            let rPr = 0;
            for (let i = 0; i < m; i++) {
                for (let j = 0; j < m; j++) {
                    rPr += r[i] * P[i][j] * r[j];
                }
            }
            
            const trP = Engine.trace(P);
            const score = -0.5 * trP + 0.5 * rPr;
            const P2 = Engine.matMul(P, P);
            const info = 0.5 * Engine.trace(P2);
            
            if (Math.abs(info) < TOLERANCE.SINGULAR) break;
            
            if (prevScore !== null && Math.sign(score) !== Math.sign(prevScore)) {
                dampingFactor *= 0.5;
            }
            prevScore = score;
            
            let step = dampingFactor * score / info;
            let newTau2 = tau2 + step;
            
            if (newTau2 < 0) {
                newTau2 = 0;
                atBoundary = true;
            }
            
            if (Math.abs(newTau2 - tau2) < TOLERANCE.CONVERGENCE) {
                tau2 = newTau2;
                converged = true;
                break;
            }
            
            tau2 = newTau2;
        }
        
        return { tau2: Math.max(0, tau2), converged, iterations: maxIter, atBoundary, message: atBoundary ? 'Converged at boundary' : (converged ? 'Converged' : 'Max iterations') };
    },
    
    tau2DL: (cs, ts, maxIter = 30) => {
        const m = cs.length, n = ts.length;
        const df = m - n + 1;
        
        if (df <= 0) return { tau2: 0, converged: true, iterations: 0, message: 'df ≤ 0' };
        
        let tau2 = 0;
        let converged = false;
        
        for (let iter = 0; iter < maxIter; iter++) {
            const est = Engine.estimate(cs, ts, tau2);
            if (!est) return { tau2: Math.max(0, tau2), converged: false, iterations: iter + 1, message: 'Estimation failed' };
            
            const W = cs.map(c => 1 / (c.se * c.se + tau2));
            const H = Engine.matMul(est.B, Engine.matMul(est.V, Engine.transpose(est.B)));
            for (let i = 0; i < H.length; i++) {
                for (let j = 0; j < H[0].length; j++) {
                    H[i][j] *= Math.sqrt(W[i] * W[j]);
                }
            }
            
            const WH = Engine.matMul(H, Engine.diagMat(W));
            const denom = W.reduce((a, b) => a + b, 0) - Engine.trace(WH);
            
            if (denom <= 0) break;
            
            const newTau2 = Math.max(0, (est.Q - df) / denom);
            
            if (Math.abs(newTau2 - tau2) < TOLERANCE.CONVERGENCE) {
                tau2 = newTau2;
                converged = true;
                break;
            }
            
            tau2 = newTau2;
        }
        
        return { tau2: Math.max(0, tau2), converged, iterations: maxIter, message: converged ? 'Converged' : 'Max iterations' };
    },
    
    tau2PM: (cs, ts, maxIter = 100) => {
        const m = cs.length, n = ts.length;
        const df = m - n + 1;
        
        if (df <= 0) return { tau2: 0, converged: true, iterations: 0, message: 'df ≤ 0' };
        
        let tau2 = 0;
        let converged = false;
        
        for (let iter = 0; iter < maxIter; iter++) {
            const est = Engine.estimate(cs, ts, tau2);
            if (!est) return { tau2: Math.max(0, tau2), converged: false, iterations: iter + 1, message: 'Estimation failed' };
            
            if (est.Q <= df) {
                tau2 = 0;
                converged = true;
                break;
            }
            
            const eps = Math.max(1e-8, tau2 * 0.001);
            const est2 = Engine.estimate(cs, ts, tau2 + eps);
            if (!est2) break;
            
            const fprime = (est2.Q - est.Q) / eps;
            if (Math.abs(fprime) < TOLERANCE.SINGULAR) break;
            
            const newTau2 = Math.max(0, tau2 - (est.Q - df) / fprime);
            
            if (Math.abs(newTau2 - tau2) < 1e-8) {
                tau2 = newTau2;
                converged = true;
                break;
            }
            
            tau2 = newTau2;
        }
        
        return { tau2: Math.max(0, tau2), converged, iterations: maxIter, message: converged ? 'Converged' : 'Max iterations' };
    },
    
    computeTau2: (cs, ts, method = 'reml') => {
        switch (method) {
            case 'reml': return Engine.tau2REML(cs, ts);
            case 'dl': return Engine.tau2DL(cs, ts);
            case 'pm': return Engine.tau2PM(cs, ts);
            default: return Engine.tau2REML(cs, ts);
        }
    },
    
    I2: (Q, df) => {
        if (df <= 0 || Q <= 0) return 0;
        return Math.min(100, Math.max(0, 100 * (Q - df) / Q));
    },
    
    pairwise: (theta, V, ts, metric, tau2 = 0, df = null, useHKSJ = false) => {
        const n = ts.length;
        const isRatio = ['OR', 'RR'].includes(metric);
        const results = [];
        const crit = useHKSJ && df && df > 0 ? Engine.qt(0.975, df) : 1.96;
        
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                if (i !== j) {
                    const es = theta[i] - theta[j];
                    const varDiff = V[i][i] + V[j][j] - 2 * V[i][j];
                    const se = varDiff > 0 ? Math.sqrt(varDiff) : 0;
                    
                    let z, pv;
                    if (se > TOLERANCE.SINGULAR) {
                        z = es / se;
                        pv = 2 * (1 - Engine.pnorm(Math.abs(z)));
                    } else {
                        z = es !== 0 ? (es > 0 ? Infinity : -Infinity) : 0;
                        pv = es !== 0 ? 0 : 1;
                    }
                    
                    const predSE = Math.sqrt(se * se + tau2);
                    const predLo = es - crit * predSE;
                    const predHi = es + crit * predSE;
                    
                    results.push({
                        t1: ts[i], t2: ts[j], es, se,
                        lo: es - crit * se, hi: es + crit * se,
                        d: isRatio ? Math.exp(es) : es,
                        dL: isRatio ? Math.exp(es - crit * se) : es - crit * se,
                        dH: isRatio ? Math.exp(es + crit * se) : es + crit * se,
                        predLo: isRatio ? Math.exp(predLo) : predLo,
                        predHi: isRatio ? Math.exp(predHi) : predHi,
                        z, pv, sig: pv < 0.05, crit, seWarning: se < TOLERANCE.SINGULAR
                    });
                }
            }
        }
        
        return results;
    },
    
    pScores: (theta, V, smallBetter) => {
        const n = theta.length;
        return theta.map((_, i) => {
            let score = 0;
            for (let j = 0; j < n; j++) {
                if (i !== j) {
                    const d = theta[i] - theta[j];
                    const varDiff = V[i][i] + V[j][j] - 2 * V[i][j];
                    const se = varDiff > 0 ? Math.sqrt(varDiff) : 0;
                    let p = se > TOLERANCE.SINGULAR ? Engine.pnorm(d / se) : (d > 0 ? 1 : (d < 0 ? 0 : 0.5));
                    if (smallBetter) p = 1 - p;
                    score += p;
                }
            }
            return score / (n - 1);
        });
    },
    
    networkDiagnostics: (cs, ts) => {
        const n = ts.length;
        const adj = {};
        ts.forEach(t => adj[t] = new Set());
        cs.forEach(c => { adj[c.t1].add(c.t2); adj[c.t2].add(c.t1); });
        
        const degrees = ts.map(t => adj[t].size);
        const avgDeg = degrees.reduce((a, b) => a + b, 0) / n;
        const maxEdges = n * (n - 1) / 2;
        const actualEdges = Engine.directComparisons(cs).length;
        const density = maxEdges > 0 ? actualEdges / maxEdges : 0;
        
        const artPts = [];
        const visited = new Set();
        const disc = {}, low = {}, parent = {};
        let time = 0;
        
        const dfs = (u) => {
            let children = 0;
            visited.add(u);
            disc[u] = low[u] = time++;
            
            for (const v of adj[u]) {
                if (!visited.has(v)) {
                    children++;
                    parent[v] = u;
                    dfs(v);
                    low[u] = Math.min(low[u], low[v]);
                    if ((parent[u] === undefined && children > 1) || (parent[u] !== undefined && low[v] >= disc[u])) {
                        if (!artPts.includes(u)) artPts.push(u);
                    }
                } else if (v !== parent[u]) {
                    low[u] = Math.min(low[u], disc[v]);
                }
            }
        };
        
        if (ts.length > 0) dfs(ts[0]);
        return { degrees, avgDegree: avgDeg, density, articulationPoints: artPts };
    },
    
    consistencyTest: (cs, ts, tau2) => {
        const m = cs.length, n = ts.length;
        const designMap = {};
        cs.forEach(c => {
            if (!designMap[c.study]) designMap[c.study] = new Set();
            designMap[c.study].add(c.t1);
            designMap[c.study].add(c.t2);
        });
        
        const designs = {};
        cs.forEach(c => {
            const key = Array.from(designMap[c.study]).sort().join('|');
            if (!designs[key]) designs[key] = [];
            designs[key].push(c);
        });
        
        const designKeys = Object.keys(designs);
        let Qwithin = 0, dfWithin = 0;
        
        designKeys.forEach(key => {
            const dc = designs[key];
            if (dc.length > 1) {
                const dt = Array.from(new Set([...dc.map(c => c.t1), ...dc.map(c => c.t2)])).sort();
                const est = Engine.estimate(dc, dt, tau2);
                if (est) { Qwithin += est.Q; dfWithin += dc.length - dt.length + 1; }
            }
        });
        
        const fullEst = Engine.estimate(cs, ts, tau2);
        if (!fullEst) return null;
        
        const Qtotal = fullEst.Q;
        const dfTotal = m - n + 1;
        const Qbetween = Qtotal - Qwithin;
        const dfBetween = dfTotal - dfWithin;
        const dfWarning = dfBetween < 0;
        const safeDfBetween = Math.max(0, dfBetween);
        
        return {
            Qtotal, dfTotal, pTotal: dfTotal > 0 ? 1 - Engine.pchisq(Qtotal, dfTotal) : 1,
            Qwithin, dfWithin, pWithin: dfWithin > 0 ? 1 - Engine.pchisq(Qwithin, dfWithin) : 1,
            Qbetween: Math.max(0, Qbetween), dfBetween: safeDfBetween,
            pBetween: safeDfBetween > 0 ? 1 - Engine.pchisq(Math.max(0, Qbetween), safeDfBetween) : 1,
            nDesigns: designKeys.length,
            inconsistent: safeDfBetween > 0 && (1 - Engine.pchisq(Math.max(0, Qbetween), safeDfBetween)) < 0.10,
            dfWarning, warningMessage: dfWarning ? 'Negative df - possible model misspecification' : null
        };
    },
    
    nodeSplit: (cs, ts, t1, t2, tau2) => {
        const direct = cs.filter(c => (c.t1 === t1 && c.t2 === t2) || (c.t1 === t2 && c.t2 === t1));
        if (!direct.length) return null;
        
        let sumW = 0, sumWY = 0;
        direct.forEach(c => {
            const w = 1 / (c.se * c.se + tau2);
            const es = c.t1 === t1 ? c.es : -c.es;
            sumW += w;
            sumWY += w * es;
        });
        
        const directES = sumWY / sumW;
        const directSE = Math.sqrt(1 / sumW);
        
        const indirect = cs.filter(c => !((c.t1 === t1 && c.t2 === t2) || (c.t1 === t2 && c.t2 === t1)));
        if (!indirect.length) return null;
        
        const indirectT = Array.from(new Set([...indirect.map(c => c.t1), ...indirect.map(c => c.t2)])).sort();
        if (!indirectT.includes(t1) || !indirectT.includes(t2)) return null;
        if (!Engine.isConnected(indirect, indirectT)) return null;
        
        const result = Engine.estimate(indirect, indirectT, tau2);
        if (!result) return null;
        
        const i1 = indirectT.indexOf(t1);
        const i2 = indirectT.indexOf(t2);
        const indirectES = result.theta[i1] - result.theta[i2];
        const varDiff = result.V[i1][i1] + result.V[i2][i2] - 2 * result.V[i1][i2];
        const indirectSE = varDiff > 0 ? Math.sqrt(varDiff) : 0;
        
        const diff = directES - indirectES;
        const diffSE = Math.sqrt(directSE * directSE + indirectSE * indirectSE);
        const z = diffSE > TOLERANCE.SINGULAR ? diff / diffSE : 0;
        const pv = 2 * (1 - Engine.pnorm(Math.abs(z)));
        
        return { t1, t2, direct: { es: directES, se: directSE, n: direct.length }, indirect: { es: indirectES, se: indirectSE }, diff, diffSE, z, pv, inconsistent: pv < 0.05 };
    },
    
    isConnected: (cs, ts) => {
        if (ts.length <= 1) return true;
        const adj = {};
        ts.forEach(t => adj[t] = new Set());
        cs.forEach(c => { adj[c.t1].add(c.t2); adj[c.t2].add(c.t1); });
        const visited = new Set([ts[0]]);
        const queue = [ts[0]];
        while (queue.length) {
            const node = queue.shift();
            adj[node].forEach(nb => {
                if (!visited.has(nb)) { visited.add(nb); queue.push(nb); }
            });
        }
        return visited.size === ts.length;
    },
    
    directComparisons: (cs) => {
        const seen = new Set();
        const comps = [];
        cs.forEach(c => {
            const key = [c.t1, c.t2].sort().join('|');
            if (!seen.has(key)) { seen.add(key); comps.push({ t1: c.t1, t2: c.t2 }); }
        });
        return comps;
    },
    
    edgeWeights: (cs) => {
        const w = {};
        cs.forEach(c => {
            const key = [c.t1, c.t2].sort().join('|');
            w[key] = (w[key] || 0) + 1;
        });
        return w;
    },
    
    // ========== NEW: High-Impact NMA Additions ==========
    
    // ========== 1. Contribution Matrix (Papakonstantinou et al. 2018) ==========
    // Computes % contribution of each direct comparison to every NMA estimate
    // NOTE: This is a simplified approximation using the hat matrix.
    // The full Papakonstantinou method uses path/stream decomposition for exact contributions.
    // This approximation is reasonable for most networks but may differ for complex topologies.
    contributionMatrix: (cs, ts, tau2 = 0) => {
        const m = cs.length, n = ts.length;
        if (m === 0 || n < 2) return null;
        
        const B = Engine.buildB(cs, ts);
        const W = cs.map(c => 1 / (c.se * c.se + tau2));
        const L = Engine.buildL(B, W);
        const Lp = Engine.pinv(L);
        if (!Lp) return null;
        
        // H = B * L+ * B' is the hat matrix (before weighting)
        const BLp = Engine.matMul(B, Lp);
        const H = Engine.matMul(BLp, Engine.transpose(B));
        
        // Weight the hat matrix: H_w[i,j] = sqrt(W[i]) * H[i,j] * sqrt(W[j])
        const Hw = Engine.zeros(m, m);
        for (let i = 0; i < m; i++) {
            for (let j = 0; j < m; j++) {
                Hw[i][j] = Math.sqrt(W[i]) * H[i][j] * Math.sqrt(W[j]);
            }
        }
        
        // Direct evidence proportion = diagonal of H (own contribution)
        const directProp = Engine.diag(Hw);
        
        // Build contribution matrix for each pairwise comparison
        // For each NMA estimate (t_i vs t_j), compute contribution from each study
        const contributions = {};
        const pairKeys = [];
        
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                if (i !== j) {
                    const key = `${ts[i]}|${ts[j]}`;
                    pairKeys.push({ t1: ts[i], t2: ts[j], key });
                    
                    // Contribution from each study to this comparison
                    const studyContrib = [];
                    let totalContrib = 0;
                    
                    for (let k = 0; k < m; k++) {
                        // Flow from study k to comparison (i,j)
                        // Uses hat matrix row corresponding to (i-j) difference
                        const contrib = Math.abs(BLp[k][i] - BLp[k][j]) * Math.sqrt(W[k]);
                        studyContrib.push({
                            study: cs[k].study,
                            t1: cs[k].t1,
                            t2: cs[k].t2,
                            contrib: contrib
                        });
                        totalContrib += contrib;
                    }
                    
                    // Normalize to percentages
                    if (totalContrib > 0) {
                        studyContrib.forEach(s => {
                            s.pct = (s.contrib / totalContrib) * 100;
                        });
                    }
                    
                    // Sort by contribution
                    studyContrib.sort((a, b) => b.pct - a.pct);
                    contributions[key] = studyContrib;
                }
            }
        }
        
        // Compute overall direct vs indirect proportion for each comparison
        const directIndirect = {};
        for (const pk of pairKeys) {
            const contribs = contributions[pk.key];
            let directSum = 0, indirectSum = 0;
            
            contribs.forEach(c => {
                const isDirect = (c.t1 === pk.t1 && c.t2 === pk.t2) || 
                                (c.t1 === pk.t2 && c.t2 === pk.t1);
                if (isDirect) directSum += c.pct;
                else indirectSum += c.pct;
            });
            
            directIndirect[pk.key] = {
                direct: directSum,
                indirect: indirectSum,
                ratio: directSum / (directSum + indirectSum) || 0
            };
        }
        
        return { 
            H: Hw, 
            directProp, 
            contributions, 
            directIndirect,
            pairKeys 
        };
    },
    
    // ========== 2. Q-Profile CI for τ² (Viechtbauer 2007) ==========
    // Gold standard confidence interval for heterogeneity
    tau2CI_QProfile: (cs, ts, tau2, alpha = 0.05) => {
        const m = cs.length, n = ts.length;
        const df = m - n + 1;
        
        if (df <= 0) return { lower: 0, upper: Infinity, df, method: 'Q-profile' };
        
        // Chi-squared quantiles
        const qLower = Engine.qchisq(alpha / 2, df);
        const qUpper = Engine.qchisq(1 - alpha / 2, df);
        
        // Function to compute Q at given tau2
        const computeQ = (t2) => {
            const est = Engine.estimate(cs, ts, Math.max(0, t2));
            return est ? est.Q : df;  // Return df if estimation fails
        };
        
        // Current Q at point estimate
        const Q0 = computeQ(tau2);
        
        // Robust bisection with bracketing
        const bisect = (targetQ, lo, hi, maxIter = 100) => {
            // Ensure we have valid brackets
            let qLo = computeQ(lo);
            let qHi = computeQ(hi);
            
            // Q decreases as tau2 increases, so qLo > qHi typically
            for (let i = 0; i < maxIter; i++) {
                const mid = (lo + hi) / 2;
                const qMid = computeQ(mid);
                
                if (Math.abs(qMid - targetQ) < 0.0001 || Math.abs(hi - lo) < 1e-8) {
                    return mid;
                }
                
                // Q decreases with increasing tau2
                if (qMid > targetQ) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            return (lo + hi) / 2;
        };
        
        // Lower bound of τ² CI: where Q = qUpper (large Q → small τ²)
        // If Q at τ²=0 is less than qUpper, lower bound is 0
        let tau2Lower = 0;
        const Q_at_zero = computeQ(0);
        if (Q_at_zero > qUpper) {
            // Need to find τ² where Q drops to qUpper
            // Search between 0 and a reasonable upper bound
            const searchHi = Math.max(tau2 * 50, Q_at_zero / df, 1);
            tau2Lower = bisect(qUpper, 0, searchHi);
        }
        
        // Upper bound of τ² CI: where Q = qLower (small Q → large τ²)
        let tau2Upper = Infinity;
        if (Q0 > qLower) {
            // Q can still decrease; find where it reaches qLower
            // Start search from current tau2
            let searchHi = Math.max(tau2 * 100, 10);
            
            // Expand search range if needed
            while (computeQ(searchHi) > qLower && searchHi < 1e6) {
                searchHi *= 2;
            }
            
            if (computeQ(searchHi) <= qLower) {
                tau2Upper = bisect(qLower, tau2, searchHi);
            }
            // If still can't find it, upper bound is truly very large/infinite
        }
        
        return { 
            lower: Math.max(0, tau2Lower), 
            upper: tau2Upper, 
            df, 
            Q: Q0,
            qLower,
            qUpper,
            method: 'Q-profile' 
        };
    },
    
    // Chi-squared quantile function (inverse CDF)
    qchisq: (p, df) => {
        if (p <= 0) return 0;
        if (p >= 1) return Infinity;
        if (df <= 0) return 0;
        
        // Newton-Raphson approximation
        // Initial guess using Wilson-Hilferty approximation
        let x = df * Math.pow(1 - 2/(9*df) + Engine.qnorm(p) * Math.sqrt(2/(9*df)), 3);
        
        // For small df, Wilson-Hilferty can give negative/poor initial guess
        if (df <= 2 || x <= 0) {
            // Use simpler initial guess for small df
            x = Math.max(0.001, df * p);
        } else {
            x = Math.max(0.001, Math.abs(x));
        }
        
        for (let i = 0; i < 30; i++) {
            const cdf = Engine.pchisq(x, df);
            const pdf = Engine.dchisq(x, df);
            if (Math.abs(pdf) < 1e-14) break;
            const delta = (cdf - p) / pdf;
            x = Math.max(0.001, x - delta);
            if (Math.abs(delta) < 1e-8) break;
        }
        
        return x;
    },
    
    // Chi-squared PDF
    dchisq: (x, df) => {
        if (x <= 0 || df <= 0) return 0;
        const k = df / 2;
        return Math.pow(x, k - 1) * Math.exp(-x / 2) / (Math.pow(2, k) * Engine.gamma(k));
    },
    
    // ========== 3. Comparison-Adjusted Funnel Plot (Chaimani et al. 2013) ==========
    // Correct funnel plot for NMA - standard funnel plots don't work
    comparisonAdjustedFunnel: (cs, ts, tau2 = 0, treatmentOrder = null) => {
        // If no treatment order specified, use alphabetical with reference last
        if (!treatmentOrder) {
            treatmentOrder = [...ts].sort();
        }
        
        // Get NMA estimates
        const est = Engine.estimate(cs, ts, tau2);
        if (!est) return null;
        
        // Build lookup for pooled estimates
        const pooledLookup = {};
        for (let i = 0; i < ts.length; i++) {
            for (let j = 0; j < ts.length; j++) {
                if (i !== j) {
                    pooledLookup[`${ts[i]}|${ts[j]}`] = est.theta[i] - est.theta[j];
                }
            }
        }
        
        // For each study, compute centered residual
        const points = cs.map((c, idx) => {
            // Get pooled estimate for this comparison
            const pooled = pooledLookup[`${c.t1}|${c.t2}`] || 0;
            
            // Determine sign based on treatment ordering
            const idx1 = treatmentOrder.indexOf(c.t1);
            const idx2 = treatmentOrder.indexOf(c.t2);
            const sign = idx1 < idx2 ? 1 : -1;
            
            // Centered effect = sign * (observed - pooled)
            const centered = sign * (c.es - pooled);
            
            return {
                x: centered,
                y: c.se,
                study: c.study,
                comparison: `${c.t1} vs ${c.t2}`,
                observed: c.es,
                pooled: pooled,
                residual: c.es - pooled
            };
        });
        
        // Compute Egger's test on centered effects
        const eggerTest = Engine.eggerTest(points);
        
        return {
            points,
            treatmentOrder,
            eggerTest,
            method: 'comparison-adjusted'
        };
    },
    
    // Egger's regression test for funnel asymmetry
    eggerTest: (points) => {
        if (points.length < 3) return { intercept: 0, se: 0, t: 0, p: 1, significant: false };
        
        // Regress (ES/SE) on (1/SE)
        // Equivalent to weighted regression of ES on SE with weights 1/SE²
        const n = points.length;
        let sumX = 0, sumY = 0, sumXX = 0, sumXY = 0;
        
        points.forEach(pt => {
            const x = 1 / pt.y;  // precision
            const y = pt.x / pt.y;  // standardized effect
            sumX += x;
            sumY += y;
            sumXX += x * x;
            sumXY += x * y;
        });
        
        const meanX = sumX / n;
        const meanY = sumY / n;
        const sxx = sumXX - n * meanX * meanX;
        
        if (Math.abs(sxx) < 1e-10) {
            return { intercept: 0, se: 0, t: 0, p: 1, significant: false };
        }
        
        const slope = (sumXY - n * meanX * meanY) / sxx;
        const intercept = meanY - slope * meanX;
        
        // SE of intercept using proper formula
        let ssRes = 0;
        points.forEach(pt => {
            const x = 1 / pt.y;
            const y = pt.x / pt.y;
            const pred = intercept + slope * x;
            ssRes += (y - pred) * (y - pred);
        });
        
        const df = n - 2;
        const mse = df > 0 ? ssRes / df : 0;
        const seIntercept = Math.sqrt(mse * (1/n + meanX * meanX / sxx));
        
        const t = seIntercept > 1e-10 ? intercept / seIntercept : 0;
        
        // Proper t-distribution p-value using beta function approximation
        const p = df > 0 ? Engine.pt(Math.abs(t), df) * 2 : 1;
        
        return {
            intercept,
            se: seIntercept,
            t,
            df,
            p: Math.min(1, p),
            significant: p < 0.1  // Traditional threshold for Egger's test
        };
    },
    
    // T-distribution CDF (two-tailed probability in upper tail)
    pt: (t, df) => {
        if (df <= 0) return 0.5;
        // Use regularized incomplete beta function relationship
        // P(T > t) = 0.5 * I_{df/(df+t²)}(df/2, 1/2) for t > 0
        const x = df / (df + t * t);
        return 0.5 * Engine.betaInc(x, df / 2, 0.5);
    },
    
    // Regularized incomplete beta function I_x(a,b)
    betaInc: (x, a, b) => {
        if (x <= 0) return 0;
        if (x >= 1) return 1;
        
        // Use continued fraction expansion for numerical stability
        const bt = Math.exp(
            Engine.lnGamma(a + b) - Engine.lnGamma(a) - Engine.lnGamma(b) +
            a * Math.log(x) + b * Math.log(1 - x)
        );
        
        // Use symmetry for better convergence
        if (x < (a + 1) / (a + b + 2)) {
            return bt * Engine.betaCF(x, a, b) / a;
        } else {
            return 1 - bt * Engine.betaCF(1 - x, b, a) / b;
        }
    },
    
    // Continued fraction for incomplete beta
    betaCF: (x, a, b) => {
        const maxIter = 100;
        const eps = 1e-10;
        
        let qab = a + b;
        let qap = a + 1;
        let qam = a - 1;
        let c = 1;
        let d = 1 - qab * x / qap;
        if (Math.abs(d) < eps) d = eps;
        d = 1 / d;
        let h = d;
        
        for (let m = 1; m <= maxIter; m++) {
            let m2 = 2 * m;
            let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1 + aa * d;
            if (Math.abs(d) < eps) d = eps;
            c = 1 + aa / c;
            if (Math.abs(c) < eps) c = eps;
            d = 1 / d;
            h *= d * c;
            
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = 1 + aa * d;
            if (Math.abs(d) < eps) d = eps;
            c = 1 + aa / c;
            if (Math.abs(c) < eps) c = eps;
            d = 1 / d;
            let del = d * c;
            h *= del;
            
            if (Math.abs(del - 1) < eps) break;
        }
        return h;
    },
    
    // Log gamma function (Lanczos approximation)
    lnGamma: (z) => {
        const g = 7;
        const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                   771.32342877765313, -176.61502916214059, 12.507343278686905,
                   -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
        
        if (z < 0.5) {
            return Math.log(Math.PI / Math.sin(Math.PI * z)) - Engine.lnGamma(1 - z);
        }
        
        z -= 1;
        let x = c[0];
        for (let i = 1; i < g + 2; i++) {
            x += c[i] / (z + i);
        }
        let t = z + g + 0.5;
        return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
    },
    
    // ========== 4. CINeMA-Style Confidence Rating ==========
    // Semi-automatic 6-domain framework (Nikolakopoulou et al. 2020)
    cinemaRating: (cs, ts, tau2, pairwiseResults, consistencyResults, robInput = null) => {
        const ratings = {};
        const n = ts.length;
        
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                if (i === j) continue;
                
                const key = `${ts[i]}|${ts[j]}`;
                const pw = pairwiseResults.find(p => p.t1 === ts[i] && p.t2 === ts[j]);
                if (!pw) continue;
                
                const domains = {
                    // 1. Within-study bias (from ROB input if available)
                    withinStudyBias: 'not assessed',
                    
                    // 2. Reporting bias (from funnel plot)
                    reportingBias: 'not assessed',
                    
                    // 3. Indirectness (requires user input)
                    indirectness: 'not assessed',
                    
                    // 4. Imprecision (from CI width and clinical threshold)
                    imprecision: Engine.assessImprecision(pw),
                    
                    // 5. Heterogeneity (from τ² and I²)
                    heterogeneity: Engine.assessHeterogeneity(tau2, cs),
                    
                    // 6. Incoherence (from consistency tests)
                    incoherence: consistencyResults ? 
                        (consistencyResults.inconsistent ? 'major concerns' : 'no concerns') : 
                        'not assessed'
                };
                
                // Assess ROB if provided
                if (robInput) {
                    const relevantStudies = cs.filter(c => 
                        (c.t1 === ts[i] && c.t2 === ts[j]) || 
                        (c.t1 === ts[j] && c.t2 === ts[i])
                    );
                    const highRob = relevantStudies.filter(s => s.rob === 'high').length;
                    const someRob = relevantStudies.filter(s => s.rob === 'some').length;
                    const total = relevantStudies.length;
                    
                    if (total > 0) {
                        const highPct = highRob / total;
                        if (highPct > 0.5) domains.withinStudyBias = 'major concerns';
                        else if (highPct > 0.2 || someRob / total > 0.5) domains.withinStudyBias = 'some concerns';
                        else domains.withinStudyBias = 'no concerns';
                    }
                }
                
                // Overall confidence level
                const concerns = Object.values(domains).filter(v => v !== 'not assessed');
                const majorConcerns = concerns.filter(v => v === 'major concerns').length;
                const someConcerns = concerns.filter(v => v === 'some concerns').length;
                
                let confidence;
                if (majorConcerns >= 2) confidence = 'very low';
                else if (majorConcerns === 1 || someConcerns >= 2) confidence = 'low';
                else if (someConcerns === 1) confidence = 'moderate';
                else confidence = 'high';
                
                ratings[key] = {
                    comparison: `${ts[i]} vs ${ts[j]}`,
                    domains,
                    confidence,
                    majorConcerns,
                    someConcerns
                };
            }
        }
        
        return ratings;
    },
    
    // Helper: Assess imprecision domain
    // Based on GRADE guidance: imprecision concerns when CI includes both benefit and harm
    // or when CI is very wide relative to the optimal information size
    assessImprecision: (pw, clinicalThreshold = null) => {
        const ciWidth = pw.hi - pw.lo;
        const absEffect = Math.abs(pw.es);
        
        // Check if CI crosses null (or clinical threshold if provided)
        const threshold = clinicalThreshold !== null ? clinicalThreshold : 0;
        const crossesThreshold = (pw.lo < threshold && pw.hi > threshold);
        
        // Check if CI crosses null AND includes clinically important effects
        // Standard: CI width > 2 times the point estimate suggests high imprecision
        const relativeWidth = absEffect > 0.001 ? ciWidth / absEffect : ciWidth;
        
        // GRADE criteria adapted:
        // - No concerns: Precise estimate, doesn't cross null or threshold
        // - Some concerns: Either crosses null OR wide CI (but not both)
        // - Major concerns: Crosses null AND very wide CI, suggesting high uncertainty
        
        if (!crossesThreshold) {
            // Doesn't cross threshold
            if (relativeWidth < 1.0) return 'no concerns';
            if (relativeWidth < 2.0) return 'some concerns';
            return 'major concerns';  // Very wide CI even without crossing null
        } else {
            // Crosses threshold (includes null or clinical threshold)
            if (relativeWidth < 1.5) return 'some concerns';  // Crosses but relatively precise
            return 'major concerns';  // Crosses and wide
        }
    },
    
    // Helper: Assess heterogeneity domain
    // Based on Cochrane guidance: I² thresholds of 40% and 60%
    assessHeterogeneity: (tau2, cs) => {
        // Compute I² properly from tau2 and average within-study variance
        const m = cs.length;
        const avgVar = cs.reduce((s, c) => s + c.se * c.se, 0) / m;
        const totalVar = avgVar + tau2;
        const I2 = totalVar > 0 ? (tau2 / totalVar) * 100 : 0;
        
        // Cochrane thresholds: 0-40% low, 40-60% moderate, 60%+ substantial
        // But also consider absolute tau2 value
        if (I2 < 40 && tau2 < 0.04) return 'no concerns';
        if (I2 < 60 && tau2 < 0.16) return 'some concerns';
        return 'major concerns';
    },
    
    // ========== 5. Threshold Analysis (Phillippo et al. 2018) ==========
    // How much would evidence need to change before recommendation changes?
    // NOTE: This is a simplified linear approximation. Full Phillippo method uses
    // invariant intervals computed via optimization. This approximation works well
    // for first-order sensitivity but may underestimate thresholds for complex networks.
    thresholdAnalysis: (cs, ts, tau2, decisionThreshold = 0) => {
        const est = Engine.estimate(cs, ts, tau2);
        if (!est) return null;
        
        const n = ts.length;
        const m = cs.length;
        const results = [];
        
        // For each comparison, compute threshold
        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                const es = est.theta[i] - est.theta[j];
                const varDiff = est.V[i][i] + est.V[j][j] - 2 * est.V[i][j];
                const se = Math.sqrt(Math.max(0, varDiff));
                
                // Distance to decision threshold (in log scale)
                const distanceToThreshold = Math.abs(es - decisionThreshold);
                
                // For each study, compute influence on this comparison
                const studyThresholds = cs.map((c, k) => {
                    // Compute partial derivative: ∂(theta_i - theta_j)/∂y_k
                    // Using the hat matrix H = B(B'WB)^+ B'W
                    // Influence = (H[k,i] - H[k,j]) for study k on comparison (i,j)
                    
                    // More accurate: use hat matrix directly
                    let influence = 0;
                    for (let p = 0; p < n; p++) {
                        const Bki = est.B[k][i];
                        const Bkj = est.B[k][j];
                        // Partial derivative through Laplacian pseudoinverse
                        influence += (est.V[i][p] * Bki - est.V[j][p] * Bki) * est.W[k];
                        influence += (est.V[i][p] * (-Bkj) - est.V[j][p] * (-Bkj)) * est.W[k];
                    }
                    
                    // Simplified influence (when above is too complex)
                    // Just use direct contribution from design matrix
                    const simpleInfluence = (est.B[k][i] - est.B[k][j]) * est.W[k];
                    influence = Math.abs(simpleInfluence);
                    
                    // Threshold = distance / influence
                    // How much study ES needs to change to flip decision
                    const threshold = influence > 1e-10 ? 
                        distanceToThreshold / influence : Infinity;
                    
                    return {
                        study: c.study,
                        t1: c.t1,
                        t2: c.t2,
                        originalES: c.es,
                        influence: influence,
                        threshold: threshold,
                        // Threshold in SE units for interpretability
                        thresholdSE: c.se > 0 ? threshold / c.se : Infinity
                    };
                });
                
                // Sort by threshold (smallest = most influential)
                studyThresholds.sort((a, b) => a.threshold - b.threshold);
                
                // Robustness: minimum threshold across all studies
                const minThreshold = Math.min(...studyThresholds.map(s => s.threshold));
                const minThresholdSE = Math.min(...studyThresholds.map(s => s.thresholdSE));
                
                // Robust if requires >1.5 SE change in most influential study
                const isRobust = minThresholdSE > 1.5;
                
                results.push({
                    t1: ts[i],
                    t2: ts[j],
                    estimate: es,
                    se: se,
                    distanceToThreshold,
                    distanceInSE: se > 0 ? distanceToThreshold / se : Infinity,
                    isRobust,
                    minThreshold,
                    minThresholdSE,
                    studyThresholds: studyThresholds.slice(0, 5)  // Top 5 most influential
                });
            }
        }
        
        return {
            decisionThreshold,
            comparisons: results,
            overallRobustness: results.filter(r => r.isRobust).length / results.length
        };
    },
    
    // ========== 6. Network Meta-Regression ==========
    // Allow covariate adjustment
    estimateWithCovariates: (cs, ts, tau2 = 0, covariates = null) => {
        if (!covariates || covariates.length !== cs.length) {
            return Engine.estimate(cs, ts, tau2);
        }
        
        const m = cs.length, n = ts.length;
        const B = Engine.buildB(cs, ts);
        const W = cs.map(c => 1 / (c.se * c.se + tau2));
        
        // Number of covariates
        const p = Array.isArray(covariates[0]) ? covariates[0].length : 1;
        
        // Build extended design matrix [B | X]
        const X = Engine.zeros(m, n + p);
        for (let i = 0; i < m; i++) {
            for (let j = 0; j < n; j++) {
                X[i][j] = B[i][j];
            }
            // Add covariates
            const cov = Array.isArray(covariates[i]) ? covariates[i] : [covariates[i]];
            for (let k = 0; k < p; k++) {
                X[i][n + k] = cov[k] || 0;
            }
        }
        
        // Weighted least squares: (X'WX)^(-1) X'Wy
        const XtW = Engine.zeros(n + p, m);
        for (let i = 0; i < n + p; i++) {
            for (let j = 0; j < m; j++) {
                XtW[i][j] = X[j][i] * W[j];
            }
        }
        
        const XtWX = Engine.matMul(XtW, X);
        const XtWXinv = Engine.pinv(XtWX);
        if (!XtWXinv) return null;
        
        const y = cs.map(c => c.es);
        const XtWy = XtW.map(row => row.reduce((s, v, j) => s + v * y[j], 0));
        const beta = Engine.matVec(XtWXinv, XtWy);
        
        // Split into treatment effects and covariate effects
        const theta = beta.slice(0, n);
        const gamma = beta.slice(n);
        
        // Fitted values
        const yhat = X.map(row => row.reduce((s, v, j) => s + v * beta[j], 0));
        
        // Q statistic
        const Q = y.reduce((s, yi, i) => s + W[i] * Math.pow(yi - yhat[i], 2), 0);
        
        // R² (proportion of heterogeneity explained)
        const estBase = Engine.estimate(cs, ts, tau2);
        const R2 = estBase ? Math.max(0, 1 - Q / estBase.Q) : 0;
        
        return {
            theta,
            gamma,
            V: XtWXinv,
            Q,
            R2,
            covariateNames: covariates.names || gamma.map((_, i) => `Covariate ${i + 1}`)
        };
    },
    
    // ========== 7. SUCRA (alongside P-scores) ==========
    // Surface Under Cumulative Ranking curve
    SUCRA: (theta, V, smallBetter) => {
        const n = theta.length;
        const pscores = Engine.pScores(theta, V, smallBetter);
        
        // Compute rank probabilities via simulation with proper covariance
        const nSim = 2000;
        const rankCounts = Array(n).fill(null).map(() => Array(n).fill(0));
        
        // Cholesky decomposition of V for correlated sampling
        const L = Engine.cholesky(V);
        
        for (let sim = 0; sim < nSim; sim++) {
            // Generate independent standard normals
            const z = Array(n).fill(0).map(() => Engine.rnorm());
            
            // Transform to multivariate normal: samples = theta + L * z
            const samples = theta.map((t, i) => {
                let sum = t;
                for (let j = 0; j <= i; j++) {
                    sum += L[i][j] * z[j];
                }
                return sum;
            });
            
            // Rank (0 = best)
            const indexed = samples.map((s, i) => ({ i, s }));
            indexed.sort((a, b) => smallBetter ? a.s - b.s : b.s - a.s);
            indexed.forEach((item, rank) => rankCounts[item.i][rank]++);
        }
        
        // Compute cumulative rank probabilities and SUCRA
        // SUCRA = (1/(n-1)) * sum_{r=1}^{n-1} P(rank <= r)
        const sucra = rankCounts.map((counts) => {
            const probs = counts.map(c => c / nSim);
            let cumSum = 0;
            let sucraSum = 0;
            for (let r = 0; r < n - 1; r++) {
                cumSum += probs[r];  // P(rank <= r)
                sucraSum += cumSum;
            }
            return sucraSum / (n - 1);
        });
        
        return {
            pscores,
            sucra,
            rankProbabilities: rankCounts.map(c => c.map(v => v / nSim)),
            agreement: pscores.map((p, i) => Math.abs(p - sucra[i]) < 0.05).every(v => v)
        };
    },
    
    // Cholesky decomposition (returns lower triangular L where V = L * L')
    cholesky: (A) => {
        const n = A.length;
        const L = Engine.zeros(n, n);
        
        for (let i = 0; i < n; i++) {
            for (let j = 0; j <= i; j++) {
                let sum = A[i][j];
                for (let k = 0; k < j; k++) {
                    sum -= L[i][k] * L[j][k];
                }
                if (i === j) {
                    // Diagonal: handle near-zero or negative (numerical issues)
                    L[i][j] = sum > 1e-10 ? Math.sqrt(sum) : 1e-10;
                } else {
                    L[i][j] = L[j][j] > 1e-10 ? sum / L[j][j] : 0;
                }
            }
        }
        return L;
    },
    
    // Standard normal random number (Box-Muller)
    rnorm: () => {
        let u1 = Math.random(), u2 = Math.random();
        while (u1 === 0) u1 = Math.random();
        return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    },
    
    // ========== 8. League Table with Contribution Shading ==========
    leagueTableWithContributions: (ts, pairwise, contributions) => {
        const n = ts.length;
        const table = Engine.zeros(n, n);
        const shading = Engine.zeros(n, n);  // 0-1 for direct proportion
        
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                if (i === j) {
                    table[i][j] = null;
                    shading[i][j] = null;
                } else {
                    const pw = pairwise.find(p => p.t1 === ts[i] && p.t2 === ts[j]);
                    table[i][j] = pw ? {
                        es: pw.d,
                        lo: pw.dL,
                        hi: pw.dH,
                        sig: pw.sig
                    } : null;
                    
                    // Shading based on direct evidence proportion
                    const key = `${ts[i]}|${ts[j]}`;
                    if (contributions && contributions.directIndirect && contributions.directIndirect[key]) {
                        shading[i][j] = contributions.directIndirect[key].ratio;
                    } else {
                        shading[i][j] = 0;
                    }
                }
            }
        }
        
        return { table, shading, treatments: ts };
    },
    
    // ========== 9. Study Importance (Rücker et al. 2020) ==========
    // Variance reduction when adding each study (0-1, 1=essential)
    // Optimized version with connectivity caching
    studyImportance: (cs, ts, tau2 = 0) => {
        const fullEst = Engine.estimate(cs, ts, tau2);
        if (!fullEst) return null;
        
        const n = cs.length;
        const importance = [];
        
        // Pre-compute pairwise variance lookup for all comparisons
        const fullVarLookup = {};
        for (let i = 0; i < ts.length; i++) {
            for (let j = i + 1; j < ts.length; j++) {
                const varDiff = fullEst.V[i][i] + fullEst.V[j][j] - 2 * fullEst.V[i][j];
                fullVarLookup[`${ts[i]}|${ts[j]}`] = varDiff;
            }
        }
        
        for (let idx = 0; idx < n; idx++) {
            const c = cs[idx];
            
            // Remove this study
            const reduced = cs.filter((_, i) => i !== idx);
            const reducedTs = [...new Set([...reduced.map(r => r.t1), ...reduced.map(r => r.t2)])];
            
            // Check if treatments still in network
            if (!reducedTs.includes(c.t1) || !reducedTs.includes(c.t2)) {
                // Study introduces new treatments - highly important
                importance.push({
                    study: c.study,
                    t1: c.t1,
                    t2: c.t2,
                    importance: 0.95,
                    critical: false,
                    reason: 'Introduces treatments to network'
                });
                continue;
            }
            
            // Check connectivity
            if (!Engine.isConnected(reduced, reducedTs)) {
                importance.push({
                    study: c.study,
                    t1: c.t1,
                    t2: c.t2,
                    importance: 1.0,
                    critical: true,
                    reason: 'Network disconnects without this study'
                });
                continue;
            }
            
            // Re-estimate without this study
            const reducedEst = Engine.estimate(reduced, reducedTs, tau2);
            if (!reducedEst) {
                importance.push({
                    study: c.study,
                    t1: c.t1,
                    t2: c.t2,
                    importance: 0.9,
                    critical: false,
                    reason: 'Estimation fails without this study'
                });
                continue;
            }
            
            // Compute variance change for the direct comparison
            const key = [c.t1, c.t2].sort().join('|');
            const i1 = ts.indexOf(c.t1);
            const i2 = ts.indexOf(c.t2);
            const ri1 = reducedTs.indexOf(c.t1);
            const ri2 = reducedTs.indexOf(c.t2);
            
            if (ri1 < 0 || ri2 < 0 || i1 < 0 || i2 < 0) {
                importance.push({
                    study: c.study,
                    t1: c.t1,
                    t2: c.t2,
                    importance: 0.5,
                    critical: false,
                    reason: 'Index mismatch'
                });
                continue;
            }
            
            const fullVar = fullVarLookup[key] || 
                (fullEst.V[i1][i1] + fullEst.V[i2][i2] - 2 * fullEst.V[i1][i2]);
            const reducedVar = reducedEst.V[ri1][ri1] + reducedEst.V[ri2][ri2] - 
                2 * reducedEst.V[ri1][ri2];
            
            // Importance = relative variance increase
            // If variance doubles when removed, importance = 0.5
            // If variance goes to infinity, importance = 1.0
            const varRatio = reducedVar / fullVar;
            const imp = Math.min(1, Math.max(0, 1 - 1 / varRatio));
            
            importance.push({
                study: c.study,
                t1: c.t1,
                t2: c.t2,
                importance: imp,
                varianceIncrease: varRatio - 1,
                varianceIncreasePercent: ((varRatio - 1) * 100).toFixed(1),
                critical: false,
                reason: varRatio > 2 ? 'High variance increase' : 'Moderate contribution'
            });
        }
        
        // Sort by importance (descending)
        importance.sort((a, b) => b.importance - a.importance);
        
        return importance;
    }
};

// ============================================================================
// Validation Module
// ============================================================================

const Validation = {
    tests: [],
    
    runAll() {
        this.tests = [];
        
        // Test 1: Matrix Operations
        this.tests.push(this.testMatrixOps());
        
        // Test 2: Triangle Network
        this.tests.push(this.testTriangleNetwork());
        
        // Test 3: REML Convergence
        this.tests.push(this.testREML());
        
        // Test 4: DL Estimation
        this.tests.push(this.testDL());
        
        // Test 5: P-scores
        this.tests.push(this.testPScores());
        
        // Test 6: I² Statistic
        this.tests.push(this.testI2());
        
        // Test 7: Confidence Intervals
        this.tests.push(this.testCIs());
        
        // Test 8: R Comparison
        this.tests.push(this.testRComparison());
        
        // Test 9: Network Connectivity
        this.tests.push(this.testConnectivity());
        
        // Test 10: Numerical Stability
        this.tests.push(this.testNumericalStability());
        
        return this.tests;
    },
    
    testMatrixOps() {
        const L = [[2, -1, -1], [-1, 2, -1], [-1, -1, 2]];
        const Lp = Engine.pinv(L);
        const rowSums = Lp.map(row => row.reduce((a, b) => a + b, 0));
        const pass = rowSums.every(s => Math.abs(s) < 1e-10);
        
        return {
            name: 'Matrix Operations',
            description: 'Laplacian pseudoinverse row sums ≈ 0',
            pass,
            details: `Row sums: [${rowSums.map(s => s.toExponential(2)).join(', ')}]`
        };
    },
    
    testTriangleNetwork() {
        const contrasts = [
            { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.1 },
            { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 0.1 },
            { study: 'S3', t1: 'A', t2: 'C', es: 0.8, se: 0.1 }
        ];
        const treatments = ['A', 'B', 'C'];
        const est = Engine.estimate(contrasts, treatments, 0);
        
        const theta = est.theta;
        const aVsB = theta[0] - theta[1];
        const bVsC = theta[1] - theta[2];
        const aVsC = theta[0] - theta[2];
        const transitivityError = Math.abs((aVsB + bVsC) - aVsC);
        const pass = transitivityError < 1e-10;
        
        return {
            name: 'Triangle Network',
            description: 'Transitivity: A-B + B-C = A-C',
            pass,
            details: `Direct A-C: ${aVsC.toFixed(4)}, Indirect: ${(aVsB + bVsC).toFixed(4)}`
        };
    },
    
    testREML() {
        const contrasts = [
            { study: 'S1', t1: 'A', t2: 'B', es: 0.50, se: 0.40 },
            { study: 'S2', t1: 'A', t2: 'C', es: 0.45, se: 0.35 },
            { study: 'S3', t1: 'B', t2: 'C', es: -0.10, se: 0.30 },
            { study: 'S4', t1: 'A', t2: 'D', es: 0.80, se: 0.45 },
            { study: 'S5', t1: 'C', t2: 'D', es: 0.25, se: 0.38 }
        ];
        const treatments = ['A', 'B', 'C', 'D'];
        const result = Engine.tau2REML(contrasts, treatments);
        const pass = result.converged && result.tau2 >= 0;
        
        return {
            name: 'REML τ² Estimation',
            description: 'Converges with non-negative τ²',
            pass,
            details: `τ² = ${result.tau2.toFixed(6)}, converged: ${result.converged}`
        };
    },
    
    testDL() {
        const contrasts = [
            { study: 'S1', t1: 'A', t2: 'B', es: 0.50, se: 0.40 },
            { study: 'S2', t1: 'A', t2: 'C', es: 0.45, se: 0.35 },
            { study: 'S3', t1: 'B', t2: 'C', es: -0.10, se: 0.30 }
        ];
        const treatments = ['A', 'B', 'C'];
        const reml = Engine.tau2REML(contrasts, treatments);
        const dl = Engine.tau2DL(contrasts, treatments);
        const diff = Math.abs(reml.tau2 - dl.tau2);
        const pass = diff < 0.1 || (reml.tau2 < 0.01 && dl.tau2 < 0.01);
        
        return {
            name: 'DL τ² Estimation',
            description: 'DL ≈ REML within tolerance',
            pass,
            details: `REML: ${reml.tau2.toFixed(6)}, DL: ${dl.tau2.toFixed(6)}`
        };
    },
    
    testPScores() {
        const contrasts = [
            { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.1 },
            { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 0.1 },
            { study: 'S3', t1: 'A', t2: 'C', es: 0.8, se: 0.1 }
        ];
        const treatments = ['A', 'B', 'C'];
        const est = Engine.estimate(contrasts, treatments, 0);
        const pscores = Engine.pScores(est.theta, est.V, false);
        const sum = pscores.reduce((a, b) => a + b, 0);
        const expectedSum = treatments.length * 0.5;
        const inRange = pscores.every(p => p >= 0 && p <= 1);
        const pass = Math.abs(sum - expectedSum) < 0.01 && inRange;
        
        return {
            name: 'P-scores Calculation',
            description: 'Sum = n×0.5, all in [0,1]',
            pass,
            details: `Sum: ${sum.toFixed(4)} (expected ${expectedSum}), range OK: ${inRange}`
        };
    },
    
    testI2() {
        const cases = [
            { Q: 10, df: 5, expected: 50 },
            { Q: 5, df: 5, expected: 0 },
            { Q: 3, df: 5, expected: 0 },
            { Q: 100, df: 10, expected: 90 }
        ];
        
        const results = cases.map(c => {
            const computed = Engine.I2(c.Q, c.df);
            return Math.abs(computed - c.expected) < 0.1;
        });
        
        const pass = results.every(r => r);
        
        return {
            name: 'I² Statistic',
            description: 'Correct I² calculation',
            pass,
            details: `All ${cases.length} test cases passed: ${pass}`
        };
    },
    
    testCIs() {
        const contrasts = [
            { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.1 },
            { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 0.1 },
            { study: 'S3', t1: 'A', t2: 'C', es: 0.8, se: 0.1 }
        ];
        const treatments = ['A', 'B', 'C'];
        const est = Engine.estimate(contrasts, treatments, 0);
        const df = contrasts.length - treatments.length + 1;
        
        const wald = Engine.pairwise(est.theta, est.V, treatments, 'MD', 0, df, false);
        const hksj = Engine.pairwise(est.theta, est.V, treatments, 'MD', 0, df, true);
        
        const waldAB = wald.find(p => p.t1 === 'A' && p.t2 === 'B');
        const hksjAB = hksj.find(p => p.t1 === 'A' && p.t2 === 'B');
        
        const waldWidth = waldAB.hi - waldAB.lo;
        const hksjWidth = hksjAB.hi - hksjAB.lo;
        const pass = hksjWidth >= waldWidth - 0.001;
        
        return {
            name: 'Confidence Intervals',
            description: 'HKSJ CIs ≥ Wald CIs',
            pass,
            details: `Wald width: ${waldWidth.toFixed(3)}, HKSJ width: ${hksjWidth.toFixed(3)}`
        };
    },
    
    testRComparison() {
        const contrasts = [
            { study: 'S1', t1: 'A', t2: 'B', es: 0.50, se: 0.40 },
            { study: 'S2', t1: 'A', t2: 'C', es: 0.45, se: 0.35 },
            { study: 'S3', t1: 'B', t2: 'C', es: -0.10, se: 0.30 },
            { study: 'S4', t1: 'A', t2: 'D', es: 0.80, se: 0.45 },
            { study: 'S5', t1: 'C', t2: 'D', es: 0.25, se: 0.38 }
        ];
        const treatments = ['A', 'B', 'C', 'D'];
        const tau = Engine.tau2REML(contrasts, treatments);
        const est = Engine.estimate(contrasts, treatments, tau.tau2);
        const pscores = Engine.pScores(est.theta, est.V, false);
        const pscoreSum = pscores.reduce((a, b) => a + b, 0);
        
        // R netmeta gives: tau2≈0, Q≈0.038, pscores sum≈2.0
        const tau2Match = tau.tau2 < 0.01;
        const qMatch = Math.abs(est.Q - 0.038) < 0.01;
        const pscoreMatch = Math.abs(pscoreSum - 2.0) < 0.01;
        const pass = tau2Match && qMatch && pscoreMatch;
        
        return {
            name: 'R netmeta Comparison',
            description: 'Matches R netmeta output',
            pass,
            details: `τ²: ${tau.tau2.toFixed(4)}, Q: ${est.Q.toFixed(4)}, P-score sum: ${pscoreSum.toFixed(4)}`
        };
    },
    
    testConnectivity() {
        const connected = [
            { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.1 },
            { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 0.1 }
        ];
        const disconnected = [
            { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.1 },
            { study: 'S2', t1: 'C', t2: 'D', es: 0.3, se: 0.1 }
        ];
        
        const test1 = Engine.isConnected(connected, ['A', 'B', 'C']);
        const test2 = !Engine.isConnected(disconnected, ['A', 'B', 'C', 'D']);
        const pass = test1 && test2;
        
        return {
            name: 'Network Connectivity',
            description: 'Detects connected/disconnected networks',
            pass,
            details: `Connected detected: ${test1}, Disconnected detected: ${test2}`
        };
    },
    
    testNumericalStability() {
        const smallSE = [
            { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.001 },
            { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 0.001 },
            { study: 'S3', t1: 'A', t2: 'C', es: 0.8, se: 0.001 }
        ];
        const largeSE = [
            { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 10.0 },
            { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 10.0 },
            { study: 'S3', t1: 'A', t2: 'C', es: 0.8, se: 10.0 }
        ];
        
        const treatments = ['A', 'B', 'C'];
        const estSmall = Engine.estimate(smallSE, treatments, 0);
        const estLarge = Engine.estimate(largeSE, treatments, 0);
        
        const smallOK = estSmall !== null && estSmall.theta.every(t => isFinite(t));
        const largeOK = estLarge !== null && estLarge.theta.every(t => isFinite(t));
        const pass = smallOK && largeOK;
        
        return {
            name: 'Numerical Stability',
            description: 'Handles extreme SE values',
            pass,
            details: `Small SE: ${smallOK ? 'OK' : 'FAIL'}, Large SE: ${largeOK ? 'OK' : 'FAIL'}`
        };
    }
};

if (typeof module!=='undefined'&&module.exports){ module.exports = { Engine, Validation, TOLERANCE }; }
