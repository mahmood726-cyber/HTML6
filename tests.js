// ============================================================================
// tests.js - Pure Node test harness for engine.js (NMA Frequentist engine)
// Run: node tests.js   ->  prints checks, exits 0 if all pass, 1 otherwise.
// Expected values hand-derived independently (see comments / report).
// ============================================================================

const { Engine, Validation } = require('./engine.js');

let passed = 0, failed = 0;
function ok(name, cond, extra) {
    if (cond) { passed++; console.log('  PASS  ' + name + (extra ? '  (' + extra + ')' : '')); }
    else { failed++; console.log('  FAIL  ' + name + (extra ? '  (' + extra + ')' : '')); }
}
function near(a, b, tol) { return Math.abs(a - b) < (tol === undefined ? 1e-9 : tol); }

console.log('=== Matrix operations ===');

// --- 1. Matrix inverse: invert a known 2x2 and verify A*Ainv = I ---
// A = [[4,7],[2,6]], det = 4*6 - 7*2 = 10, Ainv = (1/10)[[6,-7],[-2,4]]
(() => {
    const A = [[4, 7], [2, 6]];
    const Ai = Engine.inv(A);
    ok('inv 2x2 matches closed form [[0.6,-0.7],[-0.2,0.4]]',
        near(Ai[0][0], 0.6) && near(Ai[0][1], -0.7) && near(Ai[1][0], -0.2) && near(Ai[1][1], 0.4),
        `[[${Ai[0][0]},${Ai[0][1]}],[${Ai[1][0]},${Ai[1][1]}]]`);
    const P = Engine.matMul(A, Ai); // A * Ainv should be I
    ok('A * inv(A) == I (2x2)',
        near(P[0][0], 1) && near(P[1][1], 1) && near(P[0][1], 0) && near(P[1][0], 0));
})();

// --- 2. Matrix inverse 3x3: A*Ainv = I (partial-pivoting robustness) ---
(() => {
    const A = [[2, 0, 1], [1, 3, 2], [1, 0, 2]]; // det = 2*(3*2-0) -0 +1*(0-3) = 12-3 = 9
    const Ai = Engine.inv(A);
    const P = Engine.matMul(A, Ai);
    let isI = true;
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) if (!near(P[i][j], i === j ? 1 : 0, 1e-9)) isI = false;
    ok('A * inv(A) == I (3x3)', isI);
    // Singular matrix should return null
    const Sing = [[1, 2], [2, 4]];
    ok('inv(singular) === null', Engine.inv(Sing) === null);
})();

// --- 3. Pseudoinverse of complete-graph Laplacian: row sums == 0 ---
(() => {
    const M = [[2, -1, -1], [-1, 2, -1], [-1, -1, 2]];
    const Mp = Engine.pinv(M);
    const rowSums = Mp.map(r => r.reduce((a, b) => a + b, 0));
    ok('pinv(Laplacian) row sums ~ 0', rowSums.every(s => near(s, 0, 1e-8)),
        rowSums.map(s => s.toExponential(1)).join(','));
    // For this M, M^2 = 3M so M+ = M/9 (true Moore-Penrose). Check M+[0][0] == 2/9.
    // (Pre-fix the buggy SVD returned 1/3 here due to equal-singular-value degeneracy.)
    ok('pinv(complete-graph L)[0][0] == 2/9', near(Mp[0][0], 2 / 9, 1e-8), Mp[0][0].toFixed(6));
    // Defining Moore-Penrose property on a rank-1 2x2 Laplacian: A*A+*A == A.
    const R = [[1, -1], [-1, 1]];
    const Rp = Engine.pinv(R);
    const RpR = Engine.matMul(Engine.matMul(R, Rp), R);
    ok('pinv: A*A+*A == A on rank-deficient [[1,-1],[-1,1]] (SVD degeneracy fix)',
        near(RpR[0][0], 1, 1e-8) && near(RpR[0][1], -1, 1e-8),
        `Rp[0][0]=${Rp[0][0].toFixed(4)} (expect 0.25)`);
})();

console.log('=== Statistical primitives ===');

// --- 4. Normal CDF: pnorm(0) ~ 0.5 (NOT 0), monotone, tails ---
(() => {
    ok('pnorm(0) ~ 0.5', near(Engine.pnorm(0), 0.5, 1e-3), Engine.pnorm(0).toFixed(5));
    ok('pnorm(1.959964) ~ 0.975', near(Engine.pnorm(1.959964), 0.975, 1e-3), Engine.pnorm(1.959964).toFixed(5));
    ok('pnorm(-z) == 1 - pnorm(z) (symmetry)', near(Engine.pnorm(-1.3) + Engine.pnorm(1.3), 1, 1e-6));
})();

// --- 5. qnorm inverse of pnorm ---
(() => {
    ok('qnorm(0.975) ~ 1.95996', near(Engine.qnorm(0.975), 1.959964, 1e-3), Engine.qnorm(0.975).toFixed(5));
    ok('qnorm(0.5) == 0', near(Engine.qnorm(0.5), 0));
})();

console.log('=== Hand-worked triangle network (tau2 = 0) ===');
// Contrasts: A-B=0.5, B-C=0.3, A-C=0.8, all se=0.1 (w=100). Treatments [A,B,C].
// L = 100*[[2,-1,-1],[-1,2,-1],[-1,-1,2]]; B'Wy = 100*[1.3,-0.2,-1.1].
// M+ = M/9 (since M^2=3M), so theta = (1/9)*M*[1.3,-0.2,-1.1] = [0.43333,-0.06667,-0.36667].
// Contrast variances: V = M/900; var(A-B)= (2+2-2*(-1))/900 = 6/900 = 0.0066667; se = 0.081650.
(() => {
    const cs = [
        { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.1 },
        { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 0.1 },
        { study: 'S3', t1: 'A', t2: 'C', es: 0.8, se: 0.1 }
    ];
    const ts = ['A', 'B', 'C'];
    const est = Engine.estimate(cs, ts, 0);
    const th = est.theta;
    ok('theta_A == 0.433333 (hand-derived)', near(th[0], 13 / 30, 1e-6), th[0].toFixed(6));
    ok('theta_B == -0.066667 (hand-derived)', near(th[1], -2 / 30, 1e-6), th[1].toFixed(6));
    ok('theta_C == -0.366667 (hand-derived)', near(th[2], -11 / 30, 1e-6), th[2].toFixed(6));
    ok('sum(theta) == 0 (centering)', near(th[0] + th[1] + th[2], 0, 1e-9));
    // Consistency / transitivity in a consistent network: theta_A-C == (A-B)+(B-C)
    const ac = th[0] - th[2], abPlusBc = (th[0] - th[1]) + (th[1] - th[2]);
    ok('consistency theta_AC == theta_AB + theta_BC', near(ac, abPlusBc, 1e-12), ac.toFixed(6));
    // Pooled contrast recovers inputs exactly (over-determined consistent net)
    ok('pooled A-B == 0.5 (input recovered)', near(th[0] - th[1], 0.5, 1e-9));
    ok('pooled A-C == 0.8 (input recovered)', near(th[0] - th[2], 0.8, 1e-9));
    // Hand-derived contrast variance
    const vAB = est.V[0][0] + est.V[1][1] - 2 * est.V[0][1];
    ok('var(A-B) == 6/900 = 0.0066667 (hand-derived)', near(vAB, 6 / 900, 1e-7), vAB.toFixed(7));
    // V symmetric
    let sym = true;
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) if (!near(est.V[i][j], est.V[j][i], 1e-9)) sym = false;
    ok('V (variance-covariance) is symmetric', sym);
})();

console.log('=== Two 2-arm studies, single comparison pooling ===');
// Two A-B studies: es 0.4 se 0.2 (w=25), es 0.6 se 0.1 (w=100). tau2=0.
// IV pooled = (25*0.4 + 100*0.6)/125 = (10+60)/125 = 0.56 ; var = 1/125 = 0.008 ; se=0.089443.
(() => {
    const cs = [
        { study: 'S1', t1: 'A', t2: 'B', es: 0.4, se: 0.2 },
        { study: 'S2', t1: 'A', t2: 'B', es: 0.6, se: 0.1 }
    ];
    const ts = ['A', 'B'];
    const est = Engine.estimate(cs, ts, 0);
    const pooled = est.theta[0] - est.theta[1];
    ok('IV-pooled A-B == 0.56 (hand-derived)', near(pooled, 0.56, 1e-9), pooled.toFixed(6));
    const vAB = est.V[0][0] + est.V[1][1] - 2 * est.V[0][1];
    ok('var(pooled A-B) == 1/125 = 0.008 (hand-derived)', near(vAB, 0.008, 1e-9), vAB.toFixed(7));
})();

console.log('=== Q statistic / heterogeneity ===');
// Same two A-B studies. Common-effect Q = sum w_i (y_i - ybar)^2 with ybar=0.56.
// = 25*(0.4-0.56)^2 + 100*(0.6-0.56)^2 = 25*0.0256 + 100*0.0016 = 0.64 + 0.16 = 0.80.
(() => {
    const cs = [
        { study: 'S1', t1: 'A', t2: 'B', es: 0.4, se: 0.2 },
        { study: 'S2', t1: 'A', t2: 'B', es: 0.6, se: 0.1 }
    ];
    const ts = ['A', 'B'];
    const est = Engine.estimate(cs, ts, 0);
    ok('Q == 0.80 (hand-derived)', near(est.Q, 0.80, 1e-9), est.Q.toFixed(6));
    // df = m - n + 1 = 2 - 2 + 1 = 1. I2 = max(0,(Q-df)/Q) = 0 since Q<df.
    ok('I2(Q=0.8, df=1) == 0 (Q<df)', near(Engine.I2(0.80, 1), 0, 1e-9));
    ok('I2(Q=10, df=5) == 50', near(Engine.I2(10, 5), 50, 1e-9));
    ok('I2(Q=100, df=10) == 90', near(Engine.I2(100, 10), 90, 1e-9));
})();

console.log('=== P-scores (ranking) properties ===');
// In a consistent network of K treatments, sum of P-scores = K/2 and mean = 0.5 (Ruecker & Schwarzer 2015).
(() => {
    const cs = [
        { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.1 },
        { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 0.1 },
        { study: 'S3', t1: 'A', t2: 'C', es: 0.8, se: 0.1 }
    ];
    const ts = ['A', 'B', 'C'];
    const est = Engine.estimate(cs, ts, 0);
    const ps = Engine.pScores(est.theta, est.V, false);
    ok('all P-scores in [0,1]', ps.every(p => p >= 0 && p <= 1), ps.map(p => p.toFixed(3)).join(','));
    const sum = ps.reduce((a, b) => a + b, 0);
    ok('sum(P-scores) == K/2 == 1.5', near(sum, ts.length / 2, 1e-6), sum.toFixed(6));
    ok('mean(P-scores) == 0.5', near(sum / ts.length, 0.5, 1e-6));
    // With larger ES = better (default smallBetter=false), A has highest theta -> highest P-score.
    ok('A ranked best (highest P-score)', ps[0] > ps[1] && ps[0] > ps[2]);
    // smallBetter flips: P-score_small = 1 - P-score_large
    const psSmall = Engine.pScores(est.theta, est.V, true);
    ok('smallBetter flips: ps_small_i == 1 - ps_large_i',
        ps.every((p, i) => near(p + psSmall[i], 1, 1e-9)));
})();

console.log('=== P-score formula vs explicit Phi computation ===');
// P-score_i = mean_{j!=i} Phi((theta_i - theta_j)/se_ij). Verify against a manual recomputation.
(() => {
    const cs = [
        { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.1 },
        { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 0.1 },
        { study: 'S3', t1: 'A', t2: 'C', es: 0.8, se: 0.1 }
    ];
    const ts = ['A', 'B', 'C'];
    const est = Engine.estimate(cs, ts, 0);
    const ps = Engine.pScores(est.theta, est.V, false);
    // manual for A (index 0)
    const n = 3, V = est.V, th = est.theta;
    let manual = 0;
    for (let j = 0; j < n; j++) {
        if (j !== 0) {
            const d = th[0] - th[j];
            const se = Math.sqrt(V[0][0] + V[j][j] - 2 * V[0][j]);
            manual += Engine.pnorm(d / se);
        }
    }
    manual /= (n - 1);
    ok('P-score(A) matches manual Phi mean', near(ps[0], manual, 1e-9), ps[0].toFixed(6));
})();

console.log('=== SUCRA ~ P-score agreement (stochastic) ===');
(() => {
    const cs = [
        { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.1 },
        { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 0.1 },
        { study: 'S3', t1: 'A', t2: 'C', es: 0.8, se: 0.1 }
    ];
    const ts = ['A', 'B', 'C'];
    const est = Engine.estimate(cs, ts, 0);
    const out = Engine.SUCRA(est.theta, est.V, false);
    ok('SUCRA values all in [0,1]', out.sucra.every(s => s >= 0 && s <= 1), out.sucra.map(s => s.toFixed(3)).join(','));
    // SUCRA should agree with P-scores within Monte-Carlo noise (2000 sims). Use loose tol.
    const agree = out.sucra.every((s, i) => Math.abs(s - out.pscores[i]) < 0.05);
    ok('SUCRA ~ P-scores within 0.05 (MC, Ruecker 2015)', agree,
        'sucra=[' + out.sucra.map(s => s.toFixed(3)).join(',') + '] pscore=[' + out.pscores.map(p => p.toFixed(3)).join(',') + ']');
})();

console.log('=== Connectivity guard ===');
(() => {
    ok('connected network detected',
        Engine.isConnected([{ t1: 'A', t2: 'B' }, { t1: 'B', t2: 'C' }], ['A', 'B', 'C']));
    ok('disconnected network detected',
        !Engine.isConnected([{ t1: 'A', t2: 'B' }, { t1: 'C', t2: 'D' }], ['A', 'B', 'C', 'D']));
})();

console.log('=== HKSJ widens (or equals) CI vs Wald ===');
(() => {
    const cs = [
        { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.1 },
        { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 0.1 },
        { study: 'S3', t1: 'A', t2: 'C', es: 0.8, se: 0.1 }
    ];
    const ts = ['A', 'B', 'C'];
    const est = Engine.estimate(cs, ts, 0);
    const df = cs.length - ts.length + 1; // = 1
    const wald = Engine.pairwise(est.theta, est.V, ts, 'MD', 0, df, false);
    const hksj = Engine.pairwise(est.theta, est.V, ts, 'MD', 0, df, true);
    const wAB = wald.find(p => p.t1 === 'A' && p.t2 === 'B');
    const hAB = hksj.find(p => p.t1 === 'A' && p.t2 === 'B');
    ok('HKSJ (t) CI width >= Wald (z) CI width', (hAB.hi - hAB.lo) >= (wAB.hi - wAB.lo) - 1e-9,
        `wald=${(wAB.hi - wAB.lo).toFixed(3)} hksj=${(hAB.hi - hAB.lo).toFixed(3)}`);
})();

console.log('=== Ratio-scale (OR) back-transform on log scale ===');
// For OR metric, effects are pooled on log scale and exp()-d. log(A-B)=0.5 -> OR=exp(0.5)=1.64872.
(() => {
    const cs = [
        { study: 'S1', t1: 'A', t2: 'B', es: 0.5, se: 0.1 },
        { study: 'S2', t1: 'B', t2: 'C', es: 0.3, se: 0.1 },
        { study: 'S3', t1: 'A', t2: 'C', es: 0.8, se: 0.1 }
    ];
    const ts = ['A', 'B', 'C'];
    const est = Engine.estimate(cs, ts, 0);
    const pw = Engine.pairwise(est.theta, est.V, ts, 'OR', 0, null, false);
    const ab = pw.find(p => p.t1 === 'A' && p.t2 === 'B');
    ok('OR(A vs B) == exp(0.5) == 1.64872 (log-scale pooling)', near(ab.d, Math.exp(0.5), 1e-6), ab.d.toFixed(6));
    ok('OR point is exp of log effect, CI bounds also exp-d', ab.dL < ab.d && ab.d < ab.dH);
})();

console.log('=== Built-in Validation module (10 internal tests) ===');
(() => {
    const results = Validation.runAll();
    results.forEach(r => ok('Validation: ' + r.name, r.pass, r.details));
})();

console.log('\n' + passed + ' passed, ' + failed + ' failed');
process.exit(failed === 0 ? 0 : 1);
