"""Exhaustive ideal-lattice search for an extremal 7-modular J (dim 24, min 10) with Z[(1+sqrt-7)/2]-structure.
K = Q(zeta_m), [K:Q] = 24, sqrt(-7) in K (7 | m).  L = (I, Tr(lam x conj(y))) with lam in K^+ totally positive,
lam * Diff * I * conj(I) = (sqrt(-7))  =>  L^# = sqrt(-7)^(-1) L, det 7^12.  Enumerate I over ideal classes
(representatives) and lam over all totally positive unit multiples (mod squares); exact qfminim for norm <= 8."""
import sys, json, itertools, time
import cypari2
pari = cypari2.Pari(); pari.allocatemem(4 * 10**9)
m = int(sys.argv[1]); out = sys.argv[2]
pari(f"K = bnfinit(polcyclo({m},t), 1); nf = K.nf; n = poldegree(nf.pol);")
pari("cj(x) = nfalgtobasis(nf, subst(lift(nfbasistoalg(nf, x)), t, t^(%d-1)) %% nf.pol);" % m)
pari("s7 = nfroots(nf, x^2 + 7)[1]; s7 = nfalgtobasis(nf, s7);")
pari("emb(x) = nfeltembed(nf, x);")
print("degree", pari("n"), "class group", pari("K.cyc"), flush=True)
# torsion, fundamental units of K, and fundamental units of the real subfield K+ = Q(t + 1/t)
pari("w = K.tu[1]; tu = nfalgtobasis(nf, K.tu[2]); U = [nfalgtobasis(nf, u) | u <- K.fu];")
pari("Kp = bnfinit(polredabs(subst(minpoly(Mod(t + t^(%d-1), polcyclo(%d, t))), t, y)), 1);" % (m, m))
# embed K+ into K: find a root of Kp.pol in K
pari("rt = nfroots(nf, Kp.pol); emb_p(a) = nfalgtobasis(nf, subst(lift(a), x, rt[1]));")
pari("Ur = concat([nfalgtobasis(nf, -1)], [emb_p(u) | u <- Kp.fu]);")
# root-of-unity exponent a with cj(u) = tu^a * u
pari("rexp(v) = my(q = nfeltdiv(nf, cj(v), v)); for(k = 0, w - 1, if(nfeltpow(nf, tu, k) == q, return(k))); -1;")
print("unit rank", pari("#U"), "torsion order", pari("w"), "real units (incl -1)", pari("#Ur"), flush=True)
# ideal class representatives
pari("reps = [idealhnf(nf,1)]; cyc = K.cyc; gens = K.gen;")
if int(pari("#cyc")) > 0:
    pari("reps = []; forvec(e = vector(#cyc, i, [0, cyc[i]-1]), reps = concat(reps, [idealfactorback(nf, gens, e)]));")
print("ideal class reps:", pari("#reps"), flush=True)
results = []
t0 = time.time()
nreps = int(pari("#reps"))
for ri in range(1, nreps + 1):
    pari(f"Id = idealhnf(nf, reps[{ri}]); J0 = idealmul(nf, idealmul(nf, nf.diff, Id), idealhnf(nf, Mat([cj(v) | v <- Vec(Id)])));")
    # target principal ideal: (sqrt-7) / (Diff I conj(I))
    pari("Q = idealdiv(nf, idealhnf(nf, s7), J0); pr = bnfisprincipal(K, Q);")
    if any(int(x) for x in pari("pr[1]")):
        print(f"class rep {ri}: target not principal, skip", flush=True); continue
    pari("lam0 = pr[2];")
    # make lam real: cj(lam0) = tu^j lam0; multiply by eps = tu^k * prod U_i^e_i
    pari("j = rexp(lam0); A = [rexp(u) | u <- U];")
    pari("lamr = 0; forvec(e = vector(#U, i, [0, 1]), s = (j + sum(i = 1, #U, e[i] * A[i])) % w; "
         "if(s % 2 == 0, eps = nfeltpow(nf, tu, s / 2); for(i = 1, #U, if(e[i], eps = nfeltmul(nf, eps, U[i]))); "
         "l = nfeltmul(nf, lam0, eps); if(cj(l) == l, lamr = l; break)));")
    if pari("lamr") == 0:
        print(f"class rep {ri}: no unit makes lam real; skipping", flush=True); continue
    nr = int(pari("#Ur")) - 1
    print(f"class rep {ri}: real lam found; real units {nr} + sign", flush=True)
    seen = set()
    count = 0
    for signs in itertools.product([0, 1], repeat=nr + 1):
        pari("l = lamr;")
        if signs[0]: pari("l = -l;")
        for i, e in enumerate(signs[1:]):
            if e: pari(f"l = nfeltmul(nf, l, Ur[{i+2}]);")
        if not bool(pari("vecmin(real(emb(l))) > 0")):
            continue
        count += 1
        pari("Bs = Vec(Id); G = matrix(n, n, i, j, nfelttrace(nf, nfeltmul(nf, l, nfeltmul(nf, Bs[i], cj(Bs[j])))));")
        if not bool(pari("denominator(G) == 1")):
            continue
        even = bool(pari("vecsum(vector(n, i, G[i,i] % 2)) == 0"))
        det = int(pari("matdet(G)"))
        mi = pari("qfminim(G, 8, 0, 2)")
        short = int(mi[0]) // 2
        minimum = int(pari("my(V = qfminim(G, 8, , 0)[3]); vecmin(vector(#V, k, V[,k]~ * G * V[,k]))")) if short else None
        key = (det, short, minimum)
        results.append(dict(class_rep=ri, signs=list(signs), even=even, det=det, short_pairs=short, minimum=minimum))
        print(f"  lam #{count}: even {even}, det = 7^{round(__import__('math').log(det, 7), 3) if det > 0 else det}, "
              f"pairs of norm <= 8: {short}, min {minimum}  [{time.time()-t0:.0f}s]", flush=True)
        if even and det == 7 ** 12 and short == 0:
            print("FOUND extremal 7-modular ideal lattice!", flush=True)
json.dump(results, open(out, "w"), indent=1)
print("done", len(results), "lattices tested")
