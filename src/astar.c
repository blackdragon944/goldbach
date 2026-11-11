/*
 * astar.c -- fast computation of a*(n) for the Sierpinski-Goldbach conjecture.
 *
 * For even n and odd prime a, n is "covered by a" if there are odd primes p, q
 * with p + q = n, q > a, and p a binary submask of n - a  ((p & (n-a)) == p,
 * equivalently C(n-a, p) odd).  a*(n) = least such a.
 *
 * Primality: an odd-only bit sieve up to N, filled by a segmented sieve of
 * Eratosthenes (segments run in parallel).  Memory: N/16 bytes.
 *
 * Modes
 *   astar scan  N [threads] [dumpN dumpfile]
 *       a*(n) for all even 8 <= n <= N: exceptions, records, histograms per
 *       decade; optionally dumps a*(n) (uint16, index n/2) for n <= dumpN.
 *   astar count A D outfile
 *       G_A(n) = #{p : p odd prime, p subset n-A, p < n-A, n-p prime}
 *       for every even n <= D, written as uint32 (index n/2).
 *   astar carries N [threads]
 *       a = 3 only, allow up to c carries in p + (q-3); prints, for c = 0..3,
 *       the number of even 10 <= n <= N with no pair of <= c carries.
 *       (min carries are capped at 4 = "more than 3")
 *   astar ext N [threads]
 *       a*(n), Chen-type a*_C(n), binary vertex rule a_bin(n), explicit
 *       construction a_k(n) for all even 10 <= n <= N (see run_ext).
 *   astar tries X
 *       order-of-vertices statistics for the n <= X where a = 3 fails.
 *   astar mult full N0 N1 [early] | mult sample N0 N1 count seed
 *       covering multiplicity M(p, n) of Goldbach pairs (see run_mult).
 *
 * Build: see README.md  (zig cc / gcc / clang, -O2, pthreads).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <time.h>
#include <math.h>

typedef uint64_t u64;
typedef uint32_t u32;

/* ------------------------------------------------------------------ sieve */

static uint8_t *bits;    /* bit i  <->  odd number 2i+1 is prime */
static u64 LIMIT;        /* sieve covers [0, LIMIT] */
static u32 *base;        /* odd primes <= sqrt(LIMIT) */
static u32 nbase;
static int NTHREADS = 4;

static inline int isprime(u64 x)
{
    if (!(x & 1)) return x == 2;
    return (bits[x >> 4] >> ((x >> 1) & 7)) & 1;
}

#define SEG_BYTES (1u << 17) /* 128 KiB = 1M odd numbers per segment */

static u64 nbytes, nsegs;
static volatile long long seg_next;

static void sieve_segment(u64 s)
{
    u64 b0 = s * SEG_BYTES, b1 = b0 + SEG_BYTES;
    if (b1 > nbytes) b1 = nbytes;
    u64 i0 = b0 * 8, i1 = b1 * 8; /* odd indices [i0, i1) <-> numbers 2i+1 */
    uint8_t *B = bits;
    for (u32 k = 0; k < nbase; k++) {
        u64 pr = base[k];
        u64 sq = pr * pr;
        if (sq > 2 * i1 + 1) break;
        u64 lo = 2 * i0 + 1; /* first number in segment */
        u64 st;
        if (sq >= lo) st = sq;
        else {
            st = (lo + pr - 1) / pr * pr;
            if (!(st & 1)) st += pr;
        }
        for (u64 j = st >> 1; j < i1; j += pr) B[j >> 3] &= (uint8_t)~(1u << (j & 7));
    }
}

static void *sieve_worker(void *arg)
{
    (void)arg;
    for (;;) {
        long long s = __atomic_fetch_add(&seg_next, 1, __ATOMIC_RELAXED);
        if ((u64)s >= nsegs) break;
        sieve_segment((u64)s);
    }
    return NULL;
}

static void build_sieve(u64 N)
{
    LIMIT = N;
    nbytes = (N >> 4) + 1;
    bits = (uint8_t *)malloc(nbytes);
    if (!bits) { fprintf(stderr, "out of memory (%llu bytes)\n", (unsigned long long)nbytes); exit(1); }
    memset(bits, 0xFF, nbytes);
    bits[0] &= ~1u; /* 1 is not prime */

    u64 r = 1;
    while ((r + 1) * (r + 1) <= N) r++;
    uint8_t *small = (uint8_t *)calloc(r + 2, 1);
    base = (u32 *)malloc(sizeof(u32) * (r / 2 + 2));
    nbase = 0;
    for (u64 i = 3; i <= r; i += 2) {
        if (small[i]) continue;
        base[nbase++] = (u32)i;
        for (u64 j = i * i; j <= r; j += 2 * i) small[j] = 1;
    }
    free(small);

    nsegs = (nbytes + SEG_BYTES - 1) / SEG_BYTES;
    seg_next = 0;
    pthread_t th[64];
    for (int t = 0; t < NTHREADS; t++) pthread_create(&th[t], NULL, sieve_worker, NULL);
    for (int t = 0; t < NTHREADS; t++) pthread_join(th[t], NULL);
    /* clear bits beyond LIMIT */
    for (u64 x = N + 1; x < nbytes * 16; x++)
        if (x & 1) bits[x >> 4] &= (uint8_t)~(1u << ((x >> 1) & 7));
}

/* ------------------------------------------------------- the core search */

#define MAXA 8192 /* number of odd primes a tried; a* beyond this is reported */
static u32 oddp[MAXA];

static void init_oddp(void)
{
    int k = 0;
    for (u32 x = 3; k < MAXA; x += 2) {
        int pr = 1;
        for (u32 d = 3; d * d <= x; d += 2)
            if (x % d == 0) { pr = 0; break; }
        if (pr) oddp[k++] = x;
    }
}

/* does a cover n?  enumerate odd submasks p = 2s+1 of m = n-a, p < m */
static inline int covers(u64 n, u64 a)
{
    u64 m = n - a;          /* odd */
    u64 h = m >> 1;         /* p = 2s+1 with s a proper submask of h */
    if (h == 0) return 0;
    u64 s = (h - 1) & h;
    for (;;) {
        u64 p = 2 * s + 1;
        if (p >= 3 && isprime(p) && isprime(n - p)) return 1;
        if (s == 0) return 0;
        s = (s - 1) & h;
    }
}

/* index k of a*(n) = oddp[k], or -1 if none (with a + 5 <= n) */
static int astar_idx(u64 n)
{
    for (int k = 0; k < MAXA; k++) {
        u64 a = oddp[k];
        if (a + 5 > n) return -1;
        if (covers(n, a)) return k;
    }
    return -2; /* ran out of a's: report */
}

/* --------------------------------------------------------------- scan mode */

#define NDEC 12
#define CHUNK (1ull << 21) /* n values per chunk (even and odd counted) */
#define MAXREC 256

typedef struct { u64 n; int k; } rec_t;
typedef struct { rec_t r[MAXREC]; int nr; } chunkrec_t;

static u64 SCAN_N, DUMP_N;
static uint16_t *dump;
static chunkrec_t *crec;
static u64 nchunks;
static volatile long long chunk_next;
static pthread_mutex_t exc_lock = PTHREAD_MUTEX_INITIALIZER;
static u64 exc[1024];
static int nexc;

typedef struct { int tid; u64 h[NDEC][256]; u64 hbig[NDEC]; int kmax; } tstate_t;

static int decade(u64 n)
{
    int d = 0;
    u64 p = 1;
    while (p < n) { p *= 10; d++; }
    return d; /* smallest d with n <= 10^d */
}

static void *scan_worker(void *arg)
{
    tstate_t *ts = (tstate_t *)arg;
    for (;;) {
        long long c = __atomic_fetch_add(&chunk_next, 1, __ATOMIC_RELAXED);
        if ((u64)c >= nchunks) break;
        u64 lo = (u64)c * CHUNK, hi = lo + CHUNK;
        if (lo < 8) lo = 8;
        if (hi > SCAN_N + 1) hi = SCAN_N + 1;
        chunkrec_t *cr = &crec[c];
        cr->nr = 0;
        int runmax = -1;
        int d = decade(lo);
        u64 dlim = 1;
        for (int i = 0; i < d; i++) dlim *= 10;
        for (u64 n = lo; n < hi; n += 2) {
            while (n > dlim) { d++; dlim *= 10; }
            int k = astar_idx(n);
            if (k < 0) {
                pthread_mutex_lock(&exc_lock);
                if (nexc < 1024) exc[nexc++] = n;
                pthread_mutex_unlock(&exc_lock);
                if (k == -2) fprintf(stderr, "a*(%llu) > %u (MAXA exceeded)\n", (unsigned long long)n, oddp[MAXA - 1]);
                if (n <= DUMP_N) dump[n >> 1] = 0;
                continue;
            }
            if (k < 256) ts->h[d][k]++; else ts->hbig[d]++;
            if (k > ts->kmax) ts->kmax = k;
            if (k > runmax) {
                runmax = k;
                if (cr->nr < MAXREC) { cr->r[cr->nr].n = n; cr->r[cr->nr].k = k; cr->nr++; }
            }
            if (n <= DUMP_N) dump[n >> 1] = (uint16_t)oddp[k];
        }
    }
    return NULL;
}

static double now(void)
{
    struct timespec t;
    timespec_get(&t, TIME_UTC);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static int cmp_u64(const void *a, const void *b)
{
    u64 x = *(const u64 *)a, y = *(const u64 *)b;
    return x < y ? -1 : x > y;
}

static int run_scan(u64 N, u64 dumpN, const char *dumpfile)
{
    double t0 = now();
    build_sieve(N);
    double t1 = now();
    printf("# sieve up to %llu: %.2f s\n", (unsigned long long)N, t1 - t0);
    SCAN_N = N;
    DUMP_N = dumpN;
    if (dumpN) {
        dump = (uint16_t *)calloc(dumpN / 2 + 1, sizeof(uint16_t));
        if (!dump) { fprintf(stderr, "out of memory for dump\n"); return 1; }
    }
    nchunks = (N + 1 + CHUNK - 1) / CHUNK;
    crec = (chunkrec_t *)calloc(nchunks, sizeof(chunkrec_t));
    tstate_t *ts = (tstate_t *)calloc(NTHREADS, sizeof(tstate_t));
    pthread_t th[64];
    chunk_next = 0;
    for (int t = 0; t < NTHREADS; t++) { ts[t].tid = t; ts[t].kmax = -1; pthread_create(&th[t], NULL, scan_worker, &ts[t]); }
    for (int t = 0; t < NTHREADS; t++) pthread_join(th[t], NULL);
    double t2 = now();
    printf("# scan of even n in [8, %llu]: %.2f s, %d threads\n", (unsigned long long)N, t2 - t1, NTHREADS);

    /* exceptions */
    qsort(exc, nexc, sizeof(u64), cmp_u64);
    printf("exceptions %d:", nexc);
    for (int i = 0; i < nexc; i++) printf(" %llu", (unsigned long long)exc[i]);
    printf("\n");

    /* records */
    int gmax = -1;
    printf("records (n a*):\n");
    for (u64 c = 0; c < nchunks; c++)
        for (int i = 0; i < crec[c].nr; i++)
            if (crec[c].r[i].k > gmax) {
                gmax = crec[c].r[i].k;
                printf("record %llu %u\n", (unsigned long long)crec[c].r[i].n, oddp[gmax]);
            }

    /* cumulative histograms: n <= 10^d */
    int kmax = -1;
    for (int t = 0; t < NTHREADS; t++) if (ts[t].kmax > kmax) kmax = ts[t].kmax;
    u64 cum[256] = {0}, cumbig = 0;
    for (int d = 0; d < NDEC; d++) {
        u64 any = 0;
        for (int t = 0; t < NTHREADS; t++) {
            for (int k = 0; k < 256; k++) { cum[k] += ts[t].h[d][k]; any += ts[t].h[d][k]; }
            cumbig += ts[t].hbig[d];
            any += ts[t].hbig[d];
        }
        u64 p10 = 1;
        for (int i = 0; i < d; i++) p10 *= 10;
        if (d < 2) continue;
        u64 upto = p10 < N ? p10 : N;
        printf("hist n<=%llu:", (unsigned long long)upto);
        for (int k = 0; k < 256 && k <= kmax; k++) if (cum[k]) printf(" %u:%llu", oddp[k], (unsigned long long)cum[k]);
        if (cumbig) printf(" big:%llu", (unsigned long long)cumbig);
        printf("\n");
        if (p10 >= N) break;
        (void)any;
    }
    if (dumpN) {
        FILE *f = fopen(dumpfile, "wb");
        if (!f) { perror(dumpfile); return 1; }
        fwrite(dump, sizeof(uint16_t), dumpN / 2 + 1, f);
        fclose(f);
        printf("# dumped a*(n) for n <= %llu to %s\n", (unsigned long long)dumpN, dumpfile);
    }
    return 0;
}

/* -------------------------------------------------------------- count mode */

static int run_count(u64 A, u64 D, const char *outfile)
{
    build_sieve(D);
    u32 *G = (u32 *)calloc(D / 2 + 1, sizeof(u32));
    for (u64 n = A + 5; n <= D; n++) {
        if (n & 1) continue;
        u64 m = n - A, h = m >> 1;
        u32 cnt = 0;
        if (h) {
            u64 s = (h - 1) & h;
            for (;;) {
                u64 p = 2 * s + 1;
                if (p >= 3 && isprime(p) && isprime(n - p)) cnt++;
                if (s == 0) break;
                s = (s - 1) & h;
            }
        }
        G[n >> 1] = cnt;
    }
    FILE *f = fopen(outfile, "wb");
    if (!f) { perror(outfile); return 1; }
    fwrite(G, sizeof(u32), D / 2 + 1, f);
    fclose(f);
    printf("# G_%llu(n) for even n <= %llu written to %s\n", (unsigned long long)A, (unsigned long long)D, outfile);
    return 0;
}

/* ------------------------------------------------------------ carries mode */

/* is there an odd prime p <= m, with q = n - p prime, q > 3, and at most cmax
   carries in p + (m - p), m = n - 3?  Digits of p are chosen from the bottom;
   the carry out of each position is determined (it equals the borrow of m - p). */
static int dfs_carry(u64 n, u64 m, int bit, int L, u64 p, int cin, int used, int cmax)
{
    if (bit == L) {
        if (cin) return 0; /* p > m */
        if (p >= 3 && p < m && isprime(p) && isprime(n - p)) return 1;
        return 0;
    }
    int mi = (int)((m >> bit) & 1);
    for (int pi = 1; pi >= 0; pi--) {
        if (bit == 0 && pi == 0) continue; /* p odd */
        int s = mi - pi - cin;
        int cout = s < 0;
        if (used + cout > cmax) continue;
        if (dfs_carry(n, m, bit + 1, L, p | ((u64)pi << bit), cout, used + cout, cmax)) return 1;
    }
    return 0;
}

static u64 CAR_N;
static volatile long long car_next;
static u64 car_cnt[64][6];

static void *carries_worker(void *arg)
{
    int t = (int)(intptr_t)arg;
    for (;;) {
        long long c = __atomic_fetch_add(&car_next, 1, __ATOMIC_RELAXED);
        u64 lo = (u64)c * 4096, hi = lo + 4096;
        if (lo > CAR_N) break;
        if (lo < 10) lo = 10;
        if (hi > CAR_N + 1) hi = CAR_N + 1;
        for (u64 n = lo; n < hi; n += 2) {
            u64 m = n - 3;
            int L = 0;
            while ((m >> L) != 0) L++;
            int c0 = 4;
            for (int cm = 0; cm <= 3; cm++)
                if (dfs_carry(n, m, 0, L, 0, 0, 0, cm)) { c0 = cm; break; }
            car_cnt[t][c0]++;
        }
    }
    return NULL;
}

static int run_carries(u64 N)
{
    build_sieve(N);
    CAR_N = N;
    car_next = 0;
    pthread_t th[64];
    for (int t = 0; t < NTHREADS; t++) pthread_create(&th[t], NULL, carries_worker, (void *)(intptr_t)t);
    for (int t = 0; t < NTHREADS; t++) pthread_join(th[t], NULL);
    u64 tot[6] = {0};
    for (int t = 0; t < NTHREADS; t++) for (int c = 0; c < 6; c++) tot[c] += car_cnt[t][c];
    /* exceptions for c = number of n whose min carries > c */
    for (int c = 0; c <= 3; c++) {
        u64 e = 0;
        for (int j = c + 1; j <= 4; j++) e += tot[j];
        printf("carries<=%d exceptions_upto_%llu %llu\n", c, (unsigned long long)N, (unsigned long long)e);
    }
    return 0;
}

/* ================================================================ ext mode
 *
 * One pass over all even 10 <= n <= N computing, for every n,
 *   a*(n)            least covering odd prime (as in scan mode);
 *   a*_C(n)          least odd prime a covering n in the Chen sense: some odd
 *                    submask p of n-a, 3 <= p < n-a, with one of p, q = n-p
 *                    prime and the other P_2 (Omega <= 2, i.e. prime or a
 *                    product of two primes);
 *   a_bin(n)         among odd primes a < BIN_B (a <= n-5) the one with the
 *                    largest popcount(n-a), ties -> smaller a; we check that it
 *                    covers n.  Also the same rule with the bound
 *                    B = c (ln n)^2, c in CB[] (a = 3 if B <= 3);
 *   a_k(n), k=2..9   d = n mod 2^k, a_k = least odd prime > d; the best k is the
 *                    one with the largest popcount(n - a_k) (a_k < BIN_B).
 */

#define BIN_B 1000
#define NCB 4
static const double CB[NCB] = {0.5, 1.0, 2.0, 3.0};
#define NAK_LO 2
#define NAK_HI 9
static int nbinp;                 /* number of odd primes < BIN_B */
static u32 nextoddprime[1 << NAK_HI];

static inline int popc(u64 x) { return __builtin_popcountll(x); }

/* Omega(x) <= 2 for odd x >= 3 */
static inline int isP2(u64 x)
{
    if (isprime(x)) return 1;
    for (int k = 0;; k++) {
        u64 r = oddp[k];
        if (r * r * r > x) return 1;  /* no factor <= cbrt(x): at most two prime factors */
        if (x % r == 0) return isprime(x / r);
    }
}

static int chen_covers(u64 n, u64 a)
{
    u64 m = n - a, h = m >> 1;
    if (h == 0) return 0;
    u64 s = (h - 1) & h;
    for (;;) {
        u64 p = 2 * s + 1;
        if (p >= 3) {
            u64 q = n - p;
            int pp = isprime(p), qp = isprime(q);
            if (pp && qp) return 1;
            if (pp && isP2(q)) return 1;
            if (qp && isP2(p)) return 1;
        }
        if (s == 0) return 0;
        s = (s - 1) & h;
    }
}

#define NDEF 48
typedef struct {
    u64 cnt[NDEC], hard3[NDEC], chenhard3[NDEC], binfail[NDEC];
    u64 hC[NDEC][256], hCbig[NDEC];
    u64 akcov[NDEC], akeq[NDEC], akdef[NDEC][NDEF];
    u64 cfail[NDEC][NCB]; u64 cfl[NCB][8]; int ncfl[NCB];
    u64 habin[256];
    u64 lemma_bad;                 /* sanity check of the lemma, must stay 0 */
    u32 abin_max; u64 abin_max_n;
    u64 tries_dummy;
} xstate_t;

typedef struct { u64 n; int k; } xrec_t;
typedef struct { xrec_t r[MAXREC]; int nr; } xchunk_t;
static xchunk_t *xrec;             /* records of a*_C per chunk */
static u64 chenlist[4096]; static int nchenlist;
static u64 binexc[1024]; static int nbinexc;
static u64 XN;

/* covers(n, a) using what a* already tells us */
static inline int covers_k(u64 n, int k, int kstar)
{
    if (k < kstar) return 0;
    if (k == kstar) return 1;
    return covers(n, oddp[k]);
}

static void *ext_worker(void *arg)
{
    xstate_t *st = (xstate_t *)arg;
    int cp[1024];
    for (;;) {
        long long c = __atomic_fetch_add(&chunk_next, 1, __ATOMIC_RELAXED);
        if ((u64)c >= nchunks) break;
        u64 lo = (u64)c * CHUNK, hi = lo + CHUNK;
        if (lo < 10) lo = 10;
        if (hi > XN + 1) hi = XN + 1;
        xchunk_t *cr = &xrec[c];
        cr->nr = 0;
        int runmax = -1;
        int d = decade(lo);
        u64 dlim = 1;
        for (int i = 0; i < d; i++) dlim *= 10;
        for (u64 n = lo; n < hi; n += 2) {
            while (n > dlim) { d++; dlim *= 10; }
            st->cnt[d]++;
            int ks = astar_idx(n);
            if (ks < 0) { fprintf(stderr, "no a* for n=%llu\n", (unsigned long long)n); continue; }
            /* --- Chen */
            int kc = ks;
            for (int k = 0; k < ks; k++)
                if (chen_covers(n, oddp[k])) { kc = k; break; }
            if (ks > 0) st->hard3[d]++;
            if (kc > 0) {
                st->chenhard3[d]++;
                if (n <= 1000000) {
                    pthread_mutex_lock(&exc_lock);
                    if (nchenlist < 4096) chenlist[nchenlist++] = n;
                    pthread_mutex_unlock(&exc_lock);
                }
            }
            if (kc < 256) st->hC[d][kc]++; else st->hCbig[d]++;
            if (kc > runmax) {
                runmax = kc;
                if (cr->nr < MAXREC) { cr->r[cr->nr].n = n; cr->r[cr->nr].k = kc; cr->nr++; }
            }
            /* --- binary rule: change points of the running argmax over a < BMAX */
            double ln = log((double)n), l2 = ln * ln;
            double bmax = BIN_B;
            if (CB[NCB - 1] * l2 > bmax) bmax = CB[NCB - 1] * l2;
            int ncp = 0, best = -1, bestbin = -1;
            for (int k = 0;; k++) {
                u64 a = oddp[k];
                if (a + 5 > n || (double)a >= bmax) break;
                int s = popc(n - a);
                if (s > best) { best = s; cp[ncp++] = k; }
                if (a < BIN_B) bestbin = ncp - 1;
            }
            int kbin = cp[bestbin];
            best = popc(n - oddp[kbin]);
            int okbin = covers_k(n, kbin, ks);
            if (!okbin) {
                st->binfail[d]++;
                pthread_mutex_lock(&exc_lock);
                if (nbinexc < 1024) binexc[nbinexc++] = n;
                pthread_mutex_unlock(&exc_lock);
            }
            st->habin[kbin < 256 ? kbin : 255]++;
            u32 abin = oddp[kbin];
            if (abin > st->abin_max) { st->abin_max = abin; st->abin_max_n = n; }
            /* rule with B = c (ln n)^2: last change point below B (a = 3 if B <= 3) */
            for (int ci = 0; ci < NCB; ci++) {
                double B = CB[ci] * l2;
                int j = 0;
                while (j + 1 < ncp && (double)oddp[cp[j + 1]] < B) j++;
                int ok = (cp[j] == kbin) ? okbin : covers_k(n, cp[j], ks);
                if (!ok) {
                    st->cfail[d][ci]++;
                    if (st->ncfl[ci] < 8) st->cfl[ci][st->ncfl[ci]++] = n;
                }
            }
            /* --- explicit construction a_k */
            int bk = -1, bs = -1;
            u64 ba = 0;
            for (int k = NAK_LO; k <= NAK_HI; k++) {
                u64 dd = n & ((1ull << k) - 1);
                u64 a = nextoddprime[dd];
                if (a >= BIN_B || a + 5 > n) continue;
                u64 e = a - dd;
                if (e <= (1ull << k) && n >= (1ull << k)) {
                    u64 low = (n - a) & ((1ull << k) - 1);
                    if (low != (1ull << k) - e) st->lemma_bad++;
                }
                int s = popc(n - a);
                if (s > bs || (s == bs && a < ba)) { bs = s; bk = k; ba = a; }
            }
            if (bk >= 0) {
                int kk = 0;
                while (oddp[kk] != ba) kk++;
                if (covers_k(n, kk, ks)) st->akcov[d]++;
                if (kk == kbin) st->akeq[d]++;
                int def = best - bs;
                st->akdef[d][def < NDEF ? def : NDEF - 1]++;
            }
        }
    }
    return NULL;
}

static int run_ext(u64 N)
{
    double t0 = now();
    build_sieve(N);
    printf("# sieve up to %llu: %.2f s\n", (unsigned long long)N, now() - t0);
    nbinp = 0;
    while (oddp[nbinp] < BIN_B) nbinp++;
    for (u32 dd = 0; dd < (1u << NAK_HI); dd++) {
        int k = 0;
        while (oddp[k] <= dd) k++;
        nextoddprime[dd] = oddp[k];
    }
    XN = N;
    nchunks = (N + 1 + CHUNK - 1) / CHUNK;
    xrec = (xchunk_t *)calloc(nchunks, sizeof(xchunk_t));
    xstate_t *st = (xstate_t *)calloc(NTHREADS, sizeof(xstate_t));
    pthread_t th[64];
    chunk_next = 0;
    double t1 = now();
    for (int t = 0; t < NTHREADS; t++) pthread_create(&th[t], NULL, ext_worker, &st[t]);
    for (int t = 0; t < NTHREADS; t++) pthread_join(th[t], NULL);
    printf("# ext scan of even n in [10, %llu]: %.2f s, %d threads, B = %d\n", (unsigned long long)N, now() - t1,
           NTHREADS, BIN_B);

    u64 lemma_bad = 0; u32 abmax = 0; u64 abmaxn = 0;
    for (int t = 0; t < NTHREADS; t++) {
        lemma_bad += st[t].lemma_bad;
        if (st[t].abin_max > abmax || (st[t].abin_max == abmax && st[t].abin_max_n < abmaxn)) {
            abmax = st[t].abin_max; abmaxn = st[t].abin_max_n;
        }
    }
    printf("lemma_violations %llu\n", (unsigned long long)lemma_bad);
    printf("abin_max %u at %llu\n", abmax, (unsigned long long)abmaxn);
    qsort(binexc, nbinexc, sizeof(u64), cmp_u64);
    printf("abin_exceptions_list");
    for (int i = 0; i < nbinexc; i++) printf(" %llu", (unsigned long long)binexc[i]);
    printf("\n");
    qsort(chenlist, nchenlist, sizeof(u64), cmp_u64);
    printf("chen_hard3_upto_1e6 %d:", nchenlist);
    for (int i = 0; i < nchenlist; i++) printf(" %llu", (unsigned long long)chenlist[i]);
    printf("\n");
    int gmax = -1;
    printf("chen records (n a*_C):\n");
    for (u64 c = 0; c < nchunks; c++)
        for (int i = 0; i < xrec[c].nr; i++)
            if (xrec[c].r[i].k > gmax) {
                gmax = xrec[c].r[i].k;
                printf("crecord %llu %u\n", (unsigned long long)xrec[c].r[i].n, oddp[gmax]);
            }

    /* cumulative statistics, n <= 10^d */
    u64 cnt = 0, h3 = 0, ch3 = 0, bf = 0, akc = 0, ake = 0, hC[256] = {0}, hCb = 0, def[NDEF] = {0}, cf[NCB] = {0};
    u64 hab[256] = {0};
    for (int t = 0; t < NTHREADS; t++) for (int k = 0; k < 256; k++) hab[k] += st[t].habin[k];
    printf("histAbin:");
    for (int k = 0; k < 256; k++) if (hab[k]) printf(" %u:%llu", oddp[k], (unsigned long long)hab[k]);
    printf("\n");
    for (int ci = 0; ci < NCB; ci++) {
        printf("cfail_first c=%.1f:", CB[ci]);
        for (int t = 0; t < NTHREADS; t++) for (int i = 0; i < st[t].ncfl[ci]; i++) printf(" %llu", (unsigned long long)st[t].cfl[ci][i]);
        printf("\n");
    }
    for (int d = 0; d < NDEC; d++) {
        for (int t = 0; t < NTHREADS; t++) {
            xstate_t *s = &st[t];
            cnt += s->cnt[d]; h3 += s->hard3[d]; ch3 += s->chenhard3[d]; bf += s->binfail[d];
            akc += s->akcov[d]; ake += s->akeq[d]; hCb += s->hCbig[d];
            for (int k = 0; k < 256; k++) hC[k] += s->hC[d][k];
            for (int k = 0; k < NDEF; k++) def[k] += s->akdef[d][k];
            for (int ci = 0; ci < NCB; ci++) cf[ci] += s->cfail[d][ci];
        }
        u64 p10 = 1;
        for (int i = 0; i < d; i++) p10 *= 10;
        if (d < 2) continue;
        u64 upto = p10 < N ? p10 : N;
        printf("upto %llu count %llu hard3 %llu chenhard3 %llu binfail %llu akcover %llu akeqbin %llu cfail",
               (unsigned long long)upto, (unsigned long long)cnt, (unsigned long long)h3, (unsigned long long)ch3,
               (unsigned long long)bf, (unsigned long long)akc, (unsigned long long)ake);
        for (int ci = 0; ci < NCB; ci++) printf(" %llu", (unsigned long long)cf[ci]);
        printf("\n");
        printf("histC n<=%llu:", (unsigned long long)upto);
        for (int k = 0; k < 256; k++) if (hC[k]) printf(" %u:%llu", oddp[k], (unsigned long long)hC[k]);
        if (hCb) printf(" big:%llu", (unsigned long long)hCb);
        printf("\n");
        printf("akdeficit n<=%llu:", (unsigned long long)upto);
        for (int k = 0; k < NDEF; k++) if (def[k]) printf(" %d:%llu", k, (unsigned long long)def[k]);
        printf("\n");
        fflush(stdout);
        if (p10 >= N) break;
    }
    return 0;
}

/* ============================================================== tries mode
 * For the "hard" n (a = 3 does not cover), 10 <= n <= X, the number of vertices
 * tried until the first covering one, for three orders of the odd primes a < B:
 *   asc  : 3, 5, 7, ...
 *   div  : odd prime divisors of n first (ascending), then the rest ascending
 *   bin  : decreasing popcount(n-a), ties -> smaller a
 * Also: how often a*(n) divides n, against the expectation sum 1/a*(n) and
 * against a shuffled pairing (a*(n') | n for a cyclic shift n' of the hard n). */
static int run_tries(u64 X, const char *outfile)
{
    build_sieve(X);
    int nb = 0;
    while (oddp[nb] < BIN_B) nb++;
    u64 hard = 0, tasc = 0, tdiv = 0, tbin = 0, divcnt = 0, nofound = 0;
    double expct = 0;
    u64 cap = 1 << 22, nh = 0;
    u64 *hn = malloc(cap * sizeof(u64));
    u32 *ha = malloc(cap * sizeof(u32));
    int order[1024], ok[1024];
    for (u64 n = 10; n <= X; n += 2) {
        if (covers(n, 3)) continue;
        hard++;
        int na = 0;
        for (int k = 0; k < nb && oddp[k] + 5 <= n; k++) na++;
        for (int k = 0; k < na; k++) ok[k] = -1;   /* lazily evaluated */
#define OK(k) (ok[k] < 0 ? (ok[k] = covers(n, oddp[k])) : ok[k])
        int ks = -1;
        for (int k = 0; k < na; k++) if (OK(k)) { ks = k; break; }
        if (ks < 0) { nofound++; continue; }
        tasc += ks + 1;
        u32 as = oddp[ks];
        if (n % as == 0) divcnt++;
        expct += 1.0 / as;
        if (nh < cap) { hn[nh] = n; ha[nh] = as; nh++; }
        /* divisors first */
        int no = 0;
        for (int k = 0; k < na; k++) if (n % oddp[k] == 0) order[no++] = k;
        for (int k = 0; k < na; k++) if (n % oddp[k] != 0) order[no++] = k;
        for (int i = 0; i < no; i++) if (OK(order[i])) { tdiv += i + 1; break; }
        /* binary order: selection by decreasing popcount */
        int used[1024] = {0};
        for (int i = 0; i < na; i++) {
            int bk = -1, bs = -1;
            for (int k = 0; k < na; k++)
                if (!used[k]) { int s = popc(n - oddp[k]); if (s > bs) { bs = s; bk = k; } }
            used[bk] = 1;
            if (OK(bk)) { tbin += i + 1; break; }
        }
#undef OK
    }
    u64 shuf = 0;
    for (u64 i = 0; i < nh; i++) if (hn[i] % ha[(i + nh / 2) % nh] == 0) shuf++;
    if (outfile) {   /* (n, a*(n)) of the hard n, for a permutation test in analysis2.py */
        FILE *f = fopen(outfile, "wb");
        if (!f) { perror(outfile); return 1; }
        for (u64 i = 0; i < nh; i++) { fwrite(&hn[i], 8, 1, f); u64 a = ha[i]; fwrite(&a, 8, 1, f); }
        fclose(f);
    }
    printf("tries X %llu B %d hard %llu nofound %llu\n", (unsigned long long)X, BIN_B, (unsigned long long)hard,
           (unsigned long long)nofound);
    printf("mean_tries asc %.4f div %.4f bin %.4f\n", (double)tasc / hard, (double)tdiv / hard, (double)tbin / hard);
    printf("astar_divides_n %llu expected_sum_1_over_a %.1f shuffled %llu\n", (unsigned long long)divcnt, expct,
           (unsigned long long)shuf);
    return 0;
}

/* =============================================================== mult mode
 * M(p, n) = #{odd primes a < q = n - p : p subset n - a}, for Goldbach pairs
 * n = p + q (odd primes).  Equivalently the number of odd m with p subset m,
 * p < m <= n - 3 and n - m prime: we run over the supersets m of p in
 * increasing order.
 *   mult full N0 N1 [early]   all even n in [N0, N1]; with early = 1 only
 *                             whether M >= 1 (stops at the first a).
 *   mult sample N0 N1 cnt seed  cnt random even n in [N0, N1], all pairs,
 *                             M >= 1 only, plus exact M for <= 400 random
 *                             pairs per n.
 * Output: per-n lines "n <n> pairs covered covered_q>100 pairs_q>100 sumM lowE uncovered_lowE uncovered_other(q>100, E>=2)
 * uncovered_q ..." (lowE: pairs whose model expectation pi(q-1)/2^(s(p)-1) of M is < 2) and
 * aggregated histograms (see analysis). */
static u64 Mcount(u64 n, u64 p, int early)
{
    u64 L = 64 - __builtin_clzll(n);
    u64 full = (L >= 64) ? ~0ull : ((1ull << L) - 1);
    u64 comp = full & ~p;
    u64 cnt = 0;
    u64 x = 0;
    for (;;) {
        x = (x - comp) & comp;      /* next superset increment, increasing */
        if (x == 0) break;
        u64 m = p | x;
        if (m > n - 3) break;
        if (isprime(n - m)) { cnt++; if (early) return 1; }
    }
    return cnt;
}

/* number of odd primes <= x, via a prefix table over 512-byte blocks */
static u32 *piblk;
static void build_pi(void)
{
    u64 nb = nbytes / 512 + 2;
    piblk = (u32 *)calloc(nb, sizeof(u32));
    u32 c = 0;
    for (u64 b = 0; b < nb; b++) {
        piblk[b] = c;
        for (u64 j = b * 512; j < (b + 1) * 512 && j < nbytes; j++) c += __builtin_popcount(bits[j]);
    }
}
static u64 piodd(u64 x)
{
    if (x < 3) return 0;
    u64 i = (x - 1) / 2;            /* bit indices 0..i  <->  odd numbers 1..x */
    u64 byte = i >> 3, b = byte / 512;
    u64 c = piblk[b];
    for (u64 j = b * 512; j < byte; j++) c += __builtin_popcount(bits[j]);
    c += __builtin_popcount(bits[byte] & ((2u << (i & 7)) - 1));
    return c;
}

/* exact-M statistics: by popcount(p), log2 bucket of M, and (bit length of n, popcount(p)) */
static double sumM[64], sumM2[64][64], sumR2[64][64];
static u64 numM[64], histM[64], numM2[64][64];
static void accumM(u64 n, u64 p, u64 M)
{
    int s = popc(p), L = 64 - __builtin_clzll(n);
    sumM[s] += (double)M; numM[s]++;
    histM[M ? 64 - __builtin_clzll(M) : 0]++;
    sumM2[L][s] += (double)M; numM2[L][s]++;
    /* model: each odd prime a < q independently has p subset n-a with prob 2^(1-s) */
    double E = (double)(piodd(n - p - 1)) / (double)(1ull << (s - 1));
    if (E > 0) sumR2[L][s] += (double)M / E;
}

static u64 rng_state;
static u64 rnd(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

static int run_mult(int argc, char **argv)
{
    int sample = !strcmp(argv[2], "sample");
    u64 N0 = strtoull(argv[3], NULL, 10), N1 = strtoull(argv[4], NULL, 10);
    int early = 0;
    u64 cnt = 0, seed = 1;
    if (sample) { cnt = strtoull(argv[5], NULL, 10); seed = strtoull(argv[6], NULL, 10); }
    else if (argc >= 6) early = atoi(argv[5]);
    build_sieve(N1);
    rng_state = seed * 0x9E3779B97F4A7C15ull + 1;
    u64 nn = sample ? cnt : (N1 - N0) / 2 + 1;
    build_pi();
    for (u64 i = 0; i < nn; i++) {
        u64 n = sample ? (N0 + 2 * (rnd() % ((N1 - N0) / 2 + 1))) : N0 + 2 * i;
        n &= ~1ull;
        if (n < 10) continue;
        u64 pairs = 0, cov = 0, pb = 0, cb = 0, zero_q[16]; int nz = 0;
        double sM = 0;
        u64 lowE = 0, uncLowE = 0, uncOther = 0;   /* pairs (and uncovered pairs) with model E = pi(q-1)/2^(s(p)-1) < 2 */
        for (u64 p = 3; p <= n - 3; p += 2) {
            if (!isprime(p) || !isprime(n - p)) continue;
            u64 q = n - p;
            pairs++;
            u64 M = Mcount(n, p, early || sample);
            if (M) cov++; else if (nz < 16) zero_q[nz++] = q;
            {
                double E = (double)piodd(q - 1) / (double)(1ull << (popc(p) - 1));
                if (E < 2) { lowE++; if (!M) uncLowE++; }
                else if (!M && q > 100) uncOther++;
            }
            if (q > 100) { pb++; if (M) cb++; }
            if (!early && !sample) { accumM(n, p, M); sM += (double)M; }
        }
        if (sample) {
            /* exact M on up to 400 random pairs */
            int got = 0, tries = 0;
            while (got < 400 && tries < 400000) {
                tries++;
                u64 p = 3 + 2 * (rnd() % ((n - 6) / 2 + 1));
                if (!isprime(p) || !isprime(n - p)) continue;
                got++;
                u64 M = Mcount(n, p, 0);
                accumM(n, p, M);
            }
        }
        printf("n %llu %llu %llu %llu %llu %.0f %llu %llu %llu uncovered_q", (unsigned long long)n, (unsigned long long)pairs,
               (unsigned long long)cov, (unsigned long long)cb, (unsigned long long)pb, sM, (unsigned long long)lowE,
               (unsigned long long)uncLowE, (unsigned long long)uncOther);
        for (int j = 0; j < nz; j++) printf(" %llu", (unsigned long long)zero_q[j]);
        printf("\n");
    }
    printf("meanM_by_popcount_p");
    for (int s = 0; s < 64; s++) if (numM[s]) printf(" %d:%.3f:%llu", s, sumM[s] / numM[s], (unsigned long long)numM[s]);
    printf("\nhistM_log2bucket");
    for (int b = 0; b < 64; b++) if (histM[b]) printf(" %d:%llu", b, (unsigned long long)histM[b]);
    printf("\n");
    for (int L = 0; L < 64; L++) {
        int any = 0;
        for (int t = 0; t < 64; t++) any |= numM2[L][t] != 0;
        if (!any) continue;
        printf("byLs L=%d", L);
        for (int t = 0; t < 64; t++)
            if (numM2[L][t]) printf(" %d:%.4f:%.4f:%llu", t, sumM2[L][t] / numM2[L][t], sumR2[L][t] / numM2[L][t],
                                    (unsigned long long)numM2[L][t]);
        printf("\n");
    }
    return 0;
}

/* --------------------------------------------------------------------- main */

int main(int argc, char **argv)
{
    init_oddp();
    if (argc >= 3 && !strcmp(argv[1], "scan")) {
        u64 N = strtoull(argv[2], NULL, 10);
        if (argc >= 4) NTHREADS = atoi(argv[3]);
        u64 dN = 0;
        const char *df = NULL;
        if (argc >= 6) { dN = strtoull(argv[4], NULL, 10); df = argv[5]; }
        if (dN > N) dN = N;
        return run_scan(N, dN, df);
    }
    if (argc >= 5 && !strcmp(argv[1], "count"))
        return run_count(strtoull(argv[2], NULL, 10), strtoull(argv[3], NULL, 10), argv[4]);
    if (argc >= 3 && !strcmp(argv[1], "carries")) {
        if (argc >= 4) NTHREADS = atoi(argv[3]);
        return run_carries(strtoull(argv[2], NULL, 10));
    }
    if (argc >= 3 && !strcmp(argv[1], "ext")) {
        if (argc >= 4) NTHREADS = atoi(argv[3]);
        return run_ext(strtoull(argv[2], NULL, 10));
    }
    if (argc >= 3 && !strcmp(argv[1], "tries"))
        return run_tries(strtoull(argv[2], NULL, 10), argc >= 4 ? argv[3] : NULL);
    if (argc >= 5 && !strcmp(argv[1], "mult"))
        return run_mult(argc, argv);
    fprintf(stderr,
            "usage:\n  astar scan N [threads] [dumpN dumpfile]\n"
            "  astar count A D outfile\n  astar carries N [threads]\n");
    return 2;
}
