#include "pcgc.h"

/* produto matriz-vetor aproveitando banda k-diagonal */
static void matvec_band(const real_t *A, const real_t *v, real_t *out, int n, int k) {
    int d = (k - 1) / 2;
    for (int i = 0; i < n; ++i) {
        int j0 = i - d; if (j0 < 0) j0 = 0;
        int j1 = i + d; if (j1 > n-1) j1 = n-1;
        real_t s = 0.0;
        int row = i * n;
        for (int j = j0; j <= j1; ++j) s += A[row + j] * v[j];
        out[i] = s;
    }
}

/* dot product */
static real_t dot(const real_t *a, const real_t *b, int n) {
    real_t s = 0.0;
    for (int i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

/* norma L2 */
static real_t norm2(const real_t *a, int n) {
    return sqrt(dot(a,a,n));
}

/* norma infinito entre x e y */
static real_t infnorm_diff(const real_t *x, const real_t *y, int n) {
    real_t m = 0.0;
    for (int i = 0; i < n; ++i) {
        real_t v = fabs(x[i] - y[i]);
        if (v > m) m = v;
    }
    return m;
}

/* aplica pré-condicionador M^{-1} r -> z
   - se M == NULL => z = r
   - se M != NULL => assumimos M é o vetor diagonal (Jacobi): z[i] = r[i] / M[i]
*/
static int apply_precond(const real_t *r, real_t *z, int n, const real_t *M) {
    if (M == NULL) {
        for (int i = 0; i < n; ++i) z[i] = r[i];
        return 0;
    } else {
        for (int i = 0; i < n; ++i) {
            if (M[i] == 0.0) return -1;
            z[i] = r[i] / M[i];
        }
        return 0;
    }
}

int conjugateGradient(real_t *A, real_t *b, real_t *x, int n, int k,
                      double eps, int maxit,
                      real_t *residuo_out, real_t *M, real_t *norma_inf_out)
{
    /* aloca vetores */
    real_t *r = (real_t *)malloc(n * sizeof(real_t));
    real_t *z = (real_t *)malloc(n * sizeof(real_t));
    real_t *p = (real_t *)malloc(n * sizeof(real_t));
    real_t *Ap = (real_t *)malloc(n * sizeof(real_t));
    real_t *x_prev = (real_t *)malloc(n * sizeof(real_t));
    if (!r || !z || !p || !Ap || !x_prev) {
        fprintf(stderr, "ERROR: allocation failed in conjugateGradient\n");
        free(r); free(z); free(p); free(Ap); free(x_prev);
        return -1;
    }

    /* inicializa x, x_prev */
    for (int i = 0; i < n; ++i) { /* assumimos x já alocado; definir como zero se necessário */
        /* se x não estiver inicializado, o caller deve ter setado calloc */
        x_prev[i] = 0.0;
    }

    /* r = b - A*x ; se x inicial = 0 => r = b */
    matvec_band(A, x, Ap, n, k); /* Ap temporário */
    for (int i = 0; i < n; ++i) r[i] = b[i] - Ap[i];

    /* z = M^{-1} r */
    if (apply_precond(r, z, n, M) != 0) {
        fprintf(stderr, "ERROR: division by zero in preconditioner\n");
        free(r); free(z); free(p); free(Ap); free(x_prev);
        return -1;
    }

    /* p = z */
    for (int i = 0; i < n; ++i) p[i] = z[i];

    real_t rz_old = dot(r, z, n);
    if (rz_old == 0.0) {
        /* solution is already zero */
        *residuo_out = norm2(r, n);
        *norma_inf_out = 0.0;
        free(r); free(z); free(p); free(Ap); free(x_prev);
        return 0;
    }

    int iter = 0;
    real_t norma_inf = 0.0;
    for (iter = 0; iter < maxit; ++iter) {
        /* Ap = A * p */
        matvec_band(A, p, Ap, n, k);

        real_t denom = dot(p, Ap, n);
        if (denom == 0.0) {
            fprintf(stderr, "ERROR: breakdown p^T A p == 0\n");
            free(r); free(z); free(p); free(Ap); free(x_prev);
            return -1;
        }
        real_t alpha = rz_old / denom;

        /* x_new = x + alpha * p */
        for (int i = 0; i < n; ++i) x[i] += alpha * p[i];

        /* r = r - alpha * Ap */
        for (int i = 0; i < n; ++i) r[i] -= alpha * Ap[i];

        /* z = M^{-1} r */
        if (apply_precond(r, z, n, M) != 0) {
            fprintf(stderr, "ERROR: division by zero in preconditioner (iter %d)\n", iter+1);
            free(r); free(z); free(p); free(Ap); free(x_prev);
            return -1;
        }

        real_t rz_new = dot(r, z, n);
        real_t beta = rz_new / rz_old;

        /* p = z + beta * p */
        for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];

        /* norma infinita entre x e x_prev */
        norma_inf = infnorm_diff(x, x_prev, n);

        /* copia x para x_prev */
        for (int i = 0; i < n; ++i) x_prev[i] = x[i];

        rz_old = rz_new;

        if (norma_inf < eps) {
            /* convergiu */
            iter = iter + 1; /* número de iterações efetivamente realizadas */
            break;
        }
    }

    /* calcula residuo L2 final */
    /* note: r já é residuo b - A x (mantido) */
    real_t resid_l2 = norm2(r, n);
    *residuo_out = resid_l2;
    *norma_inf_out = norma_inf;

    free(r); free(z); free(p); free(Ap); free(x_prev);
    return iter;
}
