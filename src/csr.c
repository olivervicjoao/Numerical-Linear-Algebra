/* csr.c -- implementação de csr.h. Ver as convenções no cabeçalho. */
#include "../include/csr.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

CSR *csr_alloc(int m, int n, int nnz)
{
    CSR *A = malloc(sizeof *A);
    if (!A)
        return NULL;
    A->m      = m;
    A->n      = n;
    A->nnz    = nnz;
    A->rowptr = calloc((size_t)m + 1, sizeof *A->rowptr);
    A->colind = malloc((size_t)(nnz > 0 ? nnz : 1) * sizeof *A->colind);
    A->val    = malloc((size_t)(nnz > 0 ? nnz : 1) * sizeof *A->val);
    if (!A->rowptr || !A->colind || !A->val)
    {
        csr_free(A);
        return NULL;
    }
    return A;
}

void csr_free(CSR *A)
{
    if (!A)
        return;
    free(A->rowptr);
    free(A->colind);
    free(A->val);
    free(A);
}

void csr_matvec(const CSR *A, const double *x, double *y)
{
    for (int i = 0; i < A->m; i++)
    {
        double s = 0.0;
        for (int k = A->rowptr[i]; k < A->rowptr[i + 1]; k++)
            s += A->val[k] * x[A->colind[k]];
        y[i] = s;
    }
}

int csr_check(const CSR *A)
{
    if (A->rowptr[0] != 0)
    {
        fprintf(stderr, "csr_check: rowptr[0] = %d (deveria ser 0)\n", A->rowptr[0]);
        return 1;
    }
    for (int i = 0; i < A->m; i++)
    {
        if (A->rowptr[i + 1] < A->rowptr[i])
        {
            fprintf(stderr, "csr_check: rowptr decresce na linha %d\n", i);
            return 1;
        }
        for (int k = A->rowptr[i]; k < A->rowptr[i + 1]; k++)
        {
            int j = A->colind[k];
            if (j < 0 || j >= A->n)
            {
                fprintf(stderr, "csr_check: coluna %d fora de [0,%d) na linha %d\n", j, A->n, i);
                return 1;
            }
            if (k > A->rowptr[i] && j <= A->colind[k - 1])
            {
                fprintf(stderr, "csr_check: colunas fora de ordem ou repetidas na linha %d\n", i);
                return 1;
            }
        }
    }
    if (A->rowptr[A->m] != A->nnz)
    {
        fprintf(stderr, "csr_check: rowptr[m] = %d, mas nnz = %d\n", A->rowptr[A->m], A->nnz);
        return 1;
    }
    return 0;
}

/* ------------------------------------------------------------------ leitura */

static void minusculas(char *s)
{
    for (; *s; s++)
        *s = (char)tolower((unsigned char)*s);
}

typedef struct
{
    int    col;
    double val;
} Par;

static int compara_par(const void *a, const void *b)
{
    int ca = ((const Par *)a)->col, cb = ((const Par *)b)->col;
    return (ca > cb) - (ca < cb);
}

CSR *csr_read_mtx(const char *arquivo)
{
    FILE *f = fopen(arquivo, "r");
    if (!f)
    {
        fprintf(stderr, "csr_read_mtx: não foi possível abrir %s\n", arquivo);
        return NULL;
    }

    char linha[1024], banner[64], obj[64], fmt[64], tipo[64], sim[64];
    if (!fgets(linha, sizeof linha, f) || sscanf(linha, "%63s %63s %63s %63s %63s", banner, obj, fmt, tipo, sim) != 5 || strcmp(banner, "%%MatrixMarket") != 0)
    {
        fprintf(stderr, "csr_read_mtx: %s não tem o cabeçalho %%%%MatrixMarket\n", arquivo);
        fclose(f);
        return NULL;
    }
    minusculas(obj);
    minusculas(fmt);
    minusculas(tipo);
    minusculas(sim);
    int pattern   = strcmp(tipo, "pattern") == 0;
    int simetrica = strcmp(sim, "symmetric") == 0;
    int antissim  = strcmp(sim, "skew-symmetric") == 0;
    if (strcmp(obj, "matrix") != 0 || strcmp(fmt, "coordinate") != 0 || !(pattern || !strcmp(tipo, "real") || !strcmp(tipo, "integer")) ||
        !(simetrica || antissim || !strcmp(sim, "general")))
    {
        fprintf(stderr, "csr_read_mtx: formato não suportado em %s: %s %s %s %s\n", arquivo, obj, fmt, tipo, sim);
        fclose(f);
        return NULL;
    }

    /* pula comentários e lê as dimensões */
    long m, n, nz;
    do
    {
        if (!fgets(linha, sizeof linha, f))
        {
            fprintf(stderr, "csr_read_mtx: %s terminou antes das dimensões\n", arquivo);
            fclose(f);
            return NULL;
        }
    } while (linha[0] == '%' || linha[strspn(linha, " \t\r\n")] == '\0');
    if (sscanf(linha, "%ld %ld %ld", &m, &n, &nz) != 3 || m <= 0 || n <= 0 || nz < 0)
    {
        fprintf(stderr, "csr_read_mtx: linha de dimensões inválida em %s\n", arquivo);
        fclose(f);
        return NULL;
    }

    /* triplas, já com o espelhamento */
    long    cap = (simetrica || antissim) ? 2 * nz : nz;
    int    *I   = malloc((size_t)(cap > 0 ? cap : 1) * sizeof *I);
    int    *J   = malloc((size_t)(cap > 0 ? cap : 1) * sizeof *J);
    double *V   = malloc((size_t)(cap > 0 ? cap : 1) * sizeof *V);
    if (!I || !J || !V)
    {
        fprintf(stderr, "csr_read_mtx: memória insuficiente\n");
        free(I);
        free(J);
        free(V);
        fclose(f);
        return NULL;
    }
    long t = 0;
    for (long k = 0; k < nz; k++)
    {
        long   i, j;
        double v = 1.0;
        int    lidos;
        lidos = pattern ? fscanf(f, "%ld %ld", &i, &j) : fscanf(f, "%ld %ld %lf", &i, &j, &v);
        if (lidos != (pattern ? 2 : 3) || i < 1 || i > m || j < 1 || j > n)
        {
            fprintf(stderr, "csr_read_mtx: entrada %ld inválida em %s\n", k + 1, arquivo);
            free(I);
            free(J);
            free(V);
            fclose(f);
            return NULL;
        }
        I[t] = (int)(i - 1);
        J[t] = (int)(j - 1);
        V[t] = v;
        t++;
        if ((simetrica || antissim) && i != j)
        {
            I[t] = (int)(j - 1);
            J[t] = (int)(i - 1);
            V[t] = antissim ? -v : v;
            t++;
        }
    }
    fclose(f);

    /* triplas -> CSR: conta por linha, distribui, ordena e soma repetidos */
    int *cont = calloc((size_t)m + 1, sizeof *cont);
    Par *P    = malloc((size_t)(t > 0 ? t : 1) * sizeof *P);
    if (!cont || !P)
    {
        fprintf(stderr, "csr_read_mtx: memória insuficiente\n");
        free(I);
        free(J);
        free(V);
        free(cont);
        free(P);
        return NULL;
    }
    for (long k = 0; k < t; k++)
        cont[I[k] + 1]++;
    for (long i = 0; i < m; i++)
        cont[i + 1] += cont[i];
    int *pos = malloc((size_t)m * sizeof *pos);
    memcpy(pos, cont, (size_t)m * sizeof *pos);
    for (long k = 0; k < t; k++)
    {
        P[pos[I[k]]].col = J[k];
        P[pos[I[k]]].val = V[k];
        pos[I[k]]++;
    }
    free(pos);
    free(I);
    free(J);
    free(V);

    long unicos = 0;
    for (long i = 0; i < m; i++)
    {
        qsort(P + cont[i], (size_t)(cont[i + 1] - cont[i]), sizeof *P, compara_par);
        for (int k = cont[i]; k < cont[i + 1]; k++)
            if (k == cont[i] || P[k].col != P[k - 1].col)
                unicos++;
    }

    CSR *A = csr_alloc((int)m, (int)n, (int)unicos);
    if (!A)
    {
        fprintf(stderr, "csr_read_mtx: memória insuficiente\n");
        free(cont);
        free(P);
        return NULL;
    }
    int q = 0;
    for (long i = 0; i < m; i++)
    {
        A->rowptr[i] = q;
        for (int k = cont[i]; k < cont[i + 1]; k++)
        {
            if (k > cont[i] && P[k].col == P[k - 1].col)
            {
                A->val[q - 1] += P[k].val; /* repetido: soma */
            }
            else
            {
                A->colind[q] = P[k].col;
                A->val[q]    = P[k].val;
                q++;
            }
        }
    }
    A->rowptr[m] = q;
    free(cont);
    free(P);
    return A;
}
