/* csr.h -- matrizes esparsas em formato CSR (compressed sparse row).
 *
 * Módulo fornecido com o TP2 de Álgebra Linear Computacional (UFJF, 2026.3).
 *
 * Convenções:
 *   - índices a partir de 0;
 *   - rowptr tem m+1 posições; colind e val têm nnz posições;
 *   - os coeficientes da linha i estão em rowptr[i] .. rowptr[i+1]-1,
 *     com colunas em ordem crescente e sem repetição.
 * Os índices são int e os valores double, os mesmos tipos de PetscInt e
 * PetscScalar na configuração padrão do PETSc, de modo que os três vetores
 * podem ser passados diretamente a MatCreateSeqAIJWithArrays.
 */
#ifndef CSR_H
#define CSR_H

#ifdef __cplusplus
extern "C"
{
#endif

    typedef struct
    {
        int     m, n;   /* número de linhas e de colunas */
        int     nnz;    /* número de coeficientes armazenados */
        int    *rowptr; /* m+1 posições */
        int    *colind; /* nnz posições */
        double *val;    /* nnz posições */
    } CSR;

    /* Aloca uma matriz m x n com espaço para nnz coeficientes.
 * rowptr sai zerado; colind e val, sem inicialização. Para montar uma matriz,
 * preencha os três vetores. Devolve NULL se faltar memória. */
    CSR *csr_alloc(int m, int n, int nnz);

    /* Libera a matriz (aceita NULL). */
    void csr_free(CSR *A);

    /* Lê um arquivo Matrix Market (.mtx) no formato coordinate, com valores
 * real, integer ou pattern (pattern: todos os valores iguais a 1) e simetria
 * general, symmetric ou skew-symmetric. Nas matrizes symmetric e
 * skew-symmetric, que guardam só o triângulo inferior, os coeficientes fora
 * da diagonal são espelhados. Coeficientes repetidos são somados.
 * Devolve NULL e escreve a causa em stderr em caso de erro. */
    CSR *csr_read_mtx(const char *arquivo);

    /* y = A x. x tem n posições e y, m posições (x e y não podem coincidir). */
    void csr_matvec(const CSR *A, const double *x, double *y);

    /* Confere a estrutura: rowptr[0] = 0, rowptr não decrescente,
 * rowptr[m] = nnz e, em cada linha, colunas dentro de [0, n) e em ordem
 * estritamente crescente. Devolve 0 se estiver tudo certo; senão, escreve o
 * primeiro problema em stderr e devolve 1. Útil para conferir a montagem. */
    int csr_check(const CSR *A);

#ifdef __cplusplus
}
#endif

#endif /* CSR_H */
