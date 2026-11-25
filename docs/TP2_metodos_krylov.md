## Programa de Pós-graduação emModelagem Computacional Universidade Federal de Juiz de Fora

## Álgebra Linear Computacional

Trabalho prático 2: métodos de Krylov pré-condicionados e cálculo de autovalores

Prof. José J. Camata

Modalidade: individual. Entrega do código e do relatório: 07/12/2026

Apresentação em aula: 07/12/2026, 20 minutos por aluno.

2026.3

## Objetivo

Implementar três métodos de Krylov para matrizes esparsas em formato CSR — o método dos gra- dientes conjugados (GC), o GMRES com reinício e o BiCGSTAB —, com e sem o pré-condicionador diagonal (Jacobi). Aplicá-los a problemas modelo e a matrizes de um repositório público, analisar número de iterações, custo e convergência, e comparar os resultados com as implementações da biblioteca PETSc (Partes 1 a 5). Implementar o método da potência, a potência inversa com deslocamento, que usa os resolvedores de Krylov das partes anteriores, e o método de Lanczos para calcular autovalores extremos das mes- mas matrizes. Usar esses autovalores para explicar o comportamento dos resolvedores e comparar os resultados com a LAPACK e, opcionalmente, com a biblioteca SLEPc (Partes 6 a 8).

## Ambiente

- Linguagem: C ou C++. Quem preferir Python deve seguir o Apêndice A, que dá as ferramentas equivalentes.

- Matrizes em formato CSR (compressed sparse row), pelo módulo csr.c/csr.h fornecido com este enunciado: vetores rowptr (m+1 posições), colind e val (nnz posições), índices a partir de 0 e colunas em ordem crescente dentro de cada linha. O módulo também lê arquivos no formato Matrix Market (csr_read_mtx), calcula o produto y Ax (csr_matvec) e confere a estrutura de uma matriz montada (csr_check). As funções estão documentadas em csr.h.

• Biblioteca PETSc, versão 3.19 ou mais recente, apenas para a comparação da Parte 5. Em Debi-

an/Ubuntu: pacote petsc-dev; ou compilação a partir de https://petsc.org. [URL 🔗](https://petsc.org)

• Biblioteca SLEPc, da mesma versão do PETSc, apenas para a comparação opcional da Parte 8. Em

Debian/Ubuntu: pacote slepc-dev; ou compilação a partir de https://slepc.upv.es. [URL 🔗](https://slepc.upv.es)

•

pacote liblapacke-dev).

Exemplo de compilação com o PETSc instalado pelo pacote do Ubuntu:

PETSC_DIR=/usr/lib/petscdir/petsc3.19/x86_64-linux-gnu-real mpicc -O2 tp.c csr.c -o tp -I\$PETSC_DIR/include \

-L\$PETSC_DIR/lib -lpetsc_real -lm -Wl,-rpath,\$PETSC_DIR/lib

Para usar a LAPACKE, acrescente -llapacke -llapack antes de -lm. Quem fizer a comparação op- cional como SLEPc compila com:

LAPACK, para as rotinas dstev e dsyevd (pacote liblapack-dev; em C, a interface LAPACKE,


```
SLEPC_DIR=/usr/lib/slepcdir/slepc3.19/x86_64-linux-gnu-real
mpicc -O2 tp.c csr.c -o tp -I$SLEPC_DIR/include -I$PETSC_DIR/include \
-L$SLEPC_DIR/lib -lslepc_real -L$PETSC_DIR/lib -lpetsc_real \
-llapacke -llapack -lm -Wl,-rpath,$SLEPC_DIR/lib -Wl,-rpath,$PETSC_DIR/lib
```

Todas as execuções são seriais (um processo). Use diana de pelo menos 5 execuções.

clock_gettime

ou

std::chrono e reporte a me-

As operações com vetores, o pré-condicionador, os três métodos de Krylov e os métodos de autova- lores devem ser implementados por você. A estrutura CSR, a leitura de arquivos e o produto matriz- vetor vêm do módulo fornecido. Da LAPACK, só podem ser usadas as rotinas indicadas no texto. O PETSc e o SLEPc só podem ser usados como referência de comparação.

## Problemas de teste

## Problemas modelo

São dois problemas das aulas, no quadrado Ω= (0,1)2, com u 0 na fronteira e f 1, discretizados em uma malha N×N de pontos interiores, h 1/(N+ 1), m N2 e ordem natural das incógnitas

(k jN+i):

- 1. Poisson, matriz é SPD. f, com o estêncil de 5 pontos 4ui ¡ j −ui−1,j −ui+1,j j−1 j+1 ¢ /h2 fi j . A

- 2. Convecção-difusão, termos γ(ui+1,j −ui−1,j )/(2h) e γ(ui, j+1 u,= f, com γ 20. Ao estêncil de Poisson somam-se os j−1)/(2h). A matriz não é simétrica.

Com a ordem natural, cada bloco de N incógnitas corresponde a uma linha horizontal j da malha, e

as duas matrizes são tridiagonais por blocos, com blocos N×N. Com β=γh/2,

−c−I T

em que I é a identidade N×N. No problema 1, γ 0 e c+ c− 1, o que dá a matriz simétrica com T tridiag(−1,4,−1). No problema 2, β 10h e os coeficientes abaixo e acima da diagonal diferem. Cada linha tem no máximo 5 coeficientes não nulos, então nnz 5N2 −4N. A diagonal é constante, igual a 4/h2, nos dois problemas.

Monte as matrizes diretamente em CSR, sem passar pela forma densa: crie a matriz com a rotina csr_alloc, preencha os três vetores e confira o resultado com csr_check. Use N∈ {32,64,128,256},

b (1, . . . ,1)T (isto é, f 1) e x0 0.

## Matrizes do repositório

Escolha na SuiteSparse Matrix Collection, https://sparse.tamu.edu, pelo menos duas matrizes reais e quadradas, com 104 m 106: pelo menos uma simétrica positiva definida e pelo menos uma não simétrica. Justifique a escolha no relatório (origem do problema, tamanho, propriedades). Exemplos que atendem aos critérios, apenas como referência: [URL 🔗](https://sparse.tamu.edu)


| Matriz | Origem | m | nnz SPD |
| --- | --- | --- | --- |
| ACUSIM/Pres_Poisson | dinâmica dos fluidos |   | 14 822 715804 sim |
| GHS_psdef/apache1 | diferenças finitas 3D |   | 80 800 542184 sim |
| Schmid/thermal1 | elementos finitos, térmico 82 654 574458 sim |   |   |
| Averous/epb1 | trocador de calor |   | 14 734 95053 não |
| Wang/wang3 | dispositivo semicondutor 26 064 177168 não |   |   |

Nessas matrizes, tome a solução exata x (1, . . . ,1)T, defina b Ax e use x0 0. Não é garantido que todos os métodos convirjam em todas as matrizes: um resultado de não convergência, bem docu- mentado e explicado, também é um resultado.

Atenção: em matrizes marcadas como symmetric no cabeçalho do arquivo .mtx, só a parte triangu- lar inferior está armazenada. A rotina csr_read_mtx já espelha os coeficientes fora da diagonal, de modo que o nnz da matriz lida é maior que o número de entradas do arquivo.

## Parte 1 –Montagem e caracterização das matrizes

- 1. Leia as matrizes do repositório com csr_read_mtx.

- 2. Implemente a montagem dos dois problemas modelo em CSR.

- 3. Para cada matriz usada, reporte emuma tabelam, nnz, nnz/m, se é simétrica e se a diagonal tem elementos nulos. Inclua uma figura do padrão de esparsidade de cada matriz do repositório.

Valide a montagem com csr_check e, para N pequeno (por exemplo N= 4), confira algumas linhas da matriz contra o estêncil.

## Parte 2 – Implementações

Implemente os três métodos, com interface comum: entrada A (CSR), b, x0, tolerância ε, númeromá- ximo de iterações kmax e o pré-condicionador; saída xk, número de iterações, histórico de e número de produtos matriz-vetor.

- 1. GC pré-condicionado, como no algoritmo da aula.

- 2. GMRES(ℓ) com rotações de Givens, Gram-Schmidt modificado no processo de Arnoldi, pré- condicionado.

- 3. BiCGSTAB pré-condicionado, a partir do algoritmo abaixo.

- 4. Pré-condicionador de Jacobi, M= D diag(a11, . . . ,amm), aplicado como zi ri/aii . Todos os métodos devem rodar também comM= I (sem pré-condicionador).

A forma de incorporar o pré-condicionador ao GMRES e ao BiCGSTAB fica a seu critério. Descreva e justifique a escolha no relatório.

Critério de parada, igual para todos:

próprio método, sem produtos matriz-vetor adicionais, ou k kmax 20000. Se, na forma que você implementou, o método fornecer naturalmente outra norma, explique no relatório o critério usado. Ao final, calcule sempre o resíduo verdadeiro Axk∥2 /∥b∥2. Conte como iteração do BiCGSTAB cada passagem completa pelo laço abaixo (dois produtos por A).

ε com ε 10−8, usando o resíduo disponível no


## BiCGSTAB

O BiCGSTAB (van der Vorst, 1992) serve para matrizes não simétricas, tem recorrência curta (custo e memória fixos por iteração) e não minimiza nenhuma norma do resíduo, o que pode tornar a con- vergência irregular. O algoritmo abaixo está na forma sem pré-condicionador; cabe a você incluir M. Dados A, b, x0, ε e kmax:

```
r0 Ax0, rˆ r0, ρ−1 v p
para k = 0, 1, ..., kmax−1
ρk ˆrTrk (se ρk falha)
β= (ρk/ρk−1) (α/ω)
p rk
v Ap
α=ρk/( ˆrTv)
s rk
se ε: xk+1 xk parar
t As
ω= (tTs)/(tTt)
xk+1 xk
rk+1 s−ωt
se ε: parar
fim
```

Valide as três implementações em problemas pequenos (por exemplo N= 8), comparando o número de iterações com o PETSc (Parte 5) antes dos experimentos.

## Parte 3 – Experimentos com os problemas modelo

Rode cada método aplicável em cada problema e cada N, sem pré-condicionador e com Jacobi. O

GC só se aplica ao problema 1; no GMRES(ℓ), use ℓ=30.

- 1. Reporte o número de iterações, o número de produtos matriz-vetor, o tempo total, o tempo por iteração, o número de vetores de tamanho m que o método guarda e o resíduo verdadeiro ao final.

- 2. Faça gráficos do histórico de com todos os métodos. (escala logarítmica) para N 128, um por problema,

- 3. Faça um gráfico log-log do número de iterações em função de N e estime a taxa de crescimento de cada método.

- 4. Estudo do reinício: para N= 64 e N= 128, rode o GMRES(ℓ) com Jacobi para ℓ {10,30,100} nos dois problemas e compare iterações, tempo e memória.

## Parte 4 –Matrizes do repositório

Repita o item 1 da Parte 3 para as matrizes escolhidas, com os métodos aplicáveis (o GC só nas SPD), sem pré-condicionador e com Jacobi. Reporte também o erro ˆx−x∥∞, já que a solução exata é conhecida, e o histórico do resíduo de cada caso.


## Parte 5 – Comparação com o PETSc

Resolva os mesmos sistemas das Partes 3 e 4 com o KSP do PETSc, usando a mesma matriz CSR das suas implementações. A rotina MatCreateSeqAIJWithArrays cria uma matriz do PETSc a partir dos três vetores do CSR, sem copiá-los:

```
Mat A; Vec x, b; KSP ksp;
MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, m, m, rowptr, colind, val, &A);
MatCreateVecs(A, &x, &b);
VecSet(b, 1.0); VecSet(x, 0.0);
KSPCreate(PETSC_COMM_SELF, &ksp);
KSPSetOperators(ksp, A, A);
KSPSetFromOptions(ksp); /* metodo e opcoes pela linha de comando */
KSPSolve(ksp, b, x);
KSPGetIterationNumber(ksp, &its);
```

Os vetores rowptr e colind devem ser do tipo PetscInt (int, na configuração padrão) e val do tipo

PetscScalar (double).

Para que o PETSc resolva o mesmo problema, com o mesmo critério de parada que as suas imple- mentações, use as opções:

| Método | Opções |
| --- | --- |
| todos | -ksp_rtol 1e-8 -ksp_atol 0 -ksp_max_it 20000 -ksp_divtol 1e30 |
|   | -pc_type none (ou jacobi) -ksp_converged_reason |
| GC | -ksp_type cg -ksp_norm_type unpreconditioned |
| GMRES | -ksp_type gmres -ksp_gmres_restart 30 |
|   | -ksp_gmres_modifiedgramschmidt |
| BiCGSTAB | -ksp_type bcgs |
| GMRES e BiCGSTAB | -ksp_pc_side e -ksp_norm_type, conforme a sua implementação |

Confira com -ksp_view como o PETSc aplica o pré-condicionador e qual norma usa no critério de parada; nem toda combinação é aceita por todos os métodos. Se não for possível reproduzir exata- mente a sua escolha, documente a diferença. Sem as opções da tabela, o PETSc usa Gram-Schmidt clássico (com refinamento quando necessário) no GMRES e interrompe o método quando o resíduo cresce mais de 105 vezes (-ksp_divtol). Aopção -ksp_monitor_true_residual imprime o resíduo a cada iteração e -log_view mostra o tempo gasto em KSPSolve.

Para cada caso, compare com a sua implementação o número de iterações, o resíduo verdadeiro final e os tempos total e por iteração.

## Questões para discussão (Partes 1 a 5)

- 1. Por que o pré-condicionador de Jacobi não altera o número de iterações nos problemas modelo? O que acontece nas matrizes do repositório, e por quê?

- 2. No problema 1, compare GC, GMRES(30) e BiCGSTAB em iterações, produtos matriz-vetor e tempo. Por que o GMRES(30) precisa de mais iterações que o GC, se o GMRES minimiza o re- síduo?

- 3. Como o número de iterações cresce com N em cada método? Relacione com o crescimento de κ2(A) visto nas aulas.

- 4. O que muda com ℓ no GMRES(ℓ)? Houve casos em que o GMRES(ℓ) não convergiu? Explique.

- 5. Compare o histórico do resíduo do BiCGSTAB com o do GMRES. Em algumcaso o resíduo calcu- lado pelo BiCGSTAB ficou muito diferente do verdadeiro? Por que isso pode acontecer?


- 6. Em que casos o número de iterações da sua implementação difere do PETSc, e por quê? Que fração do tempo do seu código e do PETSc é gasta no produto matriz-vetor?

## Parte 6 – Autovalores: implementações

As Partes 6 a 8 tratam do problema simétrico Ax λx, com A esparsa, e do cálculo de poucos au- topares nos extremos do espectro. Elas reaproveitam a estrutura CSR, o produto matriz-vetor e os

resolvedores das Partes 1 e 2. Em todos os métodos, um autopar aproximado (θ,x), com aceito quando

1, é

## Métodos obrigatórios

- 1. Método da potência normalizado, como quociente de Rayleigh θk xT k Axk como estimativa do autovalor.

- 2. Potência inversa com deslocamento σ. Cada iteração resolve (A−σI)y xk−1 com o seu GC (quando A−σI é SPD) ou com o seu GMRES(ℓ), e não por fatoração. A tolerância relativa τ do resolvedor interno é um parâmetro.

- 3. Lanczos com reortogonalização completa, pelo algoritmo abaixo. Os autovalores de Tk (valores de Ritz θi ) e seus autovetores si são calculados pela rotina dstev da LAPACK, que usa o algo- ritmo QR implícito para matrizes tridiagonais. Os vetores de Ritz são ° ° yi Vksi . Para decidir a parada, use a estimativa °Ayi yi ° 2 βk k si que não custa nenhum produto matriz-vetor, e confira-a contra o resíduo verdadeiro ao final.

Tk

= tridiag(β1:k−1,

α1:k, β1:k−1)

Além do produto matriz-vetor, o passo j custa cerca de 9m flops e a reortogonalização acrescenta cerca de 4jm. É preciso guardar toda a base Vk. Sem a reortogonalização, os vetores vj perdem a ortogonalidade em precisão finita e aparecem cópias espúrias de autovalores já convergidos.

## Opcional (escolha um)

- (a) Iteração do quociente de Rayleigh, como GMRES como resolvedor interno.

- (b) Iteração de subespaço com p vetores, para os p maiores autovalores.

- (c) AlgoritmoQR com deslocamento deWilkinson e deflação para a tridiagonal Tk, comparado com a dstev nos valores e no número de iterações QR por autovalor.


## Parte 7 – Autovalores: experimentos

## E1. Validação com espectro conhecido

Os autovalores da matriz do problema 1 são

Para N= 32 e N= 64, calcule λmax pela potência, λmin pela potência inversa com σ 0 (GC interno com τ 10−10) e os 5 maiores e os 5 menores autovalores pelo Lanczos. Use o mesmo x0 aleató- rio em todos os métodos. Reporte iterações, produtos matriz-vetor (inclusive os do resolvedor in-

terno), tempo e o erro

k ln(tol)/lnr, sendo r a razão entre os dois autovalores relevantes. Como λjk =λkj , há autovalores duplos: diga quantas cópias de cada um os métodos encontram.

Compare as iterações da potência e da inversa com a estimativa

## E2. Tolerância do resolvedor interno

No problema 1 comN=64 eσ= 0, rode a potência inversa comτ {10−1,10−2,10−4,10−6,10−8,10−10} e nomáximo 100 iterações externas. Para cada τ, reporte as iterações externas, o total de iterações do GC, o resíduo relativo final e o erro relativo do autovalor. Relacione o resíduo atingido com τ e o erro do autovalor com o resíduo.

## E3. Condicionamento e pré-condicionamento

Este experimento liga as duas metades do trabalho. Para o problema 1comN∈ {32,64} e para a matriz SPD do repositório que você usou na Parte 4, calcule comos seus métodos:

- λmax, os 5 menores autovalores e κ2(A);

- os mesmos valores para D−1/2AD−1/2, que tem o mesmo espectro de D−1A, com D= diag(A).

Compare κ2 com o número de iterações do GC, com e sem Jacobi, obtido nas Partes 3 e 4. Explique, a partir dos espectros, por que o Jacobi não muda nada no problema 1 e o que ele faz com o espectro da matriz do repositório. O número de condição sozinho explica as iterações?

## Parte 8 – Autovalores: problemas maiores e comparação com a LAPACK

- 1. Para o problema 1 com N= 32 (m 1024), calcule o espectro completo na forma densa com a rotina dsyevd da LAPACK. Compare o tempo com o do Lanczos para 5 autovalores e o erro com a fórmula exata.

- 2. Para o problema 1 com N 128 e para a matriz SPD do repositório, calcule os 5 maiores e os 5 menores autovalores com as suas implementações. Reporte iterações, produtos matriz-vetor, tempo e resíduo.

## Opcional: comparação com o SLEPc

Repita o item 2 como EPS do SLEPc, usando a mesma matriz CSR, como na Parte 5. Rode as configu- rações da tabela abaixo e reporte iterações, produtos matriz-vetor (com -log_view, linha MatMult), tempo e resíduo.


```
EPS eps;
EPSCreate(PETSC_COMM_SELF, &eps);
EPSSetOperators(eps, A, NULL); /* A criada com MatCreateSeqAIJWithArrays */
EPSSetProblemType(eps, EPS_HEP); /* problema simetrico */
EPSSetFromOptions(eps); /* metodo e opcoes pela linha de comando */
EPSSolve(eps);
EPSGetConverged(eps, &nconv);
for (i = 0; i < nconv; i++) {
EPSGetEigenpair(eps, i, &lambda, NULL, NULL, NULL);
EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &erro);
}
Configuração Opções
todas -eps_nev 5 -eps_tol 1e-8
Krylov–Schur, maiores -eps_type krylovschur -eps_largest_real
Krylov–Schur, menores -eps_type krylovschur -eps_smallest_real
Lanczos -eps_type lanczos -eps_lanczos_reorthog full
inversa com desl., LU -st_type sinvert -eps_target 0
inversa com desl., GC -st_type sinvert -eps_target 0 -st_ksp_type cg
-st_pc_type jacobi -st_ksp_rtol 1e-10
potência -eps_type power -eps_nev 1 -eps_max_it 100000
```

Com a opção -st_type sinvert, o SLEPc resolve, por padrão, cada sistema por LU (-st_ksp_type preonly -st_pc_type lu). Ao trocar por um método iterativo, a tolerância interna tem de ser bem menor que a de autovalores, como no E2. Sem a chamada a EPSSetProblemType com EPS_HEP, o SLEPc trata a matriz como não simétrica.

## Questões para discussão (Partes 6 a 8)

- 1. Por que a potência precisa de milhares de iterações para λmax do problema 1, enquanto a inversa com σ=0 converge em cerca de 20 para λmin? Como cada número cresce com N?

- 2. Com k produtos matriz-vetor, a potência e o Lanczos trabalham no mesmo espaço de Krylov Kk(A,x0). Por que o Lanczos extrai aproximações muito melhores? Qual a relação entre o Lanc- zos e o GC?

- 3. No E2, qual tolerância interna dá o menor custo total para tol =10−8? Por que o erro do autovalor é muito menor que o resíduo?

- 4. As estimativas do E3 explicam quantitativamente as iterações do GC das Partes 3 e 4? Onde a estimativa κ2 falha e o que a distribuição dos autovalores acrescenta?

- 5. O seu Lanczos encontra as duas cópias dos autovalores duplos do problema 1? E o SLEPc, se você fez a parte opcional? Por que um método com um único vetor inicial tem dificuldade com autovalores múltiplos?

## Entregáveis

- 1. Código: arquivo compactado ou link para repositório, com Makefile (ou CMakeLists.txt; em Python, requirements.txt) e um README explicando como compilar e reproduzir cada tabela e figura do relatório. Não inclua as matrizes; indique os nomes e os links de download.


- 2. Relatório em PDF com no máximo 15 páginas, contendo: descrição e justificativa das matrizes; descrição das implementações (estruturas de dados e decisões tomadas); máquina usada (pro- cessador, memória, compilador e flags, versões do PETSc e da LAPACK, e do SLEPc, se usado); resultados das Partes 1 a 8 em tabelas e gráficos; respostas às questões para discussão; conclu- sões.

- 3. Apresentação de 15minutos em aula, focada nos resultados e na discussão, não na descrição do código.

## Critérios de avaliação

- Correção das implementações, verificada pelo resíduo verdadeiro e pela comparação com os re- sultados do PETSc e da LAPACK.

- Qualidade e reprodutibilidade dos experimentos: metodologia de medição, descrição da máquina e README.

- Análise: relação entre os resultados e a teoria (custo por iteração, memória, condicionamento, efeito do pré-condicionador e do reinício, espectro e taxas de convergência dos métodos de auto- valores), e não apenas a apresentação dos números.

- Clareza do relatório e da apresentação.

## Apêndice A – Fazendo o trabalho emPython

O trabalho pode ser feito emPython, comasmesmasPartes 1 a 8 e osmesmosentregáveis. Bibliotecas

necessárias: numpy, scipy, matplotlib, petsc4py e, para a parte opcional, slepc4py. No relatório,

informe as versões usadas.

## Matrizes

No lugar do módulo csr.c, use a csr_matrix do SciPy: leitura com scipy.io.mmread e produto matriz-vetor com A @ x. Monte os problemas modelo construindo diretamente os vetores indptr, indices e data, sem passar pela forma densa nem usar atalhos como scipy.sparse.kron, e crie a matriz com:

```
A = scipy.sparse.csr_matrix((data, indices, indptr), shape=(m, m))
```

As operações com vetores (produto interno, y y

normas) podem usar o NumPy. O pré-

condicionador e os métodos devem ser escritos por você: não use scipy.sparse.linalg.


## Equivalência das ferramentas

| C/C++ | Python |
| --- | --- |
| csr_read_mtx | scipy.io.mmread(...).tocsr() |
| csr_matvec | A @ x |
| vetores do CSR | A.indptr, A.indices, A.data |
| MatCreateSeqAIJWithArrays PETSc.Mat().createAIJ(size=..., csr=(...)) |   |
| KSPSetFromOptions e opções métodos do objeto KSP (exemplo abaixo) |   |
| EPS do SLEPc | SLEPc.EPS() (exemplo abaixo) |
| dstev | scipy.linalg.lapack.dstev |
| dsyevd | scipy.linalg.eigh(..., driver="evd") |
| std::chrono, clock_gettime time.perf_counter |   |

O petsc4py e o slepc4py podem ser instalados pelo conda,

```
conda install -c conda-forge petsc4py slepc4py
```

ou, emDebian/Ubuntu, pelos pacotes python3-petsc4py e python3-slepc4py. Pelo pip, eles com- pilam o PETSc e o SLEPc a partir do código-fonte, o que demora.

## Exemplo

```
import numpy as np, scipy.io
import petsc4py; petsc4py.init()
from petsc4py import PETSc
A = scipy.io.mmread("wang3.mtx").tocsr() # simetricas ja vem completas
A.sort_indices()
b = A @ np.ones(A.shape[0]) # so para montar b
M = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))
bb = PETSc.Vec().createWithArray(b); xx = bb.duplicate(); xx.set(0.0)
ksp = PETSc.KSP().create()
ksp.setOperators(M)
ksp.setType("bcgs") # "cg", "gmres" ou "bcgs"
ksp.getPC().setType("jacobi") # "none" ou "jacobi"
ksp.setTolerances(rtol=1e-8, atol=0.0, divtol=1e30, max_it=20000)
ksp.solve(bb, xx)
print(ksp.getIterationNumber(), ksp.getConvergedReason())
```

No GMRES e no BiCGSTAB, ajuste o lado do pré-condicionador e a norma com ksp.setPCSide e ksp.setNormType, conforme a sua implementação, antes de ksp.solve. Para o GC e para oGMRES, acrescente:

```
\# GC
ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
# GMRES
ksp.setGMRESRestart(30)
PETSc.Options()["ksp_gmres_modifiedgramschmidt"] = None
ksp.setFromOptions()
```

Parte opcional da Parte 8, autovalores como slepc4py, usando a mesma matriz M:

```
import slepc4py; slepc4py.init()
from slepc4py import SLEPc
eps = SLEPc.EPS().create()
eps.setOperators(M)
eps.setProblemType(SLEPc.EPS.ProblemType.HEP) # simetrico
```


```
eps.setDimensions(nev=5)
eps.setTolerances(tol=1e-8, max_it=100000)
eps.setType("krylovschur") # "lanczos", "power", ...
eps.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_REAL)
eps.solve()
for i in range(eps.getConverged()):
lam = eps.getEigenvalue(i).real
erro = eps.computeError(i, SLEPc.EPS.ErrorType.RELATIVE)
print(lam, erro)
```

Para a inversa com deslocamento, troque a escolha dos autopares por:

```
eps.setTarget(0.0)
eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
st = eps.getST(); st.setType("sinvert")
ksp = st.getKSP(); ksp.setType("cg"); ksp.getPC().setType("jacobi")
ksp.setTolerances(rtol=1e-10)
```
