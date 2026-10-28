2026.3 

## Universidade Federal de Juiz de Fora 

Programa de Pós-graduação em Modelagem Computacional 

# **Álgebra Linear Computacional** 

## Trabalho prático: métodos diretos e comparação com a LAPACK 

Prof. José J. Camata 

### **Modalidade:** individual. 

**Entrega do código e do relatório:** 29/10/2026 **Apresentação em aula:** 29/10/2026, 10 minutos por aluno. 

# **Objetivo** 

Implementar os métodos diretos vistos na disciplina para matrizes densas, em banda e tridiagonais, aplicá-los a matrizes reais de um repositório público e comparar exatidão e desempenho com as rotinas correspondentes da LAPACK. O trabalho também avalia o efeito de um reordenamento (Cuthill– McKee reverso) sobre a largura de banda e sobre o custo da fatoração. 

# **Ambiente** 

- Linguagem: C ou C++. Quem preferir Python deve seguir o Apêndice A, que dá as rotinas equivalentes. 

- LAPACK pela interface C **LAPACKE** , com OpenBLAS ou a BLAS de referência. Em Debian/Ubuntu: pacotes `liblapacke-dev` e `libopenblas-dev` . 

- Leitura de matrizes no formato Matrix Market com as rotinas `mmio.c` / `mmio.h` do NIST: `https: //math.nist.gov/MatrixMarket/mmio-c.html` . 

- Reordenamento RCM pela Boost Graph Library ( `boost::cuthill_mckee_ordering` , cabeçalho `boost/graph/cuthill_mckee_ordering.hpp` ; pacote `libboost-graph-dev` ). Passando um iterador reverso, a rotina devolve a ordem RCM. O RCM pode ser implementado também. 

### Exemplo de compilação: 

```
g++-O2-std=c++17tp.cppmmio.c-otp-llapacke-lopenblas
```

Todas as medições de tempo devem ser feitas com a BLAS em uma única thread, para que a comparação seja entre códigos seriais: 

```
exportOPENBLAS_NUM_THREADS=1
```

Use `std::chrono` ou `clock_gettime` e reporte a mediana de pelo menos 5 execuções. 

As bibliotecas acima só podem ser usadas para leitura dos arquivos, para o reordenamento e como referência de comparação. As fatorações e substituições das Partes 2 e 3 devem ser implementadas por você. 

# **Escolha das matrizes** 

Cada aluno escolhe as próprias matrizes na **SuiteSparse Matrix Collection** , `https://sparse.tamu. edu` . A página permite filtrar por dimensão, número de não nulos, simetria e definição positiva, e 

1 

cada matriz tem um link de download no formato Matrix Market ( `.mtx` ). O repositório Matrix Market do NIST, `https://math.nist.gov/MatrixMarket/` , também pode ser usado. 

Escolha **pelo menos três** matrizes reais e quadradas, com 500 _≤ m ≤_ 5000 (para que a forma densa caiba na memória), que juntas cubram os seguintes casos: 

- (a) uma matriz **simétrica positiva definida** (para Cholesky); 

- (b) uma matriz **não simétrica** , para a LU com pivoteamento parcial; 

- (c) uma matriz esparsa cuja largura de banda **diminua de forma significativa** com o RCM (Parte 3). Pode ser a mesma do item (a). 

Justifique a escolha no relatório (origem do problema, tamanho, propriedades). Alguns exemplos que atendem aos critérios, apenas como referência: 

|Matriz|Tipo|_m_|nnz|SPD|_κ_estimado|
|---|---|---|---|---|---|
|`HB/gr_30_30`|problema 2D|900|7 744|sim|1,9_×_10<sup>2</sup>|
|`HB/1138_bus`|rede elétrica|1 138|4 054|sim|8,6_×_10<sup>6</sup>|
|`HB/nos7`|problema 2D/3D|729|4 617|sim|2,4_×_10<sup>9</sup>|
|`HB/bcsstk14`|estrutural|1 806|63 454|sim|1,2_×_10<sup>10</sup>|
|`HB/bcsstk15`|estrutural|3 948|117 816|sim|6,5_×_10<sup>9</sup>|
|`HB/orsirr_1`|dinâmica dos fluidos|1 030|6 858|não|7,7_×_10<sup>4</sup>|
|`HB/west0989`|processo químico|989|3 518|não|9,9_×_10<sup>11</sup>|



Os valores de _κ_ são as estimativas informadas no próprio repositório. 

**Atenção:** em matrizes marcadas como `symmetric` no cabeçalho do arquivo `.mtx` , só a parte triangular inferior está armazenada. Ao montar a matriz é preciso espelhar os coeficientes fora da diagonal. 

# **Parte 1 – Leitura e caracterização** 

Para cada matriz escolhida, leia o arquivo `.mtx` , monte a forma densa (coluna-maior, como a LAPACK espera) e reporte em uma tabela: 

- dimensão _m_ , número de não nulos nnz e densidade nnz/ _m_<sup>2</sup> ; 

- se a matriz é simétrica; 

- larguras de banda _bI_ e _bS_ ; 

- estimativa de _κ_ 1( _A_ ) usando `dlange` seguida de `dgecon` (ou `dpocon` , no caso SPD); 

- uma figura do padrão de esparsidade ( _spy plot_ ), que pode ser gerada com qualquer ferramenta. 

# **Parte 2 – Implementações** 

Implemente, separando sempre a fatoração da solução: 

1. **LU com pivoteamento parcial** para matrizes densas, _PA = LU_ , com _L_ e _U_ sobrescrevendo _A_ e o vetor de pivôs guardado à parte; substituições direta e inversa. 

2. **Cholesky** para matrizes densas SPD, _A = LL_<sup>_T_</sup> , com detecção de pivô não positivo. 

2 

3. **LU em banda sem pivoteamento** , com a matriz guardada em armazenamento em banda (no máximo ( _bI +bS +_ 1) _m_ posições) e substituições restritas à banda. Para matrizes SPD, implemente também **Cholesky em banda** . 

4. **Algoritmo de Thomas** para sistemas tridiagonais, usando apenas três vetores. 

Valide cada implementação em matrizes pequenas, comparando com a solução da LAPACK, antes de usá-la nas matrizes do repositório. 

# **Parte 3 – Reordenamento RCM e fatoração em banda** 

Para a matriz do caso (c): 

1. Obtenha a permutação RCM com a Boost (ou o SciPy, em Python) a partir do grafo de adjacência de _A_ (um vértice por linha, uma aresta para cada _ai j̸ =_ 0, _i̸ = j_ ) e forme _A_<sup>˜</sup> _= PAP_<sup>_T_</sup> . O sistema reordenado é _A_<sup>˜</sup> ( _Px_ ) _= Pb_ . 

2. Reporte _bI_ , _bS_ e o armazenamento em banda antes e depois do reordenamento, e mostre os dois padrões de esparsidade lado a lado. 

3. Resolva o sistema reordenado com a sua fatoração em banda (LU ou Cholesky) e compare o tempo com a fatoração densa da Parte 2 e com a previsão dada pelas contagens de operações. 

# **Parte 4 – Comparação com a LAPACK** 

Compare cada implementação com a rotina correspondente: 

|Sua implementação|LAPACK|Flops da fatoração|
|---|---|---|
|LU densa com pivoteamento|`dgetrf`/|`dgetrs`<br>2<br>3<sup>_m_3</sup><br><br>|
|Cholesky densa|`dpotrf`/|`dpotrs`<br>1<br>3<sup>_m_3</sup>|
|LU em banda|`dgbtrf`/|`dgbtrs`<br>_≈_2_m bIbS_<br>|
|Cholesky em banda|`dpbtrf`/|`dpbtrs`<br>_≈m b_<sup>2</sup>|
|Thomas|`dgtsv`|_≈_8_m_(fatoração e solução)|



Observe que `dgbtrf` faz pivoteamento parcial e por isso exige um arranjo com 2 _bI + bS +_ 1 linhas ( `ldab` ); a sua LU em banda não pivota. 

Em todos os testes, tome a solução exata _x =_ (1,...,1)<sup>_T_</sup> e defina _b = Ax_ . Para cada matriz e cada método, reporte: 

- tempo da fatoração e da solução, e a taxa em GFLOP/s calculada com as contagens da tabela; 

- erro de retroação _∥b − A_ ˆ _x∥∞_ �� _∥A∥∞ ∥x_ ˆ _∥∞_ �; 

- erro direto relativo _∥x_ ˆ _− x∥∞_ / _∥x∥∞_ , comparado com o limite _κ_ ( _A_ ) _u_ , onde _u ≈_ 1,1 _×_ 10<sup>_−_16</sup> ; 

- para a LU, o fator de crescimento _ρ =_ max _i j |ui j |_ / max _i j |ai j |_ . 

Além das matrizes escolhidas, faça um estudo de escala com a matriz tridiagonal tridiag( _−_ 1,2, _−_ 1) para _m =_ 10<sup>3</sup> ,10<sup>4</sup> ,...,10<sup>6</sup> (Thomas e `dgtsv` ) e com matrizes densas geradas aleatoriamente para _m_ entre 250 e 4000 (LU e Cholesky). Apresente os tempos em gráficos log-log e compare as inclinações com as previstas. 

3 

## **Questões para discussão** 

1. Por que a LAPACK é mais rápida que a sua LU densa, se os dois códigos fazem o mesmo número de operações? 

2. A diferença de desempenho entre a sua implementação e a LAPACK é a mesma para matrizes densas, em banda e tridiagonais? Explique. 

3. Os erros diretos observados são compatíveis com _κ_ ( _A_ ) _u_ ? Em que matrizes o erro direto é grande mesmo com erro de retroação pequeno? 

4. Quando vale a pena reordenar e usar a fatoração em banda em vez da densa? Considere tempo e memória. 

# **Entregáveis** 

1. **Código** : arquivo compactado ou link para repositório, com `Makefile` (ou `CMakeLists.txt` ; em Python, `requirements.txt` ) e um `README` explicando como compilar e reproduzir cada tabela e figura do relatório. Não inclua as matrizes; indique os nomes e os links de download. 

2. **Relatório** em PDF com **no máximo 10 páginas** , contendo: descrição e justificativa das matrizes; descrição das implementações (estruturas de dados e decisões tomadas); máquina usada (processador, memória, compilador e flags, versão da BLAS); resultados das Partes 1 a 4 em tabelas e gráficos; respostas às questões para discussão; conclusões. 

3. **Apresentação** de **10 minutos** em aula, focada nos resultados e na discussão, não na descrição do código. 

# **Critérios de avaliação** 

- Correção das implementações, verificada pelos erros de retroação e pela comparação com a LAPACK. 

- Qualidade e reprodutibilidade dos experimentos (metodologia de medição, descrição da máquina, `README` ). 

- Análise: relação entre os resultados e a teoria (custos, condicionamento, estabilidade), e não apenas a apresentação dos números. 

- Clareza do relatório e da apresentação. 

# **Apêndice A – Fazendo o trabalho em Python** 

O trabalho pode ser feito em Python, com as mesmas Partes 1 a 4 e os mesmos entregáveis. Bibliotecas necessárias: `numpy` , `scipy` , `numba` , `threadpoolctl` e `matplotlib` . No relatório, informe as versões usadas e a BLAS do NumPy ( `np.show_config()` ). 

4 

## **Implementações com Numba** 

Laços em Python puro são executados pelo interpretador e ficam ordens de grandeza mais lentos que o mesmo laço em C. Comparar um código assim com a LAPACK mediria o custo do interpretador, e não a diferença entre algoritmos. Por isso as implementações das Partes 2 e 3 devem ser escritas com **laços explícitos em funções compiladas com o decorador** `@njit` **do Numba** . O NumPy pode ser usado para alocar vetores e matrizes, mas não para as operações das fatorações e substituições (nada de `np.dot` , `@` ou fatias vetorizadas nos núcleos). 

O Numba compila a função na primeira chamada: chame cada função uma vez com uma matriz pequena antes de medir tempos. 

## **Equivalência das ferramentas** 

|C/C++|Python|
|---|---|
|`mmio.c`|`scipy.io.mmread`|
|`boost::cuthill_mckee_ordering`|`scipy.sparse.csgraph.reverse_cuthill_mckee`|
|`LAPACKE_dgetrf`,`LAPACKE_dpotrf`, ...<br>`export OPENBLAS_NUM_THREADS=1`|`scipy.linalg.lapack.dgetrf`,`.dpotrf`, ...<br>`threadpoolctl.threadpool_limits(1)`|
|`std::chrono`,`clock_gettime`|`time.perf_counter`|



Use as rotinas de `scipy.linalg.lapack` , que chamam a LAPACK diretamente, e não as funções de alto nível ( `scipy.linalg.solve` , `lu_factor` ), que fazem verificações e cópias extras. Todas as rotinas da tabela da Parte 4, e também `dlange` , `dgecon` e `dpocon` , estão disponíveis. 

## **Exemplo** 

```
importnumpyasnp,scipy.io
fromscipy.linalgimportlapack
fromscipy.sparse.csgraphimportreverse_cuthill_mckee
fromthreadpoolctlimportthreadpool_limits
A=scipy.io.mmread("bcsstk14.mtx").tocsr()#simetricasjavemcompletas
perm=reverse_cuthill_mckee(A,symmetric_mode=True)
Ar=A[perm][:,perm]#PAP^T
Ad=Ar.toarray(order="F")#coluna-maior,comoaLAPACK
x=np.ones(Ad.shape[0]);b=Ad@x
withthreadpool_limits(limits=1,user_api="blas"):
c,info=lapack.dpotrf(Ad,lower=1)
xh,info=lapack.dpotrs(c,b,lower=1)
```

## **Cuidados** 

- Ao contrário do `mmio.c` , o `mmread` já devolve as matrizes simétricas completas: não é preciso espelhar os coeficientes. 

- Passe à LAPACK arrays em ordem de coluna ( `order="F"` ou `np.asfortranarray` ). Caso contrário, o SciPy faz uma cópia que entra na medição de tempo. 

- Arranjo em banda, com índices a partir de 0: para `dgbtrf` , `ab` tem 2 _bI + bS +_ 1 linhas e _ai j_ fica em `ab[bI+bS+i-j, j]` ; para `dpbtrf` com `lower=1` , `ab` tem _b +_ 1 linhas e _ai j_ ( _i ≥ j_ ) fica em `ab[i-j, j]` . 

5 

