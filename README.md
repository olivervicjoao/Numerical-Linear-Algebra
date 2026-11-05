# Álgebra Linear Computacional — Métodos diretos e LAPACKE

Implementações em **C++17** desenvolvidas para o Trabalho Prático 1 de Álgebra Linear Computacional do Programa de Pós-Graduação em Modelagem Computacional da UFJF (2026). O projeto compara LU, Cholesky e Thomas próprios com rotinas equivalentes da **LAPACKE**, e investiga a redução de largura de banda por **Reverse Cuthill–McKee (RCM)**.

## 1. Métodos implementados

| Estrutura | Implementação própria | Referência LAPACKE |
|---|---|---|
| Matriz densa | LU com pivoteamento parcial | `dgetrf` + `dgetrs` |
| Matriz densa SPD | Cholesky | `dpotrf` + `dpotrs` |
| Matriz em banda | LU **sem** pivoteamento | `dgbtrf` + `dgbtrs` |
| Matriz em banda SPD | Cholesky em banda | `dpbtrf` + `dpbtrs` |
| Tridiagonal | Thomas (três vetores) | `dgtsv` |

A matriz pode ser representada de modo denso ou por **armazenamento compacto em banda**, ambos em ordem por colunas (*column-major*). As fatorações próprias usam laços explícitos e separam a etapa de fatoração da resolução; Thomas também mantém as etapas separadas. A ordenação RCM é uma implementação própria baseada em listas de adjacência e BFS. As implementações de LU em banda e Thomas **não têm pivoteamento** e podem falhar em matrizes inadequadas.

## 2. Dependências e compilação

- Compilador C++ com suporte a **C++17**, `make` e Bash;
- **LAPACKE** e **OpenBLAS** (com cabeçalhos e bibliotecas para ligação);
- Python 3 com **NumPy, SciPy e Matplotlib** para validações e gráficos; para a geração das figuras reproduzidas no relatório, são necessários esses pacotes.

Na raiz deste repositório (o Makefile está em `scripts/`, não na raiz):

```bash
make -f scripts/Makefile clean all
./nla --help
```

Por padrão, a compilação usa `-O3 -std=c++17`, inclui `include/` e liga `-llapacke -lopenblas -lm`. No cluster usado nos experimentos finais, a OpenBLAS disponibilizava também a interface LAPACKE; foi utilizado:

```bash
module load openblas/0.3.21
export OPENBLAS_ROOT=/opt/ohpc/pub/libs/gnu13/openblas/0.3.21
make -f scripts/Makefile clean all \
  INC="-I./include -I${OPENBLAS_ROOT}/include" \
  LDLIBS="-L${OPENBLAS_ROOT}/lib -Wl,-rpath,${OPENBLAS_ROOT}/lib -lopenblas -lm"
```

O caminho de instalação do módulo é **específico do cluster**; em outra máquina, adapte `OPENBLAS_ROOT` e as bibliotecas disponíveis. Para realizar uma comparação **serial**:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
```

## 3. Obter as matrizes (não incluídas no repositório)

Os arquivos `.mtx` são disponibilizados pelo **NIST Matrix Market**, coleção Harwell–Boeing. Baixe o arquivo *Compressed MatrixMarket format* em cada página, descompacte-o e organize-o assim:

| Matriz | Origem | Pasta esperada | Página de download |
|---|---|---|---|
| `bcsstk11.mtx` | Rigidez de vagão de minério (SPD) | `matrix/SPD/` | <https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc1/bcsstk11.html> |
| `bcsstk14.mtx` | Cobertura do Omni Coliseum (SPD) | `matrix/SPD/` | <https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc2/bcsstk14.html> |
| `sherman4.mtx` | Reservatório de petróleo, IMPES com barreiras | `matrix/NS/` | <https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/sherman/sherman4.html> |
| `orsreg_1.mtx` | Reservatório de petróleo, grade 21×21×5 | `matrix/NS/` | <https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/oilgen/orsreg_1.html> |

Exemplo de preparação da estrutura, depois de obter os quatro arquivos `.mtx.gz`:

```bash
mkdir -p matrix/SPD matrix/NS
# Descompacte bcsstk11.mtx.gz e bcsstk14.mtx.gz em matrix/SPD/
# Descompacte sherman4.mtx.gz e orsreg_1.mtx.gz em matrix/NS/
```

O programa **não usa os vetores do lado direito dos problemas físicos originais**: constrói uma solução exata `x* = (1,...,1)^T` e define `b = A x*`, permitindo calcular erros. Para repetir a exploração opcional de 19 matrizes, obtenha também os arquivos adicionais listados em `scripts/run_experimentos_finais.bash` nas pastas `matrix/SPD/` e `matrix/NS/`. Nenhuma matriz de terceiros é redistribuída no ZIP.

## 4. Reproduzir os experimentos do relatório

A campanha final no cluster foi executada em **um processo serial**, no nó `compute-0-0` (AMD Opteron 6272, Rocky Linux 9.6, GCC 13.2.0, OpenBLAS 0.3.21), com um nó reservado pelo OpenPBS usando `place=exclhost`. As medições usam `std::chrono::steady_clock` e são resumidas por **medianas**.

### Execução completa

```bash
bash scripts/run_experimentos_finais.bash
```

O comando roda: Thomas para `n=1000,10000,100000,1000000` (**20 repetições**); LU/Cholesky densas para `n=250,500,1000,2000,4000` (**5 repetições**); e as quatro matrizes reais (**20 repetições** cada). O script prepara os gráficos básicos, valida os CSVs e registra o ambiente. A execução integral inclui casos de fatoração densa e pode exigir várias horas no hardware de referência.

Para rodar apenas partes, sem repetir o resto:

```bash
bash scripts/run_experimentos_finais.bash --only-thomas
bash scripts/run_experimentos_finais.bash --only-dense
bash scripts/run_experimentos_finais.bash --only-selected
# Opcional: campanha exploratória com as 19 matrizes (requer seus .mtx)
bash scripts/run_experimentos_finais.bash --only-nist
```

### Execução isolada do executável

```bash
./nla --matrix matrix/SPD/bcsstk11.mtx --output tests/exemplo_bcsstk11 --repeats 20 --spy
./nla --thomas-scale 1000000 --output tests/exemplo_thomas --repeats 20
./nla --dense-scale 1000 --output tests/exemplo_denso --repeats 5
```

O argumento `--spy` chama a rotina de visualização, que requer SciPy e Matplotlib. O script `scripts/spy_plot.py` lê arquivos Matrix Market por meio de `scipy.io.mmread` sem usar o argumento `spmatrix`, para compatibilidade com diferentes versões do SciPy.

## 5. Arquivos e correspondência com as tabelas do relatório

| Parte / tabela | Saída principal |
|---|---|
| Caracterização (dimensão, nnz, densidade, bandas, condição) | `tests/nist_final_20/all_matrix_metadata.csv` |
| Tempos por método e erros das matrizes reais | `tests/nist_final_20/all_summary.csv` |
| Efeito de RCM: densa versus banda reordenada | `tests/nist_final_20/all_part3_3.csv` e `all_part3_comparison.csv` |
| Construção do grafo, ordenação, exportação, leitura | `tests/nist_final_20/all_preprocessing.csv` e `all_preprocessing_raw.csv` |
| Escalabilidade de Thomas | `tests/thomas_scaling/all_summary.csv`, `scaling_slopes.csv` |
| Escalabilidade LU e Cholesky | `tests/dense_scaling/all_summary.csv`, `scaling_slopes.csv` |
| LU/Cholesky em `n=4000` | `tests/dense_scaling/dense_scale_4000/summary.csv` |
| Repetições individuais e medidas de dispersão | `all_raw.csv`, `timing_variability.csv` em cada grupo |

Em `all_summary.csv`, o campo `variant` diferencia matrizes `original` e `rcm`; `method` distingue implementações próprias e LAPACKE. Os resultados das 19 matrizes exploratórias, quando executados com `--only-nist`, são gravados em `tests/nist_repeticao_5/`, **não** em `tests/nist_final_20/`. Os exemplos exploratórios de estabilidade citados no texto foram obtidos em uma campanha anterior, não repetida pelo comando padrão.

Os CSVs trazem `factor_s`, `solve_s`, `total_s`, `factor_gflops`, `forward_error_inf`, `backward_error_inf`, `kappa1_est`, `kappa_u` e `lu_growth`. `factor_gflops` é calculado a partir de contagens **teóricas**, não de contadores de hardware. Em `dgtsv`, usa-se o tempo **total** (não há fatoração/resolução separadas na chamada). A estimativa de pré-processamento RCM inclui exportação e nova leitura e não é o tempo isolado da heurística.

## 6. Reproduzir as figuras em português

Instale NumPy, SciPy e Matplotlib em seu ambiente Python e execute, depois dos experimentos:

```bash
python3 -m pip install numpy scipy matplotlib
python3 scripts/gerar_figuras_relatorio_ptbr.py --root . --output figuras_lista_1
```

O script grava os cinco PDFs **com os mesmos nomes esperados pelo LaTeX do relatório**:

- `escalabilidade_thomas_cluster.pdf`;
- `escalabilidade_lu_cholesky_cluster.pdf`;
- `desempenho_fatoracoes_densas_ptbr.pdf`;
- `rcm_estruturais_ptbr.pdf`;
- `rcm_reservatorios_ptbr.pdf`.

Também gera PNGs para inspeção. Para as figuras RCM são necessários os quatro arquivos `.mtx` originais nas pastas acima e os `matrix_rcm.mtx` produzidos pela campanha final. Para reproduzir a tabela do relatório, use os CSVs indicados na seção 5. O script imprime mensagens claras quando algum arquivo necessário não estiver disponível.

Os scripts existentes `scripts/plot_results.py` e `scripts/spy_plot.py` continuam disponíveis para as visualizações de diagnóstico geradas durante os benchmarks. O novo script gera a versão **padronizada em português** para o relatório.

## 7. Validação e observações

```bash
python3 scripts/validate_benchmarks.py tests/thomas_scaling tests/dense_scaling tests/nist_final_20
```

O critério de correção utiliza uma solução conhecida e registra erros **direto e de retroação**. Problemas mal condicionados podem produzir erros diretos elevados mesmo com resíduos pequenos. O reordenamento RCM não garante redução de banda: no conjunto final, por exemplo, a banda de `bcsstk14` aumentou. As matrizes reais têm origens e escalas diferentes; extrapolações de tempos para outros computadores dependem da configuração de BLAS, CPU, memória e compilador.

## 8. Estrutura do projeto

```text
include/            Cabeçalhos das classes, métodos e experimentos
src/                Implementações em C++
scripts/            Makefile, execução, validação e geração de gráficos
main.cc             Interface de linha de comando
README.md           Este guia
```

Referências para teoria e análise: Golub & Van Loan, *Matrix Computations* (2013); Saad, *Iterative Methods for Sparse Linear Systems* (2003); Trefethen & Bau, *Numerical Linear Algebra* (1997).
