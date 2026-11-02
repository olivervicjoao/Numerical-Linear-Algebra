#!/bin.bash


export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

bash scripts/run_all_tests.sh \
    --output tests/thomas_scaling \
    --repeats 20 \
    --thomas-scale



bash scripts/run_all_tests.sh \
    --output tests/dense_scaling \
    --repeats 5 \
    --dense-scale



bash scripts/run_all_tests.sh \
    --output tests/nist_final_20 \
    --repeats 20 \
    matrix/SPD/bcsstk11.mtx \
    matrix/SPD/bcsstk14.mtx \
    matrix/NS/sherman4.mtx \
    matrix/NS/orsreg_1.mtx
