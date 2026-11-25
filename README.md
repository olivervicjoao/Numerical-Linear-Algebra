how to compile

cmake -S scripts -B build -DUSE_SLEPC=ON
cmake --build build

./build/nla
