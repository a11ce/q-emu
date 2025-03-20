Source cuda-env or otherwise get nvcc, then compile with `nvcc -std=c++17 gpu.cu -o q-emu`. 

If you want to compile the cpu tests, run `clang++ -std=c++17 main.cpp -o cpu-tests`. This is not needed for performance measurement.

The first argument is 0 for CPU or 1 for GPU. The second argument is 0 to tensor and run or 1 to just tensor.

To compare tensoring speeds, time `./q-emu 0 1` and `./q-emu 1 1`. We used `hyperfine './q-emu 0 1' './q-emu 1 1'`.

To compare the speeds for tensoring and running the example circuit 100 times, time `./q-emu 0 0` and `./q-emu 1 0`.

To get the time taken for just running the circuit 100 times, subtract the previous two numbers per implementation.

See main.cpp for other usage examples (for the CPU version).
