mkdir build && cd build
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -Dmaat_USE_LIEF=OFF -Dmaat_BUILD_PYTHON_BINDINGS=OFF ..
cmake --build .
