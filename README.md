# என் -Body Planet Simulation

Just of bunch of big balls on the screen, floating around the way Newton said they would.

Made with C++, SDL, and Emscripten.

To run, ensure you have Emscripten and SDL installed, and then compile with `emcc main.cpp -s WASM=1 -s USE_SDL=2 -O3 -std=c++11 -o index.js -I$(pwd)`. I apologize for this disaster, I will be a decent human being and add a Makefile sometime in the future :)
