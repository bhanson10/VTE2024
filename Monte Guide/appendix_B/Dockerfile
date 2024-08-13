FROM gcc:latest
WORKDIR /usr/src/app
COPY . .
RUN gcc -shared -o hello_world.so hello_world.c -fPIC
