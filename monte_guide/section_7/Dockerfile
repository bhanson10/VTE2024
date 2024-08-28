FROM gcc:latest
WORKDIR /usr/src/app
COPY gbees.c .
RUN gcc -shared -o gbees.so gbees.c -fPIC
