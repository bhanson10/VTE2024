FROM ubuntu:latest
RUN apt-get update && apt-get install -y build-essential liblapack-dev libblas-dev liblapacke-dev
WORKDIR /app
COPY gbees.c .
RUN gcc -shared -o gbees.so gbees.c -fPIC -llapacke -llapack -lblas
