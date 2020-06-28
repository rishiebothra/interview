FROM gcc:8.3

WORKDIR /app/lbm_mrt

RUN apt-get update && apt-get install -y make

#RUN apt-get update && apt-get install -y emacs24


