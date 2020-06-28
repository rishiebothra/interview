
# Instructions:

## 1.  Compile the code on the terminal
```
make
```
## 2.  Execute
```
./a.out input_data output_data
```

## 3.  Docker:
```
## buld the docker
docker build -f ./Dockerfile -t devops:latest .

## enter into docker environment
docker run --privileged --rm -ti -v $PWD:/app/lbm_mrt devops:latest /bin/bash

## running the docker
docker run --privileged --rm -ti -v $PWD:/app/lbm_mrt devops:latest make && ./a.out input_data output_data
```






