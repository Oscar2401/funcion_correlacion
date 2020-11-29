# Two point correlation function (with BPC)

To compile the code, write:

```json
g++ -o main -fopenmp main.cpp
```

To run the code you need to specify the address of the data files and specify a name for the histogram file. For example:

```json
./main /home/echeveste/mis_trabajos/cosmo_proyect/funcion_correlacion/data/data.dat  full
```

To build the correlation function and graph it, run the python code

```json
python3 graph_2DPF.py 
```
