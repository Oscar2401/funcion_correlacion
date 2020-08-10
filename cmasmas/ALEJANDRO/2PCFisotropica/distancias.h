#include <cmath>

// distancia euclidea
template <typename TDG1>
TDG1 dist(TDG1 x, TDG1 y, TDG1 z){
    return sqrt(x*x + y*y + z*z);
}

// distancia aproximada método a min b max
template <typename TDG>
TDG distancia(TDG x, TDG y){
    TDG min, max, aprox;
    if (x < 0){
        x*=-1;
    }
    if (y < 0){
        y*=-1;
    }

    if (x < y){
        min = x;
        max = y;
    }else{
        min = y;
        max = x;
    }

    aprox = (max*1007.0) + (min*441.0);
    //condicion para ajuste y mejora de precisión
    if (max < (16.0*min)){ //equivalente a decir que si max < 16*min
        aprox -= max*40.0; // equivalente a aprox - 5*max/128
    }
    //se añade un 512 para un apropiado redondeo de las cifras
    // se divide entre 1024 con >> 5
    return ((aprox+512.0)/1024.0);
}
