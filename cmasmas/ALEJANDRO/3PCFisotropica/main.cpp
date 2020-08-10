#include <iostream>
#include <fstream>
#include <ctime>
#include <string.h>
#include "iso3PCF.h"

using namespace std;

//Prototipos de funciones
void abrir_archivo(string,int, Puntos *);
void guardar_Histogramas(string, int, float ***);
void crear_Histogramas3D(int); // un solo parametro, pues son matrices cuadradas
void eliminar_Histogramas3D(int);
void eliminar_Datos();


//Variables globales
Puntos *DATA, *RAND;
float ***DDD, ***DDR, ***DRR, ***RRR;

// primer argumento del main es el nombre del archivo datos
// segundo argumento del main es el nombre del archivo rand
// tercer argumento es la cantidad de datos a trabajar
// cuarto argumento es el numero de bins
// quinto argumento es la distancia maxima
int main(int argc, char **argv){
    time_t to, tf;
    int N = stoi(argv[3]), nb = stoi(argv[4]), d_max = stof(argv[5]); //cantidad de puntos N, numero de bins nb, d_max
    DATA = new Puntos[N];
    RAND = new Puntos[N];
    //Creo los nombres de los archivos
    string nombre1 = "DDDiso_", nombre2 = "DDRiso_", nombre3 = "DRRiso_", nombre4 = "RRRiso_";
    nombre1.append(argv[3]);
    nombre2.append(argv[3]);
    nombre3.append(argv[3]);
    nombre4.append(argv[3]);
    nombre1 += ".dat";
    nombre2 += ".dat";
    nombre3 += ".dat";
    nombre4 += ".dat";
    //cargo los datos
    abrir_archivo(argv[1],N,DATA);
    abrir_archivo(argv[2],N,RAND);
    //Inicializo con 0 los histogramas
    crear_Histogramas3D(nb);
    //programa principal
    iso3PCF obj(DATA,RAND,N,nb,d_max); // instancio la clase
    to = time(NULL);
    obj.calcular_Histogramas_Puros(DDD,RRR);
    obj.calcular_Histogramas_Mixtos(DDR,DRR);
    tf = time(NULL);
    obj.simetrizar_Histograma(DDD);
    obj.~iso3PCF(); // destruyo el objeto
    eliminar_Datos(); // destruyo structs
    guardar_Histogramas(nombre1,nb,DDD);
    //guardar_Histogramas(nombre2,nb,DDR);
    //guardar_Histogramas(nombre3,nb,DRR);
    //guardar_Histogramas(nombre4,nb,RRR);
    eliminar_Histogramas3D(nb);
    cout << "Calculo realizado en: " << difftime(tf,to) << " seg"<< endl;
    return 0;
}

//toma los datos del archivo y los guarda en un arreglo de structuras.
void abrir_archivo(string nombre_archivo,int cantidad_puntos, Puntos *datos){
    ifstream archivo;
    archivo.open(nombre_archivo.c_str(), ios::in | ios::binary); //le indico al programa que se trata de un archivo binario con ios::binary
    if (archivo.fail()){
        cout << "Error al cargar el archivo " << endl;
        exit(1);
    }
    int c=0;
    float eliminar;
    while (!archivo.eof())
    {
        archivo >> datos[c].x >> datos[c].y >> datos[c].z >> eliminar; 
        c++;
    }
    archivo.close();
}

void guardar_Histogramas(string nombre, int dimension, float ***Matriz){
    int i,j;
    //puntero a todo el array
    float (*puntero_arreglo)[dimension][dimension*dimension] = reinterpret_cast<float(*)[dimension][dimension*dimension]>(Matriz);
    ofstream archivo;
    archivo.open(nombre.c_str(),ios::out | ios::binary);
    if (archivo.fail()){
        cout << "Error al guardar el archivo " << endl;
        exit(1);
    }
    for (i = 0; i < dimension; i++)
    {
        for ( j = 0; j < dimension*dimension; j++)
        {
            archivo << (*puntero_arreglo)[i][j] << " "; 
        }
        archivo << endl;
    }
    archivo.close();
}

//Se crean matrices cuadradas y se inicializan con 0.0 
void crear_Histogramas3D(int dimension){
    int i,j,k;
    //reservación dinamica de memoria para un array de dimensiones cubicas
    DDD = new float**[dimension];
    DDR = new float**[dimension];
    DRR = new float**[dimension];
    RRR = new float**[dimension];
    for (i = 0; i < dimension; i++)
    {
        *(DDD + i) = new float*[dimension];
        *(DDR + i) = new float*[dimension];
        *(DRR + i) = new float*[dimension];
        *(RRR + i) = new float*[dimension];
        for (j = 0; j < dimension; j++)
        {
            *(*(DDD+i)+j) = new float[dimension];
            *(*(DDR+i)+j) = new float[dimension];
            *(*(DRR+i)+j) = new float[dimension];
            *(*(RRR+i)+j) = new float[dimension];
        }
    }
    //inicializacion
    for (i = 0; i < dimension; i++)
    {
        for ( j = 0; j < dimension; j++)
        {
            for (k = 0; k < dimension; k++)
            {
                *(*(*(DDD + i) + j) + k)= 0.0;
                *(*(*(DDR + i) + j) + k)= 0.0;   
                *(*(*(DRR + i) + j) + k)= 0.0;
                *(*(*(RRR + i) + j) + k)= 0.0;
            }
        } 
    }
}

void eliminar_Histogramas3D(int dimension){
    int i,j;
    for (i = 0; i < dimension; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            delete[] *(*(DDD + i) + j);
            delete[] *(*(DDR + i) + j);
            delete[] *(*(DRR + i) + j);
            delete[] *(*(RRR + i) + j);
        }
        delete[] *(DDD + i);
        delete[] *(DDR + i);
        delete[] *(DRR + i);
        delete[] *(RRR + i);
    }
    delete[] DDD;
    delete[] DDR;
    delete[] DRR;
    delete[] RRR;
}

void eliminar_Datos(){
    delete[] DATA;
    delete[] RAND;
}
