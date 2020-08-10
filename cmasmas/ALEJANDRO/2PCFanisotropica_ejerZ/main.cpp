#include<iostream>
#include <fstream>
#include<string.h>
#include <stdlib.h> //memoria dinámia
#include <ctime> // manejo de tiempos de ejecución
#include "ani2PCF.h" // mi archivo de cabecera

using namespace std;

// Prototipos de funciones
void abrir_archivo(string,int, Puntos *);
void guardar_Histogramas(string, int, float **);
void crear_Histogramas2D(int); // un solo parametro, pues son matrices cuadradas
void eliminar_Histogramas2D(int);
void eliminar_Datos();


//variables globales
float **DD, **DR, **RR;
Puntos *DATA, *RAND;


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
    string nombre1 = "DDani_", nombre2 = "DRani_", nombre3 = "RRiso_";
    nombre1.append(argv[3]);
    nombre2.append(argv[3]);
    nombre3.append(argv[3]);
    nombre1 += ".dat";
    nombre2 += ".dat";
    nombre3 += ".dat";
    //cargo los datos
    abrir_archivo(argv[1],N,DATA);
    abrir_archivo(argv[2],N,RAND);
    //Inicializo con 0 los histogramas
    crear_Histogramas2D(nb);
    //programa principal
    ani2PCF obj(DATA,RAND,N,nb,d_max); // instancio la clase
    to = time(NULL);
    obj.calcular_histogramas_puros(DD,RR);
    obj.calcular_histogramas_mixtos(DR);
    tf = time(NULL);
    obj.~ani2PCF(); // destruyo el objeto
    eliminar_Datos(); // destruyo structs
    guardar_Histogramas(nombre1,nb,DD);
    guardar_Histogramas(nombre2,nb,DR);
    guardar_Histogramas(nombre3,nb,RR);
    eliminar_Histogramas2D(nb);
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

void guardar_Histogramas(string nombre, int dimension, float **Matriz){
    int i,j;
    ofstream archivo;
    archivo.open(nombre.c_str(),ios::out | ios::binary);
    if (archivo.fail()){
        cout << "Error al guardar el archivo " << endl;
        exit(1);
    }
    for (i = 0; i < dimension; i++)
    {
        for ( j = 0; j < dimension; j++)
        {
            archivo << *(*(Matriz + i) + j) << " "; 
        }
        archivo << endl;
    }
    archivo.close();
}

//Se crean matrices cuadradas y se inicializan con 0.0 
void crear_Histogramas2D(int dimension){
    int i,j;
    DD = new float*[dimension];
    DR = new float*[dimension];
    RR = new float*[dimension];
    for (i = 0; i < dimension; i++)
    {
        *(DD + i) = new float[dimension];
        *(DR + i) = new float[dimension];
        *(RR + i) = new float[dimension];
    }
    for (i = 0; i < dimension; i++)
    {
        for ( j = 0; j < dimension; j++)
        {
            *(*(DD + i) + j) = 0.0;
            *(*(DR + i) + j) = 0.0;   
            *(*(RR + i) + j) = 0.0;
        } 
    }
}

void eliminar_Histogramas2D(int dimension){
    int i;
    for (i = 0; i < dimension; i++)
    {
        delete[] *(DD + i);
        delete[] *(DR + i);
        delete[] *(RR + i); 
    }
    delete[] DD;
    delete[] DR;
    delete[] RR;
}

void eliminar_Datos(){
    delete[] DATA;
    delete[] RAND;
}
