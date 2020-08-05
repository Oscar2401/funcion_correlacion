
#include <iostream> // estandar de entrada y salida
#include <fstream> // manejo de archivos
#include <string.h> // manejo de nombres
#include "iso2PCF.h" // incluyo mi archivo de cabecera
//#include "distancias.h"

using namespace std;

// Prototipos de funciones
void abrir_archivo(string,int, Punto *);
void guardar_Histograma(string,int, float*);
void crear_Histogramas(int dim);
void eliminar_Histogramas();
void eliminar_Datos();

//Variables globales
float *DD, *DR, *RR;
Punto *datosD;
Punto *datosR;
Punto *mesh;
Punto *PM;

// primer argumento del main es el nombre del archivo datos
// segundo argumento del main es el nombre del archivo rand
// tercer argumento es la cantidad de datos a trabajar
// cuarto argumento es el numero de bins
// quinto argumento es la distancia maxima
int main(int argc, char **argv){
	int N = 500; 
	datosD = new Punto[N]; // creacion de N struct de puntos
	mesh = new Punto[N];
	PM = new Punto[N];
	
	abrir_archivo(argv[1],N,datosD);
	int tam = 10;
	make_mesh(datosD, mesh, PM, N, tam);
	
	return 0;
}

//toma los datos del archivo y los guarda en un arreglo de structuras.
void abrir_archivo(string nombre_archivo,int cantidad_puntos, Punto *datos){
    ifstream archivo;
    archivo.open(nombre_archivo.c_str(), ios::in | ios::binary); //le indico al programa que se trata de un archivo binario con ios::binary
    if (archivo.fail()){
        cout << "Error al cargar el archivo " << endl;
        exit(1);
    }
    int c=0,eliminar;
    while (!archivo.eof())
    {
        archivo >> datos[c].x >> datos[c].y >> datos[c].z >> eliminar;
        c++;
    }
    archivo.close();
}

void guardar_Histograma(string nombre,int dim, float*histograma){
    ofstream archivo;
    archivo.open(nombre.c_str(),ios::out | ios::binary);
    if (archivo.fail()){
        cout << "Error al guardar el archivo " << endl;
        exit(1);
    }
    for (int i = 0; i < dim; i++)
    {
        archivo << histograma[i] << endl;
    }
    archivo.close();
}


void crear_Histogramas(int dim){
    DD = new float[dim];
    DR = new float[dim];
    RR = new float[dim];
    for (int i = 0; i < dim; i++)
    {
        *(DD+i) = 0.0; // vector[i]
        *(DR+i) = 0.0;
        *(RR+i) = 0.0;
    }
}

void eliminar_Histogramas(){
    delete[] DD;
    delete[] DR;
    delete[] RR;
}

void eliminar_Datos(){
    delete[] datosD;
    delete[] datosR;
}
