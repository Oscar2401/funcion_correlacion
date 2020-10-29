#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include "2PCF.h"
#include <omp.h>
#include <cmath>

using namespace std;

void open_files(string, int, PointW3D *);
void save_histogram(string, int, unsigned int *);
void save_histogram_analitic(string, int, float *);

PointW3D *dataD;
unsigned int  *DD; 
float *RR;
Node ***nodeD;

int main(int argc, char **argv){
	//int n_pts = stoi(argv[3]), bn = stoi(argv[4]);
	//float d_max = stof(argv[5]);
	//int n_pts = 32768, bn = 10;
	int n_pts = 32*32*32, bn = 1000;
	float d_max = 60.0, size_box = 250.0, alpha = 2.176;
	float size_node = alpha*(size_box/pow((float)(n_pts),1/3.));
	dataD = new PointW3D[n_pts]; // Asignamos meoria a esta variable
	
	//Mensaje a usuario
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "Construcción de Histogramas DD, RR para calcular" << endl;
	cout << "la función de correlación de 2 puntos isotrópica" << endl;
	cout << "implementando el método de mallas con condiciones" << endl;
	cout << "periódicas de frontera" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "Parametros usados: \n" << endl;
	cout << "	Cantidad de puntos: " << n_pts << endl;
	cout << "	Bins de histogramas: " << bn << endl;
	cout << "	Distancia máxima: " << d_max << endl;
	cout << "	Tamaño de nodos: " << size_node << endl;
	cout << "\n::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	// Nombre de los archivos 
	string nameDD = "DDiso_mesh_3D_", nameRR = "RRiso_mesh_3D_", nameDR = "DRiso_mesh_3D_";
	nameDD.append(argv[2]);
	nameRR.append(argv[2]);
	nameDD += ".dat";
	nameRR += ".dat";
	nameDR += ".dat";
	
	// inicializamos los histogramas
	DD = new unsigned int[bn];
	RR = new float[bn];
	int i, j;
	for (i = 0; i < bn; ++i){
		*(DD+i) = 0; // vector[i]
		*(RR+i) = 0.0;
	}

	// Abrimos y trabajamos los datos en los histogramas
	open_files(argv[1],n_pts,dataD); // guardo los datos en los Struct
	
	// inicializamos las mallas
	int partitions = (int)(ceil(size_box/size_node));
	nodeD = new Node**[partitions];
	for ( i = 0; i < partitions; ++i){
		*(nodeD + i) = new Node*[partitions];
		for (int j = 0; j < partitions; ++j) *(*(nodeD + i)+j) = new Node[partitions];
	}	
	
	// Iniciamos clase
	NODE2P my_hist(bn, n_pts, size_box, size_node, d_max, dataD, nodeD);
	
	clock_t c_start = clock();
	
	my_hist.make_histoXX(DD, RR, my_hist.meshData()); //hace histogramas XX
	
	clock_t c_end = clock();
	float time_elapsed_s = ((float)(c_end-c_start))/CLOCKS_PER_SEC;
	
	my_hist.~NODE2P(); //destruimos objeto
	
	
	cout << "Termine de hacer todos los histogramas" << endl;

	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	save_histogram(nameDD, bn, DD);
	cout << "Guarde histograma DD..." << endl;
	save_histogram_analitic(nameRR, bn, RR);
	cout << "Guarde histograma RR..." << endl;

	// Eliminamos los hitogramas 
	delete[] DD;
	delete[] RR;
	
	printf("\nTiempo en CPU usado = %.4f seg.\n", time_elapsed_s );
	//printf("\nTiempo implementado = %.4f seg.\n", ((float))/CLOCKS_PER_SEC);
	cout << "Programa finalizado..." << endl;
	cin.get();
	return 0;
}

//====================================================================
//============ Sección de Funciones ================================== 
//====================================================================
void open_files(string name_file, int pts, PointW3D *datos){
	/* Función para abrir nuestros archivos de datos */
	ifstream file;
	file.open(name_file.c_str(), ios::in | ios::binary); //le indico al programa que se trata de un archivo binario con ios::binary
	if (file.fail()){
		cout << "Error al cargar el archivo " << endl;
		exit(1);
	}
	for (int c = 0; c < pts; ++c) file >> datos[c].x >> datos[c].y >> datos[c].z >> datos[c].w; 
	file.close();
}

//====================================================================
void save_histogram(string name, int bns, unsigned int *histo){
	/* Función para guardar nuestros archivos de histogramas */
	ofstream file2;
	file2.open(name.c_str(), ios::out | ios::binary);
	
	if (file2.fail()){
		cout << "Error al guardar el archivo " << endl;
		exit(1);
	}
	for (int i=0; i<bns; ++i) file2 << histo[i] << endl;
	file2.close();
}
//====================================================================
void save_histogram_analitic(string name, int bns, float *histo){
	/* Función para guardar nuestros archivos de histogramas */
	ofstream file2;
	file2.open(name.c_str(), ios::out | ios::binary);
	
	if (file2.fail()){
		cout << "Error al guardar el archivo " << endl;
		exit(1);
	}
	for (int i=0; i<bns; ++i) file2 << histo[i] << endl;
	file2.close();
}
