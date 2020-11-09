#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include "2PCFani.h"
#include <omp.h>
#include <cmath>

using namespace std;

void open_files(string, int, PointW3D *);
void save_histogram(string, int, float **);

PointW3D *dataD;
PointW3D *dataR;
float **DD; 
float **RR;
float **DR;
Node ***nodeD;
Node ***nodeR;

int main(int argc, char **argv){
	//int n_pts = stoi(argv[3]), bn = stoi(argv[4]);
	//float d_max = stof(argv[5]);
	//int n_pts = 32768, bn = 10;
	int n_pts = 32*32*32, bn = 50;
	float d_max = 150.0, size_box = 250.0, alpha = 2.176;
	float size_node = alpha*(size_box/pow((float)(n_pts),1/3.));
	dataD = new PointW3D[n_pts]; 
	dataR = new PointW3D[n_pts];
	
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
	string nameDD = "DDani_mesh_3D_", nameRR = "RRani_mesh_3D_", nameDR = "DRani_mesh_3D_";
	nameDD.append(argv[3]);
	nameRR.append(argv[3]);
	nameDR.append(argv[3]);
	nameDD += ".dat";
	nameRR += ".dat";
	nameDR += ".dat";
	
	// inicializamos los histogramas
	DD = new float*[bn];
	RR = new float*[bn];
	DR = new float*[bn];
	int i, j;
	for (i=0; i<bn; ++i){
		*(DD+i) = new float[bn];
		*(RR+i) = new float[bn];
		*(DR+i) = new float[bn];
	}
	for (i=0; i<bn; ++i){
	for (j=0; j<bn; ++j){
		*(*(DD+i)+j) = 0.0;
		*(*(RR+i)+j) = 0.0;
		*(*(DR+i)+j) = 0.0;
	}
	}

	// Abrimos y trabajamos los datos en los histogramas
	open_files(argv[1],n_pts,dataD);
	open_files(argv[2],n_pts,dataR);
	
	// inicializamos las mallas
	int partitions = (int)(ceil(size_box/size_node));
	nodeD = new Node**[partitions];
	nodeR = new Node**[partitions];
	for (i=0; i<partitions; ++i){
	*(nodeD+i) = new Node*[partitions];
	*(nodeR+i) = new Node*[partitions];
		for (j=0; j<partitions; ++j){
		*(*(nodeD+i)+j) = new Node[partitions];
		*(*(nodeR+i)+j) = new Node[partitions];
		}
	}	
	
	// Iniciamos clase
	NODE2P my_hist(bn, n_pts, size_box, size_node, d_max, dataD, nodeD, dataR, nodeR);
	delete[] dataD;
	delete[] dataR;
	clock_t c_start = clock();
	
	my_hist.make_histoXX(DD, my_hist.meshData()); 
	save_histogram(nameDD, bn, DD);
	cout << "Guarde histograma DD..." << endl;
	for (i=0; i<bn; ++i) delete[] *(DD+i);
	
	my_hist.make_histoXX(RR, my_hist.meshRand()); 
	save_histogram(nameRR, bn, RR);
	cout << "Guarde histograma RR..." << endl;
	for (i=0; i<bn; ++i) delete[] *(RR+i);
	
	my_hist.make_histoXY(DR, my_hist.meshData(), my_hist.meshRand()); 
	save_histogram(nameDR, bn, DR);
	cout << "Guarde histograma DR..." << endl;
	for (i=0; i<bn; ++i) delete[] *(DR+i);
	
	clock_t c_end = clock();
	float time_elapsed_s = ((float)(c_end-c_start))/CLOCKS_PER_SEC;
	
	my_hist.~NODE2P(); //destruimos objeto
	
	
	cout << "Termine de hacer todos los histogramas" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	
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
void save_histogram(string name, int bns, float **histo){
	int i, j;
	ofstream file;
	file.open(name.c_str(),ios::out | ios::binary);
	if (file.fail()){
		cout << "Error al guardar el archivo " << endl;
		exit(1);
	}
	for (i=0; i<bns; ++i){
		for (j=0; j<bns; ++j) file << histo[i][j] << " ";
		file << "\n";
	}
	file.close();
}
