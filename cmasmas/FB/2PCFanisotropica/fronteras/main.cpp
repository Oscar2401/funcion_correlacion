#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include "NODE.h"
#include <omp.h>

using namespace std;

void open_files(string, int, Point3D *);
void save_histogram(string, int, unsigned int *);

Point3D *dataD, *dataR;
unsigned int  *DD, *RR, *DR;
Node ***nodeD;
Node ***nodeR;

int main(int argc, char **argv){
	//int n_pts = stoi(argv[3]), bn = stoi(argv[4]);
	//float d_max = stof(argv[5]);
	int n_pts = 32768, bn = 10;
	float d_max = 100, size_box = 250, size_node = 14;
	dataD = new Point3D[n_pts]; // Asignamos meoria a esta variable
	dataR = new Point3D[n_pts];
	
	//Mensaje a usuario
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "Construcción de Histogramas DD, RR y DR para calcular" << endl;
	cout << "la función de correlación de 2 puntos isotrópica" << endl;
	cout << "implementando el método de mallas con condiciones" << endl;
	cout << "periódicas de frontera" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "Parametros usados: \n" << endl;
	cout << "Cantidad de puntos: " << n_pts << endl;
	cout << "Bins de histogramas: " << bn << endl;
	cout << "Distancia máxima: " << d_max << endl;
	cout << "Tamaño de nodos: " << size_node << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" << endl;
	// Nombre de los archivos 
	string nameDD = "DDaiso_mesh_3D_", nameRR = "RRiso_mesh_3D_", nameDR = "DRiso_mesh_3D_";
	nameDD.append(argv[3]);
	nameRR.append(argv[3]);
	nameDR.append(argv[3]);
	nameDD += ".dat";
	nameRR += ".dat";
	nameDR += ".dat";
	
	// inicializamos los histogramas
	DD = new unsigned int*[bn];
	RR = new unsigned int*[bn];
	DR = new unsigned int*[bn];
	int i,j;
	for (i = 0; i < bn; i++){
		*(DD+i) = new float[bn]; // vector[i]
		*(RR+i) = new float[bn];
		*(DR+i) = new float[bn];
	}
	for (i = 0; i < bn; i++){
		for ( j = 0; j < bn; j++){
			*(*(DD + i) + j) = 0.0;
			*(*(DR + i) + j) = 0.0;   
			*(*(RR + i) + j) = 0.0;
		} 
	}

	// Abrimos y trabajamos los datos en los histogramas
	open_files(argv[1],n_pts,dataD);
	open_files(argv[2],n_pts,dataR); // guardo los datos en los Struct
	
	// inicializamos las mallas
	int partitions = (int)(ceil(size_box/size_node));
	nodeD = new Node**[partitions];
	nodeR = new Node**[partitions];
	for ( i = 0; i < partitions; i++){
		*(nodeD + i) = new Node*[partitions];
		*(nodeR + i) = new Node*[partitions];
		for (int j = 0; j < partitions; j++){
			*(*(nodeD + i)+j) = new Node[partitions];
			*(*(nodeR + i)+j) = new Node[partitions];
		}
	}	
	
	// Iniciamos clase
	NODE my_hist(bn, n_pts, size_box, size_node, d_max, dataD, dataR, nodeD, nodeR);
	
	clock_t c_start = clock();
	
	my_hist.make_histoXX(DD, my_hist.meshData()); //hace histogramas XX
	
	clock_t c_end = clock();
	float time_elapsed_s = ((float)(c_end-c_start))/CLOCKS_PER_SEC;
	//my_hist.make_histoXX(RR, my_hist.meshRand());
	//my_hist.make_histoXY(DR, my_hist.meshData(), my_hist.meshRand()); //hace historamas XY
	my_hist.~NODE(); //destruimos objeto
	
	
	cout << "Termine de hacer todos los histogramas\n" << endl;
	
	// Mostramos los histogramas 
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" << endl;
	cout << "HITOGRAMA DD:" << endl;
	
	for (i = 0; i<bn; i++){
		printf("%d \t",DD[i]);
	}
	cout << "\n::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" << endl;
	cout << "HITOGRAMA RR:" << endl;
	for (i = 0; i<bn; i++){
		printf("%d \t",RR[i]);
	}
	cout << "\n::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" << endl;
	cout << "HITOGRAMA DR:" << endl;
	for (i = 0; i<bn; i++){
		printf("%d \t",DR[i]);
	}
	
	cout << "\n::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" << endl;
	save_histogram(nameDD, bn, DD);
	cout << "\nGuarde histograma DD..." << endl;
	save_histogram(nameRR, bn, RR);
	cout << "\nGuarde histograma RR..." << endl;
	save_histogram(nameDR, bn, DR);
	cout << "\nGuarde histograma DR..." << endl;
	
	// Eliminamos los hitogramas 
	//delete[] DD;
	//delete[] DR;
	//delete[] RR;
	
	printf("\nTiempo en CPU usado = %.4f seg.\n", time_elapsed_s );
	//printf("\nTiempo implementado = %.4f seg.\n", ((float))/CLOCKS_PER_SEC);
	cout << "Programa finalizado..." << endl;
	cin.get();
	return 0;
}

//====================================================================
//============ Sección de Funciones ================================== 
//====================================================================

void open_files(string name_file, int pts, Point3D *datos){
	/* Función para abrir nuestros archivos de datos */
	ifstream file;
	file.open(name_file.c_str(), ios::in | ios::binary); //le indico al programa que se trata de un archivo binario con ios::binary
	if (file.fail()){
		cout << "Error al cargar el archivo " << endl;
		exit(1);
	}
	
	//int c=0,remove;
	float remove;
	//while (!file.eof())
	for ( int c = 0; c < pts; c++)
	{
		file >> datos[c].x >> datos[c].y >> datos[c].z >> remove; 
		//c++;
	}
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
	for (int i = 0; i < bns; i++){
		file2 << histo[i] << endl;
	}
	file2.close();
}
