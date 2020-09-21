#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include "NODE3PCF_frontV3.h"
#include <omp.h>
#include <cmath>

using namespace std;

void open_files(string, int, Point3D *);
void save_histogram(string, int, unsigned int ***);
void delete_histos(int);
void delete_dat();

Point3D *dataD, *dataR;
unsigned int  ***DDD, ***RRR, ***DDR, ***DRR;
Node ***nodeD;

int main(int argc, char **argv){
	//int n_pts = stoi(argv[3]), bn = stoi(argv[4]);
	//float d_max = stof(argv[5]);
	//int n_pts = 32768, bn = 10;
	int n_pts = 1000, bn = 10;
	float d_max = 100.0, size_box = 250.0, size_node = 2*250/10;
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
	string nameDDD = "DDDiso_mesh_3D_", nameRRR = "RRRiso_mesh_3D_", nameDDR = "DDRiso_mesh_3D_", nameDRR = "DRRiso_mesh_3D_";
	nameDDD.append(argv[2]);
	nameRRR.append(argv[2]);
	nameDDR.append(argv[2]);
	nameDRR.append(argv[2]);
	nameDDD += ".dat";
	nameRRR += ".dat";
	nameDDR += ".dat";
	nameDRR += ".dat";
	
	// inicializamos los histogramas
	DDD = new unsigned int**[bn];
	RRR = new unsigned int**[bn];
	DDR = new unsigned int**[bn];
	DRR = new unsigned int**[bn];
	int i,j,k;
	for (i=0; i<bn; i++){
		*(DDD+i) = new unsigned int*[bn];
		*(RRR+i) = new unsigned int*[bn];
		*(DDR+i) = new unsigned int*[bn];
		*(DRR+i) = new unsigned int*[bn];
		for (j = 0; j < bn; j++){
			*(*(DDD+i)+j) = new unsigned int[bn];
			*(*(RRR+i)+j) = new unsigned int[bn];
			*(*(DDR+i)+j) = new unsigned int[bn];
			*(*(DRR+i)+j) = new unsigned int[bn];
		}
	}
	
	//inicialización
	 for (i=0; i<bn; i++){
	 	for (j=0; j<bn; j++){
	 		for (k = 0; k < bn; k++){
	 			*(*(*(DDD+i)+j)+k)= 0;
	 			*(*(*(DDR+i)+j)+k)= 0;   
	 			*(*(*(DRR+i)+j)+k)= 0;
	 			*(*(*(RRR+i)+j)+k)= 0;
	 		}
	 	} 
	 }

	// Abrimos y trabajamos los datos en los histogramas
	open_files(argv[1],n_pts,dataD);
	
	// inicializamos las mallas
	int partitions = (int)(ceil(size_box/size_node));
	nodeD = new Node**[partitions];
	for (i=0; i<partitions; i++){
		*(nodeD+i) = new Node*[partitions];
		for (j=0; j<partitions; j++){
			*(*(nodeD+i)+j) = new Node[partitions];
		}
	}	
	
	// Iniciamos clase
	NODE3P my_hist(bn, n_pts, size_box, size_node, d_max, dataD, nodeD);
	clock_t c_start = clock();
	my_hist.make_histoXXX(DDD, my_hist.meshData()); //hace histogramas
	clock_t c_end = clock();
	float time_elapsed_s = ((float)(c_end-c_start))/CLOCKS_PER_SEC;
	my_hist.~NODE3P(); //destruimos objeto
	
	//Eliminamos Datos 
	delete_dat();
	
	cout << "Termine de hacer todos los histogramas\n" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" << endl;
	cout << "HITOGRAMA DDD:" << endl;
	
	k = 0;
	for (i=0; i<bn; i++){
		 printf("\n");
		for (j=0; j<bn; j++){
			string num = to_string(DDD[i][j][k]);
			cout << DDD[i][j][k] << std::string(7-num.size(), ' ');
		}
	}
	
	// Mostramos los histogramas 
	cout << "\n::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" << endl;
	save_histogram(nameDDD, bn, DDD);
	cout << "\nGuarde histograma DDD..." << endl;
	save_histogram(nameRRR, bn, RRR);
	cout << "\nGuarde histograma RRR..." << endl;
	save_histogram(nameDDR, bn, DDR);
	cout << "\nGuarde histograma DDR..." << endl;
	save_histogram(nameDRR, bn, DRR);
	cout << "\nGuarde histograma DRR..." << endl;
	
	// Eliminamos los hitogramas
	delete_histos(bn);
	
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
void save_histogram(string name, int bns, unsigned int ***histo){
	int i, j, k, d=0;
	unsigned int **reshape = new unsigned int*[bns];
	for (i=0; i<bns; i++){
		*(reshape+i) = new unsigned int[bns*bns];
        }
	for (i=0; i<bns; i++){
	for (j=0; j<bns; j++){
	for (k=0; k<bns; k++){
		reshape[i][bns*j+k] = histo[i][j][k];
	}
	}
	}
	ofstream file;
	file.open(name.c_str(),ios::out | ios::binary);
	if (file.fail()){
		cout << "Error al guardar el archivo " << endl;
		exit(1);
	}
	for (i=0; i<bns; i++){
		for (j=0; j<bns*bns; j++){
			file << reshape[i][j] << " "; 
		}
		file << endl;
	}
	file.close();
}
//====================================================================
void delete_histos(int dim){
	int i,j;
	for (i = 0; i < dim; i++){
	for (j = 0; j < dim; j++){
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
//====================================================================
void delete_dat(){
    delete[] dataD;
}
