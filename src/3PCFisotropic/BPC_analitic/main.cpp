#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include "3PCFiso.h"
#include <omp.h>
#include <cmath>

using namespace std;

void open_files(string, int, PointW3D *, float size_box, float d_max);
void save_histogram(string, int, float ***);
void delete_histos(int);
void delete_dat();

PointW3D *dataD, *dataR;
float  ***DDD, ***RRR, ***DDR, ***DRR;
Node ***nodeD;

int main(int argc, char **argv){
	//int n_pts = stoi(argv[3]), bn = stoi(argv[4]);
	//float d_max = stof(argv[5]);
	//int n_pts = 32768, bn = 10;
	int n_pts = 5000, bn = 20;
	float d_max = 60.0, size_box = 250.0, size_node =  2.17 * 250/pow(n_pts, (double)1/3);
	dataD = new PointW3D[n_pts]; 
	
	//Mensaje a usuario
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "Construcción de Histogramas DDD, RRR, DDR y DRR para " << endl;
	cout << "calcular la función de correlación de 3 puntos " << endl;
	cout << "isotrópica implementando el método de mallas con " << endl;
	cout << "condiciones periódicas de frontera" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "	Parametros usados: \n" << endl;
	cout << "	Cantidad de puntos: " << n_pts << endl;
	cout << "	Bins de histogramas: " << bn << endl;
	cout << "	Distancia máxima: " << d_max << endl;
	cout << "	Tamaño de nodos: " << size_node << endl;
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
	DDD = new float**[bn];
	RRR = new float**[bn];
	DDR = new float**[bn];
	DRR = new float**[bn];
	int i,j,k;
	for (i=0; i<bn; i++){
	*(DDD+i) = new float*[bn];
	*(RRR+i) = new float*[bn];
	*(DDR+i) = new float*[bn];
	*(DRR+i) = new float*[bn];
		for (j = 0; j < bn; j++){
		*(*(DDD+i)+j) = new float[bn];
		*(*(RRR+i)+j) = new float[bn];
		*(*(DDR+i)+j) = new float[bn];
		*(*(DRR+i)+j) = new float[bn];
		}
	}
	
	//inicialización
	for (i=0; i<bn; i++){
	for (j=0; j<bn; j++){
	for (k=0; k<bn; k++){
		*(*(*(DDD+i)+j)+k)= 0.0;
		*(*(*(DDR+i)+j)+k)= 0.0;   
		*(*(*(DRR+i)+j)+k)= 0.0;
		*(*(*(RRR+i)+j)+k)= 0.0;
	}
	} 
	}
	
	// Abrimos y trabajamos los datos en los histogramas
	open_files(argv[1],n_pts,dataD,size_box,d_max);
	
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
	
	my_hist.make_histoXXX(DDD, my_hist.meshData());
	my_hist.make_histo_analitic(DDR, RRR, my_hist.meshData() );
	
	clock_t c_end = clock();
	float time_elapsed_s = ((float)(c_end-c_start))/CLOCKS_PER_SEC;
	
	my_hist.~NODE3P(); //destruimos objeto
	
	//Eliminamos Datos 
	delete_dat();
	
	cout << "Termine de hacer todos los histogramas\n" << endl;
	cout << "\n::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" << endl;
	save_histogram(nameDDD, bn, DDD);
	cout << "\nGuarde histograma DDD..." << endl;
	save_histogram(nameRRR, bn, RRR);
	cout << "\nGuarde histograma RRR..." << endl;
	save_histogram(nameDDR, bn, DDR);
	cout << "\nGuarde histograma DDR..." << endl;
	save_histogram(nameDRR, bn, DDR);
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
void open_files(string name_file, int pts, PointW3D *datos, float size_box, float d_max){
	/* Función para abrir nuestros archivos de datos */
	ifstream file;
	file.open(name_file.c_str(), ios::in | ios::binary); //le indico al programa que se trata de un archivo binario con ios::binary
	if (file.fail()){
		cout << "Error al cargar el archivo " << endl;
		exit(1);
	}
	int c;
	for (c=0; c<pts; ++c) file >> datos[c].x >> datos[c].y >> datos[c].z >> datos[c].w; 
	file.close();
	
	float front = size_box - d_max;
	// Eticquetamos puntos frontera
	for (c=0; c<pts; ++c){
		// x:
		if (datos[c].x<d_max) datos[c].fx=1;
		else if(datos[c].x>front) datos[c].fx=-1;
		else datos[c].fx=0;
		// y:
		if (datos[c].y<d_max) datos[c].fy=1;
		else if(datos[c].y>front) datos[c].fy=-1;
		else datos[c].fy=0;
		// z:
		if (datos[c].z<d_max) datos[c].fz=1;
		else if(datos[c].z>front) datos[c].fz=-1;
		else datos[c].fz=0;
	}
}
//====================================================================
void save_histogram(string name, int bns, float ***histo){
	int i, j, k;
	ofstream file;
	file.open(name.c_str(),ios::out | ios::binary);
	if (file.fail()){
		cout << "Error al guardar el archivo " << endl;
		exit(1);
	}
	for (i=0; i<bns; i++){
		for (j=0; j<bns; j++){
			for (k=0; k<bns; k++){
				file << histo[i][j][k] << " "; 
			}
			file << "\n";
		}
		file << "\n" << endl;
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
