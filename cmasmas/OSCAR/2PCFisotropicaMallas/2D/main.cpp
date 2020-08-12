#include <iostream>
#include <fstream>
#include <string.h>
#include <chrono>
#include "NODE.h"

using namespace std;

void open_files(string, int, Point2D *);
void save_histogram(string, int, unsigned int *);

Point2D *dataD, *dataR;
unsigned int  *DD, *RR, *DR;
Node **nodeD;
Node **nodeR;

int main(int argc, char **argv){
	//int n_pts = stoi(argv[3]), bn = stoi(argv[4]);
	//float d_max = stof(argv[5]);
	int n_pts = 32768, bn = 30;
	float d_max = 180, size_box = 250, size_node = 10;
	dataD = new Point2D[n_pts]; // Asignamos meoria a esta variable
	dataR = new Point2D[n_pts];
	
	// Nombre de los archivos 
	string nameDD = "DDaiso_mesh_2D_", nameRR = "RRiso_mesh_2D_", nameDR = "DRiso_mesh_2D_";
	nameDD.append(argv[3]);
	nameRR.append(argv[3]);
	nameDR.append(argv[3]);
	nameDD += ".dat";
	nameRR += ".dat";
	nameDR += ".dat";
	
	// inicializamos los histogramas
	DD = new unsigned int[bn];
	RR = new unsigned int[bn];
	DR = new unsigned int[bn];
	int i;
	for (i = 0; i < bn; i++){
		*(DD+i) = 0.0; // vector[i]
		*(RR+i) = 0.0;
		*(DR+i) = 0.0;
	}

	// Abrimos y trabajamos los datos en los histogramas
	open_files(argv[1],n_pts,dataD);
	open_files(argv[2],n_pts,dataR); // guardo los datos en los Struct
	
	// inicializamos las mallas
	int partitions = (int)(ceil(size_box/size_node));
	nodeD = new Node*[partitions];
	nodeR = new Node*[partitions];
	for ( i = 0; i < partitions; i++){
		*(nodeD + i) = new Node[partitions];
		*(nodeR + i) = new Node[partitions];
	}
	
	// Iniciamos clase
	NODE my_hist(bn, n_pts, size_box, size_node, d_max, dataD, dataR, nodeD, nodeR);
	
	auto start = std::chrono::system_clock::now();
	
	my_hist.make_histoXX(DD, RR); //hace histogramas XX
	my_hist.make_histoXY(DR); //hace historamas XY
	my_hist.~NODE(); //destruimos objeto
	
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::seconds>((end - start)); //mostramos los segundos que corre el programa
	printf("Time = %lld s\n", static_cast<long long int>(elapsed.count()));
	cout << "Termine de hacer todos los histogramas" << endl;
	
	// Mostramos los histogramas 
	cout << "\nHistograma DD:" << endl;
	
	for (i = 0; i<bn; i++){
		printf("%d \t",DD[i]);
	}
	cout << "\nHistograma RR:" << endl;
	for (i = 0; i<bn; i++){
		printf("%d \t",RR[i]);
	}
	cout << "\nHistograma DR:" << endl;
	for (i = 0; i<bn; i++){
		printf("%d \t",DR[i]);
	}
	
	
	save_histogram(nameDD, bn, DD);
	cout << "\nGuarde histograma DD..." << endl;
	save_histogram(nameRR, bn, RR);
	cout << "Guarde histograma RR..." << endl;
	save_histogram(nameDR, bn, DR);
	cout << "Guarde histograma DR..." << endl;
	
	// Eliminamos los hitogramas 
	delete[] DD;
	delete[] DR;
	delete[] RR;
	cout << "Programa finalizado..." << endl;
	cin.get();
	return 0;
}

//====================================================================
//============ Sección de Funciones ================================== 
//====================================================================

void open_files(string name_file, int pts, Point2D *datos){
	/* Función para abrir nuestros archivos de datos */
	ifstream file;
	file.open(name_file.c_str(), ios::in | ios::binary); //le indico al programa que se trata de un archivo binario con ios::binary
	if (file.fail()){
		cout << "Error al cargar el archivo " << endl;
		exit(1);
	}
	
	//int c=0,remove;
	float remove1, remove2;
	//while (!file.eof())
	for ( int c = 0; c < pts; c++)
	{
		file >> datos[c].x >> datos[c].y >> remove1 >> remove2; 
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
