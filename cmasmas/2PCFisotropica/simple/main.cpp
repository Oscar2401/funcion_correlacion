#include <iostream>
#include <fstream> //manejo de archivos
#include <string.h>
#include "iso2histo.h"
#include <chrono>

using namespace std;

void open_files(string, int, Point3D *);
void save_histogram(string, int, unsigned int*);

// Variable globales
unsigned int *DD, *RR, *DR;
Point3D *dataD;
Point3D *dataR;

int main(int argc, char **argv){
	
	//int n_pts = stoi(argv[3]), bn = stoi(argv[4]);
	//float d_max = stof(argv[5]);
	int np = 32768, bn = 10;
	float dmax = 180.0;
	dataD = new Point3D[np]; // Asignamos meoria a esta variable
	dataR = new Point3D[np];
	
	// Nombre de los archivos 
	string nameDD = "DDiso_", nameRR = "RRiso_", nameDR = "DRiso_";
	nameDD.append(argv[3]);
	nameRR.append(argv[3]);
	nameDR.append(argv[3]);
	nameDD += ".dat";
	nameRR += ".dat";
	nameDR += ".dat";
	
	// Creamos los histogramas
	DD = new unsigned int[bn];
	RR = new unsigned int[bn];
	DR = new unsigned int[bn];
	for (int i = 0; i < bn; i++){
		*(DD+i) = 0.0; // vector[i]
		*(RR+i) = 0.0;
		*(DR+i) = 0.0;
	}
	
	// Abrimos y trabajamos los datos en los histogramas
	open_files(argv[1], np, dataD);
	open_files(argv[2], np, dataR); // guardo los datos en los Struct
	iso2hist my_hist(bn, np, dmax, dataD, dataR);
	
	auto start = std::chrono::system_clock::now();
	my_hist.make_histoXX(DD,RR); //hace histogramas XX
	my_hist.make_histoXY(DR); //hace historamas XY
	my_hist.~iso2hist(); //destruimos objeto
	
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>((end - start)); //mostramos los segundos que corre el programa
	printf("Time = %lld ms\n", static_cast<long long int>(elapsed.count()));
	
	// Eliminamos datos 
	delete[] dataD;
    	delete[] dataR;
	
	cout << "Termine de hacer todos los histogramas" << endl;
	// Mostramos los histogramas 
	cout << "\nHistograma DD:" << endl;
	for (int k = 0; k<bn; k++){
		printf("%d \t",DD[k]);
	}
	cout << "\nHistograma RR:" << endl;
	for (int k = 0; k<bn; k++){
		printf("%d \t",RR[k]);
	}
	cout << "\nHistograma DR:" << endl;
	for (int k = 0; k<bn; k++){
		printf("%d \t",DR[k]);
	}
	
	
	// Guardamos los histogramas
	save_histogram(nameDD, bn, DD);
	cout << "Guarde histograma DD..." << endl;
	save_histogram(nameRR, bn, RR);
	cout << "Guarde histograma RR..." << endl;
	save_histogram(nameDR, bn, DR);
	cout << "Guarde histograma DR..." << endl;
	
	//delete[] DD;
	//delete[] DR;
	//delete[] RR;
	cout << "Programa Terminado..." << endl;
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
	int remove;
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

