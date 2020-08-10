#include <iostream>
#include <fstream> //manejo de archivos
#include <string.h>
#include "iso2histo.h"
#include <chrono>

using namespace std;

void open_files(string,int, Point3D *);
void save_histogram(string,int, float*);

// Variable globales
float *DD, *RR, *DR;
Point3D *dataD;
Point3D *dataR;

int main(int argc, char **argv){
	
	//int n_pts = stoi(argv[3]), bn = stoi(argv[4]);
	//float d_max = stof(argv[5]);
	int n_pts = 32768, bn = 30;
	float d_max = 180;
	dataD = new Point3D[n_pts]; // Asignamos meoria a esta variable
	dataR = new Point3D[n_pts];
	
	// Nombre de los archivos 
	string nameDD = "DDiso_", nameRR = "RRiso_", nameDR = "DRiso_";
	nameDD.append(argv[3]);
	nameRR.append(argv[3]);
	nameDR.append(argv[3]);
	nameDD += ".dat";
	nameRR += ".dat";
	nameDR += ".dat";
	
	// Creamos los histogramas
	DD = new float[bn];
	RR = new float[bn];
	DR = new float[bn];
	for (int i = 0; i<bn; i++){
		*(DD+i) = 0.0; // vector[i]
		*(RR+i) = 0.0;
		*(DR+i) = 0.0;
	}
	
	// Abrimos y trabajamos los datos en los histogramas
	open_files(argv[1],n_pts,dataD);
	open_files(argv[2],n_pts,dataR); // guardo los datos en los Struct
	iso2hist my_hist(bn,n_pts,d_max,dataD,dataR);
	
	auto start = std::chrono::system_clock::now();
	my_hist.make_histoXX(DD,RR); //hace histogramas XX
	my_hist.make_histoXY(DR); //hace historamas XY
	my_hist.~iso2hist(); //destruimos objeto
	
	// Guardamos los histogramas
	save_histogram(nameDD, bn, DD);
	save_histogram(nameRR, bn, RR);
	save_histogram(nameDR, bn, DR);
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>((end - start)/1000); //mostramos los segundos que corre el programa
	printf("Time = %lld ms\n", static_cast<long long int>(elapsed.count()));
	// Eliminamos los hitogramas 
	delete[] DD;
	delete[] DR;
	delete[] RR;
	cout << "listo" << endl;
	cin.get();
	return 0;
	
}

//====================================================================
//============ Sección de Funciones ================================== 
//====================================================================

void open_files(string name_file, int n_pts, Point3D *datos){
	/* Función para abrir nuestros archivos de datos */
	ifstream file;
	file.open(name_file.c_str(), ios::in | ios::binary); //le indico al programa que se trata de un archivo binario con ios::binary
	if (file.fail()){
		cout << "Error al cargar el archivo " << endl;
		exit(1);
	}
	
	int c=0,remove;
	while (!file.eof())
	{
		file >> datos[c].x >> datos[c].y >> datos[c].z >> remove; 
		c++;
	}
	file.close();
}

//====================================================================


void save_histogram(string name, int bn, float *histo){
	/* Función para guardar nuestros archivos de histogramas */
	ofstream file;
	file.open(name.c_str(),ios::out | ios::binary);
	if (file.fail()){
		cout << "Error al guardar el archivo " << endl;
		exit(1);
	}
	for (int i = 0; i < bn; i++){
		file << histo[i] << endl;
	}
	file.close();
}

