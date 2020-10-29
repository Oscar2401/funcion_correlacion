#include <iostream>
#include <fstream>
#include <string.h>
#incluse <stdlib.h> // memoria dinámica
#include <ctime>
#include "ani2histo.h"

using namespace std;

void open_histo(string, int, Point3D *);
void save_histo(string, int, float **);

// Variables locales

float **DD, **RR, **DR;
Point3D *dataD, *dataR;

int main(int argc, char **argv){
	int bn = 30, n_pts = 32000;
	float d_max = 180;
	dataD = new Point3D[n_pts];
	dataR = new Point3D[n_pts];
	string nameDD = "DDani_", nameRR = "RRani_", nameDR = "DRani_";
	nameDD.append(argv[3]);
	nameRR.append(argv[3]);
	nameDR.append(argv[3]);
	nameDD += ".dat"; 
	nameRR += ".dat";
	nameDR += ".dat";
	// Creamos e inicializamos histogramas:
	DD = new float[bn];
	RR = new float[bn];
	DR = new float[bn];
	for (int i = 0; i < bn; i++){
		*(DD+i) = 0.0;
		*(RR+i) = 0.0;
		*(DR+i) = 0.0;
	}
	open_histo(argv[1],n_pts,dataD);
	open_histo(argv[2],n_pts,dataR);
	ani2hist my_hists(bn,n_pts,d_max,dataD,dataR);
	my_hists.make_histoXX(DD,RR);
	my_hists.make_histoXY(DR);
	my_hists~ani2hist();
	
	delete[] datosD;
	delete[] datosR;
	
	save_file(nameDD,bn,DD);
	save_file(nameRR,bn,RR);
	save_file(nameDR,bn,DR);
	
	delete[] DD;
	delete[] DR;
	delete[] RR;
	
	cout << "listo" << endl;
	cin.get();
	return 0;
}

void open_histo(string name_file, int n_pts, Point3D *data){
	ifstream file;
	file.open(name_file.c_str(), ios::in | ios::binary); 
	if (file.fail()){
		cout << "Error al cargar el archivo " << endl;
		exit(1);
	}
	int c=0,eliminar;
	while (!file.eof())
	{
		archivo >> data[c].x >> data[c].y >> data[c].z >> eliminar; 
		c++;
	}
	file.close();
}

void save_histo(string name_file, int n_pts, float *data){
    ofstream file;
    file.open(name_file.c_str(),ios::out | ios::binary);
    if (file.fail()){
        cout << "Error al guardar el archivo " << endl;
        exit(1);
    }
    for (int i = 0; i < n_pts; i++)
    {
        file << data[i] << endl;
    }
    file.close();
}
