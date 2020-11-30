#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include "2PCF.h"
#include <omp.h>
#include <cmath>

using namespace std;

void open_files(string, int, PointW3D *);
void save_histogram(string, int, double *);
void save_histogram_analitic(string, int, double *);

PointW3D *dataD;
PointW3D *dataR;
double *DD; 
double *RR;
Node ***nodeD;
Node ***nodeR;

int main(int argc, char **argv){

	int n_pts = 32*32*32, bn = 400;
	float d_max = 180.0, size_box = 250.0, alpha = 2.176;
	float size_node = alpha*(size_box/pow((float)(n_pts),1/3.));
	dataD = new PointW3D[n_pts]; 
	dataR = new PointW3D[n_pts]; 
	
	cout << "\n        ISOTROPIC 2-POINT CORRELATION FUNCTION        \n" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "Construction of Histograms DD, RR to calculate" << endl;
	cout << "the isotropic 2-point correlation function" << endl;
	cout << "implementing the grid method and periodic boundary " << endl;
	cout << "conditions (BPC)." << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "Parameters used: \n" << endl;
	cout << "	Amount of points:     " << n_pts << endl;
	cout << "	Histogram Bins:       " << bn << endl;
	cout << "	Maximum distance:     " << d_max << endl;
	cout << "	Node size:            " << size_node << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	
	// File names
	string nameDD = "DDiso_mesh_3D_", nameRR = "RRiso_mesh_3D_";
	nameDD.append(argv[3]);
	nameRR.append(argv[3]);
	nameDD += ".dat";
	nameRR += ".dat";
	
	// Initialize the histograms
	DD = new double[bn];
	RR = new double[bn];
	int i, j;
	for (i = 0; i < bn; ++i){
		*(DD+i) = 0.0;
		*(RR+i) = 0.0;
	}

	open_files(argv[1],n_pts,dataD);
	open_files(argv[2],n_pts,dataR);
	
	// Initialize the grid
	int partitions = (int)(ceil(size_box/size_node));
	nodeD = new Node**[partitions];
	nodeR = new Node**[partitions];
	for ( i = 0; i < partitions; ++i){
		*(nodeD + i) = new Node*[partitions];
		*(nodeR + i) = new Node*[partitions];
		for (int j = 0; j < partitions; ++j){
			*(*(nodeD + i)+j) = new Node[partitions];
			*(*(nodeR + i)+j) = new Node[partitions];
		}
	}	
	
	// Start class
	NODE2P my_hist(bn, n_pts, size_box, size_node, d_max, dataD, nodeD, dataR, nodeR);
	delete[] dataD;
	delete[] dataR;
	
	clock_t c_start = clock();
	
	//construct histograms
	//===============================================
	my_hist.make_histoXX(DD, my_hist.meshData()); 
	save_histogram(nameDD, bn, DD);
	cout << "Save histogram DD ..." << endl;
	delete[] DD;
	//===============================================
	my_hist.make_histoXX(RR, my_hist.meshRand()); 
	save_histogram_analitic(nameRR, bn, RR);
	cout << "Save histogram RR ..." << endl;
	delete[] RR;
	//===============================================
	clock_t c_end = clock();
	float time_elapsed_s = ((float)(c_end-c_start))/CLOCKS_PER_SEC;
	
	my_hist.~NODE2P();
	
	
	cout << "Finish making all histograms" << endl;
	printf("\nCPU time used = %.4f seg.\n", time_elapsed_s );
	cout << "Program completed SUCCESSFULLY!" << endl;
	cin.get();
	return 0;
}

//====================================================================
//===================== Functions Section ============================
//====================================================================
void open_files(string name_file, int pts, PointW3D *datos){
	/* 
	Function to open our data files 
	*/
	ifstream file;
	file.open(name_file.c_str(), ios::in | ios::binary); 
	if (file.fail()){
		cout << "Error loading file! " << endl;
		exit(1);
	}
	for (int c = 0; c < pts; ++c) file >> datos[c].x >> datos[c].y >> datos[c].z >> datos[c].w; 
	file.close();
}

//====================================================================
void save_histogram(string name, int bns, double *histo){
	/* Función para guardar nuestros archivos de histogramas */
	ofstream file2;
	file2.open(name.c_str(), ios::out | ios::binary);
	
	if (file2.fail()){
		cout << "Failed to save file! " << endl;
		exit(1);
	}
	for (int i=0; i<bns; ++i) file2 << histo[i] << endl;
	file2.close();
}
//====================================================================
void save_histogram_analitic(string name, int bns, double *histo){
	/* Function to save our histogram files */
	ofstream file2;
	file2.open(name.c_str(), ios::out | ios::binary);
	
	if (file2.fail()){
		cout << "Failed to save file! " << endl;
		exit(1);
	}
	for (int i=0; i<bns; ++i) file2 << histo[i] << endl;
	file2.close();
}
