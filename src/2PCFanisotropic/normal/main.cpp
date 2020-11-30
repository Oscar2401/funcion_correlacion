#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include "2PCFani.h"
#include <omp.h>
#include <cmath>

using namespace std;

void open_files(string, int, PointW3D *);
void save_histogram(string, int, double **);

PointW3D *dataD;
PointW3D *dataR;
double **DD; 
double **RR;
double **DR;
Node ***nodeD;
Node ***nodeR;

int main(int argc, char **argv){

	int n_pts = 32*32*32, bn = 50;
	float d_max = 150.0, size_box = 250.0, alpha = 2.176;
	float size_node = alpha*(size_box/pow((float)(n_pts),1/3.));
	dataD = new PointW3D[n_pts]; 
	dataR = new PointW3D[n_pts];
	
	cout << "\n      ANISOTROPIC 2-POINT CORRELATION FUNCTION        \n" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "Construction of Histograms DD, RR and DR to calculate" << endl;
	cout << "the anisotropic 2-point correlation function" << endl;
	cout << "implementing the grid method." << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "Parameters used: \n" << endl;
	cout << "	Amount of points:     " << n_pts << endl;
	cout << "	Histogram Bins:       " << bn << endl;
	cout << "	Maximum distance:     " << d_max << endl;
	cout << "	Node size:            " << size_node << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	
	// File names
	string nameDD = "DDani_mesh_3D_", nameRR = "RRani_mesh_3D_", nameDR = "DRani_mesh_3D_";
	nameDD.append(argv[3]);
	nameRR.append(argv[3]);
	nameDR.append(argv[3]);
	nameDD += ".dat";
	nameRR += ".dat";
	nameDR += ".dat";
	
	// Initialize the histograms
	DD = new double*[bn];
	RR = new double*[bn];
	DR = new double*[bn];
	int i, j;
	for (i=0; i<bn; ++i){
		*(DD+i) = new double[bn];
		*(RR+i) = new double[bn];
		*(DR+i) = new double[bn];
	}
	for (i=0; i<bn; ++i){
	for (j=0; j<bn; ++j){
		*(*(DD+i)+j) = 0.0;
		*(*(RR+i)+j) = 0.0;
		*(*(DR+i)+j) = 0.0;
	}
	}
	
	open_files(argv[1],n_pts,dataD);
	open_files(argv[2],n_pts,dataR);
	
	// Initialize the grid
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
	
	// Start class
	NODE2P my_hist(bn, n_pts, size_box, size_node, d_max, dataD, nodeD, dataR, nodeR);
	delete[] dataD;
	delete[] dataR;
	clock_t c_start = clock();
	
	//construct histograms
	//==============================================
	my_hist.make_histoXX(DD, my_hist.meshData()); 
	save_histogram(nameDD, bn, DD);
	cout << "Save histogram DD ..." << endl;
	for (i=0; i<bn; ++i) delete[] *(DD+i);
	//==============================================
	my_hist.make_histoXX(RR, my_hist.meshRand()); 
	save_histogram(nameRR, bn, RR);
	cout << "Save histogram RR ..." << endl;
	for (i=0; i<bn; ++i) delete[] *(RR+i);
	//==============================================
	my_hist.make_histoXY(DR, my_hist.meshData(), my_hist.meshRand()); 
	save_histogram(nameDR, bn, DR);
	cout << "Save histogram DR ..." << endl;;
	for (i=0; i<bn; ++i) delete[] *(DR+i);
	//==============================================
	
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
//===================== Functions Section ===========================
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
void save_histogram(string name, int bns, double **histo){
	/* Function to save our histogram files */
	int i, j;
	ofstream file;
	file.open(name.c_str(),ios::out | ios::binary);
	if (file.fail()){
		cout << "Failed to save file! " << endl;
		exit(1);
	}
	for (i=0; i<bns; ++i){
		for (j=0; j<bns; ++j) file << histo[i][j] << " ";
		file << "\n";
	}
	file.close();
}
