#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include "3PCFiso.h"
#include <omp.h>
#include <cmath>

using namespace std;

void open_files(string, int, PointW3D *, float size_box, float d_max);
void save_histogram(string, int, double ***);
void delete_histos(int);
void delete_dat();

PointW3D *dataD, *dataR;
double ***DDD;
double ***RRR;
double ***DDR;
double ***DRR;
Node ***nodeD;

int main(int argc, char **argv){

	int n_pts = 10000, bn = 40;
	float d_max = 60.0, size_box = 250.0, size_node =  2.17 * 250/pow(n_pts, (double)1/3);
	dataD = new PointW3D[n_pts]; 
	
	cout << "\n      ANISOTROPIC 3-POINT CORRELATION FUNCTION        \n" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "Construction of Histograms DDD, RRR, DDR and DRR to " << endl;
	cout << "calculate the isotropic 3-point correlation function" << endl;
	cout << "implementing the grid method and BPC." << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "Parameters used: \n" << endl;
	cout << "	Amount of points:     " << n_pts << endl;
	cout << "	Histogram Bins:       " << bn << endl;
	cout << "	Maximum distance:     " << d_max << endl;
	cout << "	Node size:            " << size_node << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	
	// File names
	string nameDDD = "DDDiso_mesh_3D_", nameRRR = "RRRiso_mesh_3D_", nameDDR = "DDRiso_mesh_3D_", nameDRR = "DRRiso_mesh_3D_";
	nameDDD.append(argv[2]);
	nameRRR.append(argv[2]);
	nameDDR.append(argv[2]);
	nameDRR.append(argv[2]);
	nameDDD += ".dat";
	nameRRR += ".dat";
	nameDDR += ".dat";
	nameDRR += ".dat";
	
	// Initialize the histograms
	DDD = new double**[bn];
	RRR = new double**[bn];
	DDR = new double**[bn];
	DRR = new double**[bn];
	int i,j,k;
	for (i=0; i<bn; i++){
	*(DDD+i) = new double*[bn];
	*(RRR+i) = new double*[bn];
	*(DDR+i) = new double*[bn];
	*(DRR+i) = new double*[bn];
		for (j = 0; j < bn; j++){
		*(*(DDD+i)+j) = new double[bn];
		*(*(RRR+i)+j) = new double[bn];
		*(*(DDR+i)+j) = new double[bn];
		*(*(DRR+i)+j) = new double[bn];
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
	
	open_files(argv[1],n_pts,dataD,size_box,d_max);
	
	// Initialize the grid
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
	delete_dat();
	clock_t c_start = clock();
	//====================================================================
	my_hist.make_histoXXX(DDD, my_hist.meshData());
	//====================================================================
	my_hist.make_histo_analitic(DDR, RRR, my_hist.meshData() );
	//====================================================================
	clock_t c_end = clock();
	float time_elapsed_s = ((float)(c_end-c_start))/CLOCKS_PER_SEC;
	
	my_hist.~NODE3P();
	
	cout << "\n::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" << endl;
	save_histogram(nameDDD, bn, DDD);
	cout << "Save histogram DDD ..." << endl;
	save_histogram(nameRRR, bn, RRR);
	cout << "Save histogram RRR ..." << endl;
	save_histogram(nameDDR, bn, DDR);
	cout << "Save histogram DDR ..." << endl;
	save_histogram(nameDRR, bn, RRR);
	cout << "Save histogram DRR ..." << endl;
	
	// Eliminamos los hitogramas
	delete_histos(bn);
	
	cout << "Finish making all histograms" << endl;
	printf("\nCPU time used = %.4f seg.\n", time_elapsed_s );
	cout << "Program completed SUCCESSFULLY!" << endl;
	cin.get();
	return 0;
}

//====================================================================
//===================== Functions Section ============================
//====================================================================
void open_files(string name_file, int pts, PointW3D *datos, float size_box, float d_max){
	/* 
	Function to open our data files 
	*/
	ifstream file;
	file.open(name_file.c_str(), ios::in | ios::binary);
	if (file.fail()){
		cout << "Error loading file! " << endl;
		exit(1);
	}
	int c;
	for (c=0; c<pts; ++c) file >> datos[c].x >> datos[c].y >> datos[c].z >> datos[c].w; 
	file.close();
	
	float front = size_box - d_max;
	// We tag border points
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
void save_histogram(string name, int bns, double ***histo){
	/* 
	Function to save our histogram files
	*/
	int i, j, k;
	ofstream file;
	file.open(name.c_str(),ios::out | ios::binary);
	if (file.fail()){
		cout << "Failed to save file! " << endl;
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
