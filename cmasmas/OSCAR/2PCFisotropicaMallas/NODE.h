#include <stdlib.h>
#include <stdio.h>
#include <cmath>

struct Point2D{
	float x;
	float y;	
};

struct Point3D{
	float x;
	float y;
	float z;
};

struct Node{
	Point2D nodepos;		// Coordenadas del nodo (posición del nodo).
	int len;				// Cantidad de elementos en el nodo.
	Point2D *elements;	// Elementos del nodo.
};

class NODE{
	//Atributos de clase:
	private:
		int bn;
		int n_pts;
		float size_box;
		float size_node;
		float d_max;
		Node **nodeD;
		Node **nodeR;
		Point2D *dataD;
		Point2D *dataR;

	private: 
		void make_nodos();
	
	// Métodos de Clase:
	public:
		//Constructor de clase:
		NODE(int _bn, int _n_pts, float _size_box, float _size_node, float _d_max, Point2D *_dataD, Point2D *_dataR){
			bn = _bn;
			n_pts = _n_pts;
			size_box = _size_box;
			size_node = _size_node;
			d_max = _d_max;
			dataD = _dataD;
			dataR = _dataR;
			std::cout << "\nHaciendo nodos..." << std::endl;
			make_nodos();

		}
		void setBins(int _bn){
			bn = _bn;
		}
		void setNpts(int _n_pts){
			n_pts = _n_pts;
		}
		void setSizeBox(float _size_box){
			size_box = _size_box;
		}
		void setSizeNode(float _size_node){
			size_node = _size_node;
		}
		void setDmax(float _d_max){
			d_max = _d_max;
		}
		void setDataD(Point2D *_dataD){
			dataD = _dataD;
		}
		void setDataR(Point2D *_dataR){
			dataR = _dataR;
		}

		// Implementamos Método de mallas:
		void make_histoXX(unsigned int *, unsigned int*);
		void make_histoXY(unsigned int*);
		~NODE();
};

//=================================================================== 

void NODE::make_nodos(){
	/*
	Función para crear los nodos con los datos y puntos random
	*/
	int i, j, rowD, colD, rowR, colR, lengthD, lengthR, partitions = (int)(ceil(size_box/size_node));
	float pmed = size_node/2;
	
	nodeD = new Node*[partitions];
	nodeR = new Node*[partitions];
	for ( i = 0; i < partitions; i++){
		*(nodeD+i) = new Node[partitions];
		*(nodeR+i) = new Node[partitions];
	}
	
	// Inicializamos los nodos vacíos:
	for (int row = 0; row < partitions; row++){
		for (int col = 0; col < partitions; col++){
			nodeD[row][col].nodepos.y = ((float)(row)*(size_node))+pmed;
			nodeD[row][col].nodepos.x = ((float)(col)*(size_node))+pmed;
			nodeD[row][col].len = 0;
			nodeD[row][col].elements = new Point2D[0];
			
			nodeR[row][col].nodepos.y = ((float)(row)*(size_node))+pmed;
			nodeR[row][col].nodepos.x = ((float)(col)*(size_node))+pmed;
			nodeR[row][col].len = 0;
			nodeR[row][col].elements = new Point2D[0];
		}
	}

	// Llenamos los nodos con los puntos de data:
	for ( i = 0; i < n_pts; i++){
		
		// Dirección de punto en nodos
		colD = (int)(dataD[i].x/size_node);
        	rowD = (int)(dataD[i].y/size_node);
		
		colR = (int)(dataR[i].x/size_node);
        	rowR = (int)(dataR[i].y/size_node);
		
		nodeD[rowD][colD].len++; // incrementamos la cantidad de puntos en el nodo.
		lengthD = nodeD[rowD][colD].len;
		
		nodeR[rowR][colR].len++; // incrementamos la cantidad de puntos en el nodo.
		lengthR = nodeR[rowR][colR].len;

		Point2D *array_auxD = new Point2D[lengthD];
		Point2D *array_auxR = new Point2D[lengthR];

		for ( j = 0; j < lengthD-1; j++){
			array_auxD[j].x = nodeD[rowD][colD].elements[j].x;
			array_auxD[j].y = nodeD[rowD][colD].elements[j].y;	
		}
		for ( j = 0; j < lengthR-1; j++){	
			array_auxR[j].x = nodeR[rowR][colR].elements[j].x;
			array_auxR[j].y = nodeR[rowR][colR].elements[j].y;	
		}

		delete[] nodeD[rowD][colD].elements;
		delete[] nodeR[rowR][colR].elements;
		nodeD[rowD][colD].elements = array_auxD;
		nodeR[rowR][colR].elements = array_auxR;
		
		// Guardamos punto en nodo
		nodeD[rowD][colD].elements[lengthD-1].x = dataD[i].x;
		nodeD[rowD][colD].elements[lengthD-1].y = dataD[i].y;
		
		nodeR[rowR][colR].elements[lengthR-1].x = dataR[i].x;
		nodeR[rowR][colR].elements[lengthR-1].y = dataR[i].y;
			
	}
}


//=================================================================== 

void NODE::make_histoXX(unsigned int *DD, unsigned int *RR){
	int i, j, u, v, row, col, pos, partitions = (int)(ceil(size_box/size_node));
	float x1D, x2D, y1D, y2D, x1R, x2R, y1R, y2R;
	float dx, dy, dx_nod, dy_nod, corr = size_node*sqrt(2);
	float dis, dis_nod;
	float ds = ((float)(bn))/d_max;
	std::cout << "Estoy haciendo histogramas DD y RR..." << std::endl;
	for (row = 0; row < partitions; row++){
		for (col = 0; col < partitions; col++){
			// Distancias entre puntos del mismo nodo:
			
			// Histograma DD
			for ( i= 0; i < nodeD[row][col].len - 1; i++){
				for ( j = i+1; j < nodeD[row][col].len; j++){	
					dx =  nodeD[row][col].elements[i].x-nodeD[row][col].elements[j].x;
					dy =  nodeD[row][col].elements[i].y-nodeD[row][col].elements[j].y;
					dis = sqrt(dx*dx + dy*dy);
					if (dis < d_max){
						pos = (int)(dis*ds);
						DD[pos] += 2;
					}
				}
			}

			// Histograma RR
			for ( i= 0; i < nodeR[row][col].len - 1; i++){
				for ( j = i+1; j < nodeR[row][col].len; j++){	
					dx =  nodeR[row][col].elements[i].x-nodeR[row][col].elements[j].x;
					dy =  nodeR[row][col].elements[i].y-nodeR[row][col].elements[j].y;
					dis = sqrt(dx*dx + dy*dy);
					if (dis < d_max){
						pos = (int)(dis*ds);
						RR[pos] += 2;
					}
				}
			}

			// Distancias entre puntos del diferente nodo:
			x1D = nodeD[row][col].nodepos.x;
			y1D = nodeD[row][col].nodepos.y;
			
			x1R = nodeR[row][col].nodepos.x;
			y1R = nodeR[row][col].nodepos.y;
			for ( u = row; u < partitions; u++){
				for ( v = col; v < partitions; v++){
					if (u !=row || v !=col){ 
						
						// Histograma DD
						x2D = nodeD[u][v].nodepos.x;
						y2D = nodeD[u][v].nodepos.y;
						dx_nod = x1D-x2D;
						dy_nod = y1D-y2D;
						dis_nod = sqrt(dx_nod*dx_nod + dy_nod*dy_nod)-corr;
						if (dis_nod < d_max){
							for ( i = 0; i < nodeD[row][col].len; i++){
								for ( j = 0; j < nodeD[u][v].len; j++){	
									dx =  nodeD[row][col].elements[i].x-nodeD[u][v].elements[j].x;
									dy =  nodeD[row][col].elements[i].y-nodeD[u][v].elements[j].y;
									dis = sqrt(dx*dx + dy*dy);
									if (dis < d_max){
										pos = (int)(dis*ds);
										DD[pos] += 2;
									}
								}
							}
						}
						
						// Histograma RR
						x2R = nodeR[u][v].nodepos.x;
						y2R = nodeR[u][v].nodepos.y;
						dx_nod = x1R-x2R;
						dy_nod = y1R-y2R;
						dis_nod = sqrt(dx_nod*dx_nod + dy_nod*dy_nod)-corr;
						if (dis_nod<d_max){
							for ( i = 0; i < nodeR[row][col].len; i++){
								for ( j = 0; j < nodeR[u][v].len; j++){	
									dx =  nodeR[row][col].elements[i].x-nodeR[u][v].elements[j].x;
									dy =  nodeR[row][col].elements[i].y-nodeR[u][v].elements[j].y;
									dis = sqrt(dx*dx + dy*dy);
									if (dis < d_max){
										pos = (int)(dis*ds);
										RR[pos] += 2;
									}
								}
							}
						}
					}
				}	
			}
		}
	}
}

//=================================================================== 

void NODE::make_histoXY(unsigned int *DR){
	int i, j, u, v, row, col, pos, partitions = (int)(ceil(size_box/size_node));
	float x1, x2, y1, y2;
	float dx, dy, dx_nod, dy_nod, corr = size_node*sqrt(2);
	float dis, dis_nod;
	float ds = ((float)(bn))/d_max;
	std::cout << "Estoy haciendo histograma DR..." << std::endl;
	for (row = 0; row < partitions; row++){
		for (col = 0; col < partitions; col++){
			// Distancias entre puntos del diferente nodo:
			
			x1 = nodeD[row][col].nodepos.x;
			y1 = nodeD[row][col].nodepos.y;
			for ( u = 0; u < partitions; u++){
				for ( v = 0; v < partitions; v++){
					// Histograma DR
					x2 = nodeR[u][v].nodepos.x;
					y2 = nodeR[u][v].nodepos.y;
					dx_nod = x1-x2;
					dy_nod = y1-y2;
					dis_nod = sqrt(dx_nod*dx_nod + dy_nod*dy_nod)-corr;
					if (dis_nod < d_max){
						for ( i = 0; i < nodeD[row][col].len; i++){
							for ( j = 0; j < nodeR[u][v].len; j++){	
								dx =  nodeD[row][col].elements[i].x-nodeR[u][v].elements[j].x;
								dy =  nodeD[row][col].elements[i].y-nodeR[u][v].elements[j].y;
								dis = sqrt(dx*dx + dy*dy);
								if (dis < d_max){
									pos = (int)(dis*ds);
									DR[pos] += 1;
								}
							}
						}
					}
				}	
			}
		}
	}
}

//=================================================================== 

NODE::~NODE(){
	
}
