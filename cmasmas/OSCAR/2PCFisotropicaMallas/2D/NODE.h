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
		void make_nodos(Node **, Point2D *);
		void add(Point2D *&, int&, float, float);
	
	// Métodos de Clase:
	public:
		//Constructor de clase:
		NODE(int _bn, int _n_pts, float _size_box, float _size_node, float _d_max, Point2D *_dataD, Point2D *_dataR, Node **_nodeD, Node  **_nodeR){
			bn = _bn;
			n_pts = _n_pts;
			size_box = _size_box;
			size_node = _size_node;
			d_max = _d_max;
			dataD = _dataD;
			dataR = _dataR;
			nodeD = _nodeD;
			nodeR = _nodeR;
			std::cout << "\nHaciendo nodos..." << std::endl;
			make_nodos(nodeD,dataD);
			make_nodos(nodeR,dataR);

		}
		
		// Implementamos Método de mallas:
		void make_histoXX(unsigned int *, unsigned int*);
		void make_histoXY(unsigned int*);
		~NODE();
};

//=================================================================== 

void NODE::make_nodos(Node ** nod, Point2D *dat){
	/*
	Función para crear los nodos con los datos y puntos random
	*/
	int i, row, col, partitions = (int)(ceil(size_box/size_node));
	float p_med = size_node/2;
	
	// Inicializamos los nodos vacíos:
	for ( row = 0; row < partitions; row++){
		for ( col = 0; col < partitions; col++){
			nod[row][col].nodepos.y = ((float)(row)+p_med)*(size_node);
			nod[row][col].nodepos.x = ((float)(col)+p_med)*(size_node);
			nod[row][col].len = 0;
			nod[row][col].elements = new Point2D[0];
		}
	}
	// Llenamos los nodos con los puntos de data:
	for ( i = 0; i < n_pts; i++){
		col = (int)(dat[i].x/size_node);
        	row = (int)(dat[i].y/size_node);
		add( nod[row][col].elements, nod[row][col].len, dat[i].x, dat[i].y);
	}
}

//=================================================================== 

void NODE::add(Point2D *&array, int &lon, float _x, float _y){
	lon++;
	Point2D *array_aux = new Point2D[lon];
	for (int i = 0; i < lon-1; i++){
		array_aux[i].x = array[i].x;
		array_aux[i].y = array[i].y;
	}
	delete[] array;
	array = array_aux;
	array[lon-1].x = _x;
	array[lon-1].y = _y; 
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
					if (dis <= d_max){
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
					if (dis <= d_max){
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
			
			// Diatancias entre nodos del mismo renglon
			for (int m = 1; col + m < partitions ; m ++){
				u = row;
				v = col+m;
				
				x2D = nodeD[u][v].nodepos.x;
				y2D = nodeD[u][v].nodepos.y;
				dx_nod = x1D-x2D;
				dy_nod = y1D-y2D;
				dis_nod = sqrt(dx_nod*dx_nod + dy_nod*dy_nod)-corr;
				if (dis_nod <= d_max){
					for ( i = 0; i < nodeD[row][col].len; i++){
						for ( j = 0; j < nodeD[u][v].len; j++){	
							dx =  nodeD[row][col].elements[i].x-nodeD[u][v].elements[j].x;
							dy =  nodeD[row][col].elements[i].y-nodeD[u][v].elements[j].y;
							dis = sqrt(dx*dx + dy*dy);
							if (dis <= d_max){
								pos = (int)(dis*ds);
								DD[pos] += 2;
							}
						}
					}
				}
				
				x2R = nodeR[u][v].nodepos.x;
				y2R = nodeR[u][v].nodepos.y;
				dx_nod = x1R-x2R;
				dy_nod = y1R-y2R;
				dis_nod = sqrt(dx_nod*dx_nod + dy_nod*dy_nod)-corr;
				if (dis_nod <= d_max){
					for ( i = 0; i < nodeR[row][col].len; i++){
						for ( j = 0; j < nodeR[u][v].len; j++){	
							dx =  nodeR[row][col].elements[i].x-nodeR[u][v].elements[j].x;
							dy =  nodeR[row][col].elements[i].y-nodeR[u][v].elements[j].y;
							dis = sqrt(dx*dx + dy*dy);
							if (dis <= d_max){
								pos = (int)(dis*ds);
								RR[pos] += 2;
							}
						}
					}
				}
			}
			
			// Diatancias entre nodos de renglon diferente
			for ( u = row + 1; u < partitions; u++){
				for ( v = 0; v < partitions; v++){
					// Histograma DD
					x2D = nodeD[u][v].nodepos.x;
					y2D = nodeD[u][v].nodepos.y;
					dx_nod = x1D-x2D;
					dy_nod = y1D-y2D;
					dis_nod = sqrt(dx_nod*dx_nod + dy_nod*dy_nod)-corr;
					if (dis_nod <= d_max){
						for ( i = 0; i < nodeD[row][col].len; i++){
							for ( j = 0; j < nodeD[u][v].len; j++){	
								dx =  nodeD[row][col].elements[i].x-nodeD[u][v].elements[j].x;
								dy =  nodeD[row][col].elements[i].y-nodeD[u][v].elements[j].y;
								dis = sqrt(dx*dx + dy*dy);
								if (dis <= d_max){
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
					if (dis_nod <= d_max){
						for ( i = 0; i < nodeR[row][col].len; i++){
							for ( j = 0; j < nodeR[u][v].len; j++){	
								dx =  nodeR[row][col].elements[i].x-nodeR[u][v].elements[j].x;
								dy =  nodeR[row][col].elements[i].y-nodeR[u][v].elements[j].y;
								dis = sqrt(dx*dx + dy*dy);
								if (dis <= d_max){
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
					if (dis_nod <= d_max){
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
