
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <omp.h>

struct Point3D{
	float x;
	float y; 
	float z;
};

struct Node{
	Point3D nodepos;	// Coordenadas del nodo (posición del nodo).
	int len;		// Cantidad de elementos en el nodo.
	Point3D *elements;	// Elementos del nodo.
};

//=================================================================== 
//======================== Clase ==================================== 
//=================================================================== 

class NODE{
	//Atributos de clase:
	private:
		// Asignados
		int bn;
		int n_pts;
		float size_box;
		float size_node;
		float d_max;
		Node ***nodeD;
		Node ***nodeR;
		Point3D *dataD;
		Point3D *dataR;
		// Derivados
		float ll;
		float dd_max;
		float corr;
		float front;
		float ds;
		float ddmax_nod;
		
	private: 
		void make_nodos(Node ***, Point3D *);
		void add(Point3D *&, int&, float, float, float);
	
	// Métodos de Clase:
	public:
		//Constructor de clase:
		NODE(int _bn, int _n_pts, float _size_box, float _size_node, float _d_max, Point3D *_dataD, Point3D *_dataR, Node ***_nodeD, Node  ***_nodeR){
			bn = _bn;
			n_pts = _n_pts;
			size_box = _size_box;
			size_node = _size_node;
			d_max = _d_max;
			dataD = _dataD;
			dataR = _dataR;
			nodeD = _nodeD;
			nodeR = _nodeR;
			ll = size_box*size_box;
			dd_max = d_max*d_max;
			front = size_box - d_max;
			corr = size_node*sqrt(3);
			ds = ((float)(bn))/d_max;
			ddmax_nod = (d_max+corr)*(d_max+corr);
			make_nodos(nodeD,dataD); 
			make_nodos(nodeR,dataR);
			std::cout << "Terminé de construir nodos..." << std::endl;
		}
		
		Node ***meshData(){
			return nodeD;
		};
		Node ***meshRand(){
			return nodeR;
		};
		
		// Implementamos Método de mallas:
		void make_histoXXX(unsigned int ***, unsigned int ***, Node ***);
		void count_3_N111(int, int, int, unsigned int ***, Node ***);
		void count_3_N112(int, int, int, int, int, int, unsigned int ***, Node ***);
		void count_3_N122(int, int, int, int, int, int, unsigned int ***, Node ***);
		void count_3_N123(int, int, int, int, int, int, int, int, int, unsigned int ***, Node ***);
		//void histo_front_XX(unsigned int *, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int);
		//void histo_front_XY(unsigned int *, Node ***, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int);
		~NODE();
};

//=================================================================== 
//==================== Funciones ==================================== 
//=================================================================== 
//=================================================================== 

void NODE::make_nodos(Node ***nod, Point3D *dat){
	/*
	Función para crear los nodos con los datos y puntos random
	
	Argumentos
	nod: arreglo donde se crean los nodos.
	dat: datos a dividir en nodos.
	
	*/
	int i, row, col, mom, partitions = (int)((size_box/size_node)+1);
	float p_med = size_node/2;
	
	// Inicializamos los nodos vacíos:
	for (row=0; row<partitions; row++){
		for (col=0; col<partitions; col++){
			for (mom=0; mom<partitions; mom++){
				nod[row][col][mom].nodepos.z = ((float)(mom)*(size_node))+p_med;
				nod[row][col][mom].nodepos.y = ((float)(col)*(size_node))+p_med;
				nod[row][col][mom].nodepos.x = ((float)(row)*(size_node))+p_med;
				nod[row][col][mom].len = 0;
				nod[row][col][mom].elements = new Point3D[0];
			}
		}
	}
	// Llenamos los nodos con los puntos de dat:
	for (i=0; i<n_pts; i++){
		row = (int)(dat[i].x/size_node);
        	col = (int)(dat[i].y/size_node);
        	mom = (int)(dat[i].z/size_node);
		add( nod[row][col][mom].elements, nod[row][col][mom].len, dat[i].x, dat[i].y, dat[i].z);
	}
}

//=================================================================== 

void NODE::add(Point3D *&array, int &lon, float _x, float _y, float _z){
	lon++;
	Point3D *array_aux = new Point3D[lon];
	for (int i=0; i<lon-1; i++){
		array_aux[i].x = array[i].x;
		array_aux[i].y = array[i].y;
		array_aux[i].z = array[i].z;
	}
	delete[] array;
	array = array_aux;
	array[lon-1].x = _x;
	array[lon-1].y = _y; 
	array[lon-1].z = _z; 
}

//=================================================================== 

void NODE::make_histoXXX(unsigned int ***XXX, unsigned int ***YYY, Node ***nodeX){
	/*
	Función para crear los histogramas DD y RR.
	
	Argumentos
	DD: arreglo donde se creará el histograma DD.
	RR: arreglo donde se creará el histograma RR.
	
	*/
	//Variables compartidas en hilos: 
	int partitions = (int)((size_box/size_node)+1);
	int i, j, row, col, mom, u, v, w, a ,b, c;
	float dis, dis_nod;
	float x1N, y1N, z1N, x2N, y2N, z2N, x3N, y3N, z3N;
	float x, y, z;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod, dx_nod2, dy_nod2;
	bool con_x, con_y, con_z;
	float d_max_pm = d_max + size_node/2, front_pm = front - size_node/2;
	
	std::cout << "-> Estoy haciendo histograma XXX..." << std::endl;
	
	// x1N, y1N, z1N => Nodo pivote
	for (row = 0; row < partitions; ++row){
	x1N = nodeX[row][0][0].nodepos.x;
	for (col = 0; col < partitions; ++col){
	y1N = nodeX[row][col][0].nodepos.y;
	for (mom = 0; mom < partitions; ++mom){
	z1N = nodeX[row][col][mom].nodepos.z;	
			
	//==================================================
	// Triángulos entre puntos del mismo nodo:
	//==================================================
	count_3_N111(row, col, mom, XXX, nodeX);
	//==================================================
	// Triángulos entre puntos del diferente nodo:
	//==================================================
	u = row;
	v = col;
	//=======================
	// Nodo 2 movil en Z:
	//=======================
	for (w=mom+1;  w<partitions ; ++w){	
		z2N = nodeX[u][v][w].nodepos.z;
		dz_nod = z1N-z2N;
		dis_nod = dz_nod*dz_nod;
		if (dis_nod <= ddmax_nod){
		//==============================================
		// 2 puntos en N1 y 1 punto en N2
		//==============================================
		count_3_N112(row, col, mom, u, v, w, XXX, nodeX);
		//==============================================
		// 2 puntos en N1 y 1 punto en N2
		//==============================================
		count_3_N122(row, col, mom, u, v, w, XXX, nodeX);
		//==============================================
		// 1 punto en N1, 1 punto en N2 y un punto en N3
		//==============================================
		a = u;
		b = v;
		//=======================
		// Nodo 3 movil en Z:
		//=======================
		for (c=w+1;  c<partitions; ++c){
			z3N = nodeX[a][b][c].nodepos.z;
			dz_nod = z1N-z3N;
			dis_nod = dz_nod*dz_nod;
			if (dis_nod <= ddmax_nod){
			dz_nod = z2N-z3N;
			dis_nod = dz_nod*dz_nod;
			if (dis_nod <= ddmax_nod){
			count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);///////////////////////////////////
			}
			}
		}
		//=======================
		// Nodo 3 movil en ZY:
		//=======================
		for (b=v+1; b<partitions; ++b){
			y3N = nodeX[a][b][0].nodepos.y;
			dy_nod = y1N-y3N;
			for (c=0;  c<partitions; ++c){
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod = z1N-z3N;
				dis_nod = dy_nod*dy_nod + dz_nod*dz_nod;
				if (dis_nod <= ddmax_nod){
				dz_nod = z2N-z3N;
				dis_nod = dy_nod*dy_nod + dz_nod*dz_nod;
				if (dis_nod <= ddmax_nod){
				count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
				}
				}
		}
		}
		//=======================
		// Nodo 3 movil en ZYX:
		//=======================
		for (a=u+1; a<partitions; ++a){
			x3N = nodeX[a][0][0].nodepos.y;
			dx_nod = x1N-x3N;
			for (b=0; b<partitions; ++b){
				y3N = nodeX[a][b][0].nodepos.y;
				dy_nod = y1N-y3N;
				for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod = z1N-z3N;
					dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
					if (dis_nod <= ddmax_nod){
					dz_nod = z2N-z3N;
					dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
					if (dis_nod <= ddmax_nod){
					count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
					}
					}
				}
			}
		}
		}
	}
	u = row;
	//=======================
	// Nodo 2 movil en ZY:
	//=======================
	for (v=col+1; v<partitions ; ++v){
		y2N = nodeX[u][v][0].nodepos.y;
		dy_nod = y1N-y2N;
		for (w=0; w<partitions ; ++w){		
			z2N = nodeX[u][v][w].nodepos.z;
			dz_nod = z1N-z2N;
			dis_nod = dy_nod*dy_nod + dz_nod*dz_nod;
			if (dis_nod <= ddmax_nod){
			//==============================================
			// 2 puntos en N1 y 1 punto en N2
			//==============================================
			count_3_N112(row, col, mom, u, v, w, XXX, nodeX);
			//==============================================
			// 2 puntos en N1 y 1 punto en N2
			//==============================================
			count_3_N122(row, col, mom, u, v, w, XXX, nodeX);
			//==============================================
			// 1 punto en N1, 1 punto en N2 y un punto en N3
			//==============================================
			a = u;
			b = v;
			//=======================
			// Nodo 3 movil en Z:
			//=======================
			y3N = nodeX[a][b][0].nodepos.y;
			dy_nod2 = y1N-y3N;
			for (c=w+1;  c<partitions; ++c){
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod = z1N-z3N;
				dis_nod = dy_nod2*dy_nod2 + dz_nod*dz_nod;
				if (dis_nod <= ddmax_nod){
				dz_nod = z2N-z3N;
				dis_nod = dz_nod*dz_nod;
				if (dis_nod <= ddmax_nod){
				count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
				}
				}
			}
			//=======================
			// Nodo 3 movil en ZY:
			//=======================	
			for (b=v+1; b<partitions; ++b){
				y3N = nodeX[a][b][0].nodepos.y;
				dy_nod2 = y1N-y3N;
				for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod = z1N-z3N;
					dis_nod = dy_nod2*dy_nod2 + dz_nod*dz_nod;
					if (dis_nod <= ddmax_nod){
					dy_nod = y2N-y3N;
					dz_nod = z2N-z3N;
					dis_nod = dy_nod2*dy_nod2 + dz_nod*dz_nod;
					if (dis_nod <= ddmax_nod){
					count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
					}
					}
				}
			}
			//=======================
			// Nodo 3 movil en ZYX:
			//=======================
			for (a=u+1; a<partitions; ++a){
				x3N = nodeX[a][0][0].nodepos.y;
				dx_nod = x1N-x3N;
				for (b=0; b<partitions; ++b){
					y3N = nodeX[a][b][0].nodepos.y;
					dy_nod2 = y1N-y3N;
					for (c=0;  c<partitions; ++c){
						z3N = nodeX[a][b][c].nodepos.z;
						dz_nod = z1N-z3N;
						dis_nod = dx_nod*dx_nod + dy_nod2*dy_nod2 + dz_nod*dz_nod;
						if (dis_nod <= ddmax_nod){
						dx_nod = x2N-x3N;
						dy_nod = y2N-y3N;
						dz_nod = z2N-z3N;
						dis_nod = dx_nod*dx_nod + dy_nod2*dy_nod2 + dz_nod*dz_nod;
						if (dis_nod <= ddmax_nod){
						count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
						}
						}
					}
				}
			}
			}
		}	
	}			
	//=======================
	// Nodo 2 movil en ZYX:
	//=======================
	for (u=row+1; u<partitions; ++u){
		x2N = nodeX[u][0][0].nodepos.x;
		dx_nod = x1N-x2N;
		for (v=0; v<partitions; ++v){
			y2N = nodeX[u][v][0].nodepos.y;
			dy_nod = y1N-y2N;
			for (w=0; w<partitions; ++w){
				z2N = nodeX[u][v][w].nodepos.z;
				dz_nod = z1N-z2N;
				dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
				if (dis_nod <= ddmax_nod){
				//==============================================
				// 2 puntos en N1 y 1 punto en N2
				//==============================================
				count_3_N112(row, col, mom, u, v, w, XXX, nodeX);
				//==============================================
				// 2 puntos en N1 y 1 punto en N2
				//==============================================
				count_3_N122(row, col, mom, u, v, w, XXX, nodeX);
				//==============================================
				// 1 punto en N1, 1 punto en N2 y un punto en N3
				//==============================================
				a = u;
				b = v;
				//=======================
				// Nodo 3 movil en Z:
				//=======================
				x3N = nodeX[a][b][0].nodepos.x;
				y3N = nodeX[a][b][0].nodepos.y;
				dx_nod2 = x1N-x3N;
				dy_nod2 = y1N-y3N;
				for (c=w+1;  c<partitions; ++c){	
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod = z1N-z3N;
					dis_nod = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod*dz_nod;
					if (dis_nod <= ddmax_nod){
					dz_nod = z2N-z3N;
					dis_nod = dz_nod*dz_nod;
					if (dis_nod <= ddmax_nod){
					count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
					}
					}
				}
				//=======================
				// Nodo 3 movil en ZY:
				//=======================
				x3N = nodeX[a][b][0].nodepos.x;
				dx_nod2 = x1N-x3N;
				for (b=v+1; b<partitions; ++b){
					y3N = nodeX[a][b][0].nodepos.y;
					dy_nod2 = y1N-y3N;
					for (c=0;  c<partitions; ++c){
						z3N = nodeX[a][b][c].nodepos.z;
						dz_nod = z1N-z3N;
						dis_nod = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod*dz_nod;
						if (dis_nod <= ddmax_nod){
						dy_nod2 = y2N-y3N;
						dz_nod = z2N-z3N;
						dis_nod = dy_nod2*dy_nod2 + dz_nod*dz_nod;
						if (dis_nod <= ddmax_nod){
						count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
						}
						}
					}
				}
				//=======================
				// Nodo 3 movil en ZYX:
				//=======================		
				for (a=u+1; a<partitions; ++a){
					x3N = nodeX[a][0][0].nodepos.y;
					dx_nod2 = x1N-x3N;
					for (b=0; b<partitions; ++b){
						y3N = nodeX[a][b][0].nodepos.y;
						dy_nod2 = y1N-y3N;
						for (c=0;  c<partitions; ++c){
							z3N = nodeX[a][b][c].nodepos.z;
							dz_nod = z1N-z3N;
							dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
							if (dis_nod <= ddmax_nod){
							dx_nod2 = x2N-x3N;
							dy_nod2 = y2N-y3N;
							dz_nod = z2N-z3N;
							dis_nod = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod*dz_nod;
							if (dis_nod <= ddmax_nod){
							count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
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
	}
	}

	// Histograma RR (ANALITICA)
	//======================================
	/*
	float dr = (d_max/bn);
	float alph = 4*(3.14159265359)*(n_pts*n_pts)*dr*dr*dr/(3*size_box*size_box*size_box);
	float r1;
	for(int a=0; a<bn; a++) {
		r1 = (a+1);
        	*(YY+a) += alph*(r1*r1*r1-a*a*a);
	}
	*/
}
//=================================================================== 

//=================================================================== 
void NODE::count_3_N111(int row, int col, int mom, unsigned int ***XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en un mismo Nodo.
	
	row, col, mom => posición del Nodo. Esto define al Nodo.
	
	*/
	
	int i,j,k;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	float dx,dy,dz;
	float d12,d13,d23;
	for (i= 0; i<nodeS[row][col][mom].len-2; ++i){
		x1 = nodeS[row][col][mom].elements[i].x;
		y1 = nodeS[row][col][mom].elements[i].y;
		z1 = nodeS[row][col][mom].elements[i].z;
		
		for (j=i+1; j<nodeS[row][col][mom].len-1; ++j){
			x2 = nodeS[row][col][mom].elements[j].x;
			y2 = nodeS[row][col][mom].elements[j].y;
			z2 = nodeS[row][col][mom].elements[j].z;
			dx = x2-x1;
			dy = y2-y1;
			dz = z2-z1;
			d12 = dx*dx+dy*dy+dz*dz;
			if (d12<=dd_max){
			for (k=j+1; k<nodeS[row][col][mom].len; ++k){ 
				x3 = nodeS[row][col][mom].elements[k].x;
				y3 = nodeS[row][col][mom].elements[k].y;
				z3 = nodeS[row][col][mom].elements[k].z;
				dx = x3-x1;
				dy = y3-y1;
				dz = z3-z1;
				d13 = dx*dx+dy*dy+dz*dz;
				if (d13 <= dd_max){
				dx = x2-x3;
				dy = y2-y3;
				dz = z2-z3;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23 <= dd_max){
				*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=1;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
void NODE::count_3_N112(int row, int col, int mom, int u, int v, int w, unsigned int ***XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en dos 
	nodos con dos puntos en N1 y un punto en N2.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	
	*/
	int i,j,k;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	float dx,dy,dz;
	float d12,d13,d23;
	for (i= 0; i<nodeS[row][col][mom].len-1; ++i){
		// 1er punto en N1
		x1 = nodeS[row][col][mom].elements[i].x;
		y1 = nodeS[row][col][mom].elements[i].y;
		z1 = nodeS[row][col][mom].elements[i].z;
		for (j=0; j<nodeS[u][v][w].len; ++j){
			// 2do punto en N2
			x2 = nodeS[u][v][w].elements[j].x;
			y2 = nodeS[u][v][w].elements[j].y;
			z2 = nodeS[u][v][w].elements[j].z;
			dx = x2-x1;
			dy = y2-y1;
			dz = z2-z1;
			d12 = dx*dx+dy*dy+dz*dz;
			if (d12<=dd_max){
			for (k=i+1; k<nodeS[row][col][mom].len; ++k){
				// 3er punto en N2
				x3 = nodeS[row][col][mom].elements[k].x;
				y3 = nodeS[row][col][mom].elements[k].y;
				z3 = nodeS[row][col][mom].elements[k].z;
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23 <= dd_max){
				dx = x3-x1;
				dy = y3-y1;
				dz = z3-z1;
				d13 = dx*dx+dy*dy+dz*dz;
				if (d13 <= dd_max){
				*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=1;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
void NODE::count_3_N122(int row, int col, int mom, int u, int v, int w, unsigned int ***XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en dos 
	nodos con dos puntos en N1 y un punto en N2.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	
	*/
	int i,j,k;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	float dx,dy,dz;
	float d12,d13,d23;
	for (i= 0; i<nodeS[row][col][mom].len; ++i){
		// 1er punto en N1
		x1 = nodeS[row][col][mom].elements[i].x;
		y1 = nodeS[row][col][mom].elements[i].y;
		z1 = nodeS[row][col][mom].elements[i].z;
		for (j=0; j<nodeS[u][v][w].len-1; ++j){
			// 2do punto en N2
			x2 = nodeS[u][v][w].elements[j].x;
			y2 = nodeS[u][v][w].elements[j].y;
			z2 = nodeS[u][v][w].elements[j].z;
			dx = x2-x1;
			dy = y2-y1;
			dz = z2-z1;
			d12 = dx*dx+dy*dy+dz*dz;
			if (d12<=dd_max){
			for (k=j+1; k<nodeS[row][col][mom].len; ++k){
				// 3er punto en N2
				x3 = nodeS[row][col][mom].elements[k].x;
				y3 = nodeS[row][col][mom].elements[k].y;
				z3 = nodeS[row][col][mom].elements[k].z;
				dx = x3-x1;
				dy = y3-y1;
				dz = z3-z1;
				d13 = dx*dx+dy*dy+dz*dz;
				if (d13<=dd_max){
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23 <= dd_max){
				*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=1;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
void NODE::count_3_N123(int row, int col, int mom, int u, int v, int w, int a, int b, int c, unsigned int ***XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en tres 
	nodos con un puntos en N1, un punto en N2
	y un punto en N3.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	a, b, c => posición de N3.
	
	*/
	int i,j,k;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	float dx,dy,dz;
	float d12,d13,d23;
	for (i=0; i<nodeS[row][col][mom].len; ++i){
		// 1er punto en N1
		x1 = nodeS[row][col][mom].elements[i].x;
		y1 = nodeS[row][col][mom].elements[i].y;
		z1 = nodeS[row][col][mom].elements[i].z;
		for (j=0; j<nodeS[a][b][c].len; ++j){
			// 3do punto en N3
			x3 = nodeS[a][b][c].elements[j].x;
			y3 = nodeS[a][b][c].elements[j].y;
			z3 = nodeS[a][b][c].elements[j].z;
			dx = x3-x1;
			dy = y3-y1;
			dz = z3-z1;
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<=dd_max){
			for (k=0; k<nodeS[u][v][w].len; ++k){
				// 2er punto en N2
				x2 = nodeS[u][v][w].elements[k].x;
				y2 = nodeS[u][v][w].elements[k].y;
				z2 = nodeS[u][v][w].elements[k].z;
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23<=dd_max){
				dx = x2-x1;
				dy = y2-y1;
				dz = z2-z1;
				d12 = dx*dx+dy*dy+dz*dz;
				if (d12 <= dd_max){
				*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=1;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
/*
void NODE::histo_front_XXX(unsigned int *PP, Node ***dat, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int _row, int _col, int _mom, int _u, int _v, int _w){
	int i, j;
	float dis_f,_dis,_d_x,_d_y,_d_z;
	float _x,_y,_z;
	//======================================================================
	// Si los puentos estás en las paredes laterales de X
	if( con_in_x ){
		// forma de calcular la distancia a las proyecciones usando la distancia entre puntos dentro de la caja
		dis_f = disn + ll - 2*dn_x*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; ++i){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; ++j){
					_d_x = fabs(_x-dat[_u][_v][_w].elements[j].x)-size_box;
					_d_y = _y-dat[_u][_v][_w].elements[j].y;
					_d_z = _z-dat[_u][_v][_w].elements[j].z;
					_dis = (_d_x*_d_x) + (_d_y*_d_y) + (_d_z*_d_z); 
					if (_dis <= dd_max){
						*(PP + (int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las paredes laterales de Y		
	if( con_in_y ){
		dis_f = disn + ll - 2*dn_y*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; ++i){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; ++j){
					_d_x = _x-dat[_u][_v][_w].elements[j].x;
					_d_y = fabs(_y-dat[_u][_v][_w].elements[j].y)-size_box;
					_d_z = _z-dat[_u][_v][_w].elements[j].z;
					_dis = (_d_x*_d_x) + (_d_y*_d_y) + (_d_z*_d_z); 
					if (_dis <= dd_max){
						*(PP + (int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las paredes laterales de Z
	if( con_in_z ){
		dis_f = disn + ll - 2*dn_z*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; ++i){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; ++j){
					_d_x = _x-dat[_u][_v][_w].elements[j].x;
					_d_y = _y-dat[_u][_v][_w].elements[j].y;
					_d_z = fabs(_z-dat[_u][_v][_w].elements[j].z)-size_box;
					_dis = (_d_x*_d_x) + (_d_y*_d_y) + (_d_z*_d_z); 
					if (_dis <= dd_max){
						*(PP + (int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X y Y			
	if( con_in_x && con_in_y ){
		dis_f = disn + 2*ll - 2*(dn_x+dn_y)*size_box;
		if (dis_f < ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; ++i){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; ++j){
					_d_x = fabs(_x-dat[_u][_v][_w].elements[j].x)-size_box;
					_d_y = fabs(_y-dat[_u][_v][_w].elements[j].y)-size_box;
					_d_z = _z-dat[_u][_v][_w].elements[j].z;
					_dis = (_d_x*_d_x) + (_d_y*_d_y) + (_d_z*_d_z); 
					if (_dis <= dd_max){
						*(PP + (int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X y Z				
	if( con_in_x && con_in_z ){
		dis_f = disn + 2*ll - 2*(dn_x+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; ++i){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; ++j){
					_d_x = fabs(_x-dat[_u][_v][_w].elements[j].x)-size_box;
					_d_y = _y-dat[_u][_v][_w].elements[j].y;
					_d_z = fabs(_z-dat[_u][_v][_w].elements[j].z)-size_box;
					_dis = (_d_x*_d_x) + (_d_y*_d_y) + (_d_z*_d_z); 
					if (_dis <= dd_max){
						*(PP + (int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de Y y Z			
	if( con_in_y && con_in_z ){
		dis_f = disn + 2*ll - 2*(dn_y+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; ++i){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; ++j){
					_d_x = _x-dat[_u][_v][_w].elements[j].x;
					_d_y = fabs(_y-dat[_u][_v][_w].elements[j].y)-size_box;
					_d_z = fabs(_z-dat[_u][_v][_w].elements[j].z)-size_box;
					_dis = (_d_x*_d_x) + (_d_y*_d_y) + (_d_z*_d_z); 
					if (_dis <= dd_max){
						*(PP + (int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X, Y y Z		
	if( con_in_x && con_in_y && con_in_z ){
		dis_f = disn + 3*ll - 2*(dn_x+dn_y+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; ++i){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; ++j){
					_d_x = fabs(_x-dat[_u][_v][_w].elements[j].x)-size_box;
					_d_y = fabs(_y-dat[_u][_v][_w].elements[j].y)-size_box;
					_d_z = fabs(_z-dat[_u][_v][_w].elements[j].z)-size_box;
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(PP + (int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
}
*/
//====================================================================
/*
void NODE::histo_front_XY(unsigned int *XY, Node ***dat, Node ***ran, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int _row, int _col, int _mom, int _u, int _v, int _w){
	int _pos, i, j;
	float dis_f, _dis, _d_x, _d_y, _d_z;
	
	//======================================================================
	// Si los puentos estás en las paredes laterales de X
	if( con_in_x ){
		// forma de calcular la distancia a las proyecciones usando la distancia entre puntos dentro de la caja
		dis_f = disn + ll - 2*dn_x*size_box; 
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-fabs(dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x);
					_d_y =  dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y;
					_d_z =  dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z;
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z; 
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las paredes laterales de Y				
	if( con_in_y ){
		dis_f = disn + ll - 2*dn_y*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x;
					_d_y =  size_box-fabs(dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y);
					_d_z =  dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z;
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las paredes laterales de Z	
	if( con_in_z ){
		dis_f = disn + ll - 2*dn_z*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x;
					_d_y =  dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y;
					_d_z =  size_box-fabs(dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z);
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X y Y				
	if( con_in_x && con_in_y ){
		dis_f = disn + 2*ll - 2*(dn_x+dn_y)*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-fabs(dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x);
					_d_y =  size_box-fabs(dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y);
					_d_z =  dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z;
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X y Z	
	if( con_in_x && con_in_z ){
		dis_f = disn + 2*ll - 2*(dn_x+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-fabs(dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x);
					_d_y =  dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y;
					_d_z =  size_box-fabs(dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z);
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de Y y Z		
	if( con_in_y && con_in_z ){
		dis_f = disn + 2*ll - 2*(dn_y+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x;
					_d_y =  size_box-fabs(dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y);
					_d_z =  size_box-fabs(dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z);
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X, Y y Z		
	if( con_in_x && con_in_y && con_in_z ){
		dis_f = disn + 3*ll - 2*(dn_x+dn_y+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-fabs(dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x);
					_d_y =  size_box-fabs(dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y);
					_d_z =  size_box-fabs(dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z);
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
}
*/
NODE::~NODE(){
	
}
