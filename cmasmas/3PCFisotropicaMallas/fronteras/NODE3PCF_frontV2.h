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

class NODE3P{
	//Atributos de clase:
	private:
		// Asignados
		int bn;
		int n_pts;
		float size_box;
		float size_node;
		float d_max;
		Node ***nodeD;
		Point3D *dataD;
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
		NODE3P(int _bn, int _n_pts, float _size_box, float _size_node, float _d_max, Point3D *_dataD, Node ***_nodeD){
			bn = _bn;
			n_pts = _n_pts;
			size_box = _size_box;
			size_node = _size_node;
			d_max = _d_max;
			dataD = _dataD;
			nodeD = _nodeD;
			ll = size_box*size_box;
			dd_max = d_max*d_max;
			front = size_box - d_max;
			corr = size_node*sqrt(3);
			ds = ((float)(bn))/d_max;
			ddmax_nod = (d_max+corr)*(d_max+corr);
			make_nodos(nodeD,dataD);
			std::cout << "Terminé de construir nodos..." << std::endl;
		}
		
		Node ***meshData(){
			return nodeD;
		};
		
		// Implementamos Método de mallas:
		void make_histoXXX(unsigned int ***, Node ***);
		void make_histo_analitic(unsigned int ***, unsigned int ***, unsigned int, Node ***);
		void count_3_N111(int, int, int, unsigned int ***, Node ***);
		void count_3_N112(int, int, int, int, int, int, unsigned int ***, Node ***);
		void count_3_N123(int, int, int, int, int, int, int, int, int, unsigned int ***, Node ***);
		void symmetrize(unsigned int ***);
		void make_fron(unsigned int *** , Node *&, Node ***, int &);
		void add_front(Node *&, Point3D *&, int, float, float, float, short int, short int, short int, int&);
		void count_3_N112_front(int, int, int, int, int, int, unsigned int ***, Node ***, short int, short int, short int);
		void count_3_N123_front(int, int, int, int, int, int, int, int, int, unsigned int ***, Node ***, short int, short int, short int);
		void count_3_N123_front_FFN(int, int, int, int, int, unsigned int ***, Node ***, Node *);
		~NODE3P();
};

//=================================================================== 
//==================== Funciones ==================================== 
//=================================================================== 

void NODE3P::make_histoXXX(unsigned int ***XXX, Node ***nodeX){
	/*
	Función para crear los histogramas DDD y RRR.
	
	Argumentos
	DDD: arreglo donde se creará el histograma DDD.
	
	*/ 
	int partitions = (int)(ceil(size_box/size_node));
	int i, j, row, col, mom, u, v, w, a ,b, c;
	float dis, dis_nod, dis_nod2, dis_nod3;
	float x1N, y1N, z1N, x2N, y2N, z2N, x3N, y3N, z3N;
	float x, y, z;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod, dx_nod2, dy_nod2, dz_nod2, dx_nod3, dy_nod3, dz_nod3;
	bool con_x, con_y, con_z;
	float d_max_pm = d_max + size_node/2, front_pm = front - size_node/2;
	
	std::cout << "-> Estoy haciendo histograma XXX..." << std::endl;
	
	// x1N, y1N, z1N => Nodo pivote
	for (row=0; row<partitions; ++row){
	x1N = nodeX[row][0][0].nodepos.x;
	for (col=0; col<partitions; ++col){
	y1N = nodeX[row][col][0].nodepos.y;
	for (mom=0; mom<partitions; ++mom){
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
	x2N = nodeX[u][0][0].nodepos.x;
	y2N = nodeX[u][v][0].nodepos.y;
	for (w=mom+1;  w<partitions; ++w){	
		z2N = nodeX[u][v][w].nodepos.z;
		dz_nod = z2N-z1N;
		dis_nod = dz_nod*dz_nod;
		if (dis_nod <= ddmax_nod){
		//==============================================
		// 2 puntos en N y 1 punto en N'
		//==============================================
		count_3_N112(row, col, mom, u, v, w, XXX, nodeX);
		//==============================================
		// 1 punto en N1, 1 punto en N2 y 1 punto en N3
		//==============================================
		a = u;
		b = v;
		//=======================
		// Nodo 3 movil en Z:
		//=======================
		for (c=w+1;  c<partitions; ++c){
			z3N = nodeX[a][b][c].nodepos.z; 
			dz_nod2 = z3N-z1N;
			dis_nod2 = dz_nod2*dz_nod2;
			if (dis_nod2 <= ddmax_nod){
			count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
			}
		}
		
		//=======================
		// Nodo 3 movil en ZY:
		//=======================
		for (b=v+1; b<partitions; ++b){
			y3N = nodeX[a][b][0].nodepos.y;
			dy_nod2 = y3N-y1N;
			for (c=0;  c<partitions; ++c){
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod2 = z3N-z1N;
				dis_nod2 = dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
				if (dis_nod2 <= ddmax_nod){
				dy_nod3 = y3N-y2N;
				dz_nod3 = z3N-z2N;
				dis_nod3 = dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
				if (dis_nod3 <= ddmax_nod){
				count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
				}
				}
			}
		}
		//=======================
		// Nodo 3 movil en ZYX:
		//=======================
		for (a=u+1; a<partitions; ++a){
			x3N = nodeX[a][0][0].nodepos.x;
			dx_nod2 = x3N-x1N;
			for (b=0; b<partitions; ++b){
				y3N = nodeX[a][b][0].nodepos.y;
				dy_nod2 = y3N-y1N;
				for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
					if (dis_nod2 <= ddmax_nod){
					dx_nod3 = x3N-x2N;
					dy_nod3 = y3N-y2N;
					dz_nod3 = z3N-z2N;
					dis_nod3 = dx_nod3*dx_nod3 + dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
					if (dis_nod3 <= ddmax_nod){
					count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
					}
					}
				}
			}
		}
		}
	}
	//=======================
	// Nodo 2 movil en ZY:
	//=======================
	x2N = nodeX[u][0][0].nodepos.x;
	for (v=col+1; v<partitions ; ++v){
		y2N = nodeX[u][v][0].nodepos.y;
		dy_nod = y2N-y1N;
		for (w=0; w<partitions ; ++w){		
			z2N = nodeX[u][v][w].nodepos.z;
			dz_nod = z2N-z1N;
			dis_nod = dy_nod*dy_nod + dz_nod*dz_nod;
			if (dis_nod <= ddmax_nod){
			//==============================================
			// 2 puntos en N y 1 punto en N'
			//==============================================
			count_3_N112(row, col, mom, u, v, w, XXX, nodeX);
			//==============================================
			// 1 punto en N1, 1 punto en N2 y un punto en N3
			//==============================================
			a = u;
			b = v;
			//=======================
			// Nodo 3 movil en Z:
			//=======================
			y3N = nodeX[a][b][0].nodepos.y;
			dy_nod2 = y3N-y1N;
			for (c=w+1;  c<partitions; ++c){
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod2 = z3N-z1N;
				dis_nod2 = dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
				if (dis_nod2 <= ddmax_nod){
				dz_nod3 = z3N-z2N;
				dis_nod3 = dz_nod3*dz_nod3;
				if (dis_nod3 <= ddmax_nod){
				count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
				}
				}
			}
			//=======================
			// Nodo 3 movil en ZY:
			//=======================	
			for (b=v+1; b<partitions; ++b){
				y3N = nodeX[a][b][0].nodepos.y;
				dy_nod2 = y3N-y1N;
				for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					dis_nod2 = dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
					if (dis_nod2 <= ddmax_nod){
					dy_nod3 = y3N-y2N;
					dz_nod3 = z3N-z2N;
					dis_nod3 = dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
					if (dis_nod3 <= ddmax_nod){
					count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
					}
					}
				}
			}
			//=======================
			// Nodo 3 movil en ZYX:
			//=======================
			for (a=u+1; a<partitions; ++a){
				x3N = nodeX[a][0][0].nodepos.x;
				dx_nod2 = x3N-x1N;
				for (b=0; b<partitions; ++b){
					y3N = nodeX[a][b][0].nodepos.y;
					dy_nod2 = y3N-y1N;
					for (c=0;  c<partitions; ++c){
						z3N = nodeX[a][b][c].nodepos.z;
						dz_nod2 = z3N-z1N;
						dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
						if (dis_nod2 <= ddmax_nod){
						dx_nod3 = x3N-x2N;
						dy_nod3 = y3N-y2N;
						dz_nod3 = z3N-z2N;
						dis_nod3 = dx_nod3*dx_nod3 + dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
						if (dis_nod3 <= ddmax_nod){
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
		dx_nod = x2N-x1N;
		for (v=0; v<partitions; ++v){
			y2N = nodeX[u][v][0].nodepos.y;
			dy_nod = y2N-y1N;
			for (w=0; w<partitions; ++w){
				z2N = nodeX[u][v][w].nodepos.z;
				dz_nod = z2N-z1N;
				dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
				if (dis_nod <= ddmax_nod){
				//==============================================
				// 2 puntos en N y 1 punto en N'
				//==============================================
				count_3_N112(row, col, mom, u, v, w, XXX, nodeX);
				//==============================================
				// 1 punto en N1, 1 punto en N2 y 1 punto en N3
				//==============================================
				a = u;
				b = v;
				//=======================
				// Nodo 3 movil en Z:
				//=======================
				x3N = nodeX[a][0][0].nodepos.x;
				y3N = nodeX[a][b][0].nodepos.y;
				dx_nod2 = x3N-x1N;
				dy_nod2 = y3N-y1N;
				for (c=w+1;  c<partitions; ++c){	
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
					if (dis_nod2 <= ddmax_nod){
					dz_nod3 = z3N-z2N;
					dis_nod3 = dz_nod3*dz_nod3;
					if (dis_nod3 <= ddmax_nod){
					count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
					}
					}
				}
				//=======================
				// Nodo 3 movil en ZY:
				//=======================
				x3N = nodeX[a][0][0].nodepos.x;
				dx_nod2 = x3N-x1N;
				for (b=v+1; b<partitions; ++b){
					y3N = nodeX[a][b][0].nodepos.y;
					dy_nod2 = y3N-y1N;
					for (c=0;  c<partitions; ++c){
						z3N = nodeX[a][b][c].nodepos.z;
						dz_nod2 = z3N-z1N;
						dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
						if (dis_nod2 <= ddmax_nod){
						dy_nod3 = y3N-y2N;
						dz_nod3 = z3N-z2N;
						dis_nod3 = dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
						if (dis_nod3 <= ddmax_nod){
						count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
						}
						}
					}
				}
				//=======================
				// Nodo 3 movil en ZYX:
				//=======================		
				for (a=u+1; a<partitions; ++a){
					x3N = nodeX[a][0][0].nodepos.x;
					dx_nod2 = x3N-x1N;
					for (b=0; b<partitions; ++b){
						y3N = nodeX[a][b][0].nodepos.y;
						dy_nod2 = y3N-y1N;
						for (c=0;  c<partitions; ++c){
							z3N = nodeX[a][b][c].nodepos.z;
							dz_nod2 = z3N-z1N;
							dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
							if (dis_nod2 <= ddmax_nod){
							dx_nod3 = x3N-x2N;
							dy_nod3 = y3N-y2N;
							dz_nod3 = z3N-z2N;
							dis_nod3 = dx_nod3*dx_nod3 + dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
							if (dis_nod3 <= ddmax_nod){
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
	//================================
	// Simetrización:
	//================================
	symmetrize(XXX); 
	//===============================================================================
	// Condiciones Preriódicas de frontera
	//===============================================================================
	
	int part_front = 0;
	Node *nodeFront;
	nodeFront = new Node[875];
	int n, m;
	make_fron(XXX, nodeFront, nodeX, part_front);
	//std::cout <<  nodeFront[0].len<< std::endl;
	for (n=0; n<part_front; ++n){
		// 1er punto en N1 (frontera)
		x1N = nodeFront[n].nodepos.x;
		y1N = nodeFront[n].nodepos.y;
		z1N = nodeFront[n].nodepos.z;
		for (m=n+1; m<part_front; ++m){
			// 2do punto en N3 (caja)
			x2N = nodeFront[m].nodepos.x;
			y2N = nodeFront[m].nodepos.y;
			z2N = nodeFront[m].nodepos.z;
			dx = x2N-x1N;
			dy = y2N-y1N;
			dz = z2N-z1N;
			dis_nod = dx*dx + dy*dy + dz*dz;
			
			if (dis_nod <= ddmax_nod){
			for (u=0; u<partitions; ++u){
			x3N = nodeX[u][0][0].nodepos.x; 
			for (v=0; v<partitions; ++v){
			y3N = nodeX[u][v][0].nodepos.y;
			for (w=0; w<partitions; ++w){
			z3N = nodeX[u][v][w].nodepos.z;
			dx = x3N-x1N;
			dy = y3N-y1N;
			dz = z3N-z1N;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = x3N-x2N;
			dy = y3N-y2N;
			dz = z3N-z2N;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front_FFN( n, m, u, v, w, XXX, nodeX, nodeFront);
			}
			}	
			} } }
			}
		}
	}	
}
//=================================================================== 
void NODE3P::make_histo_analitic(unsigned int ***XXX, unsigned int ***XYY, unsigned int XXY, Node ***nodeX){
	
	// Histograma RR (ANALITICA)
	//======================================
	
	int i, j, k;
	
	float dr = (d_max/bn);
	float alph = 8*(3.14159265359)*(3.14159265359)*(n_pts*n_pts*n_pts)*dr*dr*dr*dr*dr*dr/(size_box*size_box);
	float r1;
	for(i=0; i<bn; i++) {
	for(j=0; j<bn; j++) {
	for(k=0; k<bn; k++) {
        	*(*(*(XXX+i)+j)+k)+=alph*(0.5+i)*(0.5+j)*(0.5+k);
	}
	}
	}
}
//=================================================================== 
void NODE3P::symmetrize(unsigned int ***XXX){
	int i,j,k;
	float elem;
	for (i=0; i<bn; i++){
	for (j=i; j<bn; j++){
	for (k=j; k<bn; k++){
		elem = XXX[i][j][k] + XXX[k][i][j] + XXX[j][k][i] + XXX[j][i][k] + XXX[k][j][i] + XXX[i][k][j];
		XXX[i][j][k] = elem;
		XXX[k][i][j] = elem;
		XXX[j][k][i] = elem;
		XXX[j][i][k] = elem;
		XXX[k][j][i] = elem;
		XXX[i][k][j] = elem;
	}   
	}
	}
}
//=================================================================== 
void NODE3P::count_3_N111(int row, int col, int mom, unsigned int ***XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en un mismo Nodo.
	
	row, col, mom => posición del Nodo. Esto define al Nodo.
	
	*/
	int i,j,k;
	float dx,dy,dz;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	int a, b, c;
	for (i=0; i<nodeS[row][col][mom].len-2; ++i){
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
				if (d13<=dd_max){
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23<=dd_max){
				*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=1;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N112(int row, int col, int mom, int u, int v, int w, unsigned int ***XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en dos 
	nodos con dos puntos en N1 y un punto en N2.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	
	*/
	int i,j,k;
	float dx,dy,dz;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	int a, b, c;
	for (i=0; i<nodeS[u][v][w].len; ++i){
		// 1er punto en N2
		x1 = nodeS[u][v][w].elements[i].x;
		y1 = nodeS[u][v][w].elements[i].y;
		z1 = nodeS[u][v][w].elements[i].z;
		for (j=0; j<nodeS[row][col][mom].len; ++j){
			// 2do punto en N1
			x2 = nodeS[row][col][mom].elements[j].x;
			y2 = nodeS[row][col][mom].elements[j].y;
			z2 = nodeS[row][col][mom].elements[j].z;
			dx = x2-x1;
			dy = y2-y1;
			dz = z2-z1;
			d12 = dx*dx+dy*dy+dz*dz;
			if (d12<=dd_max){
			for (k=j+1; k<nodeS[row][col][mom].len; ++k){
				// 3er punto en N1
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
				if (d23<=dd_max){
				*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=1;
				}
				}
			}
			for (k=i+1; k<nodeS[u][v][w].len; ++k){
				// 3er punto en N2
				x3 = nodeS[u][v][w].elements[k].x;
				y3 = nodeS[u][v][w].elements[k].y;
				z3 = nodeS[u][v][w].elements[k].z;
				dx = x3-x1;
				dy = y3-y1;
				dz = z3-z1;
				d13 = dx*dx+dy*dy+dz*dz;
				if (d13<=dd_max){
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23<=dd_max){
				*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=1;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N123(int row, int col, int mom, int u, int v, int w, int a, int b, int c, unsigned int ***XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en tres 
	nodos con un puntos en N1, un punto en N2
	y un punto en N3.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	a, b, c => posición de N3.
	
	*/
	int i,j,k, n, m, l;
	float dx,dy,dz;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	for (i=0; i<nodeS[row][col][mom].len; ++i){
		// 1er punto en N1
		x1 = nodeS[row][col][mom].elements[i].x;
		y1 = nodeS[row][col][mom].elements[i].y;
		z1 = nodeS[row][col][mom].elements[i].z;
		for (j=0; j<nodeS[a][b][c].len; ++j){
			// 2do punto en N3
			x3 = nodeS[a][b][c].elements[j].x;
			y3 = nodeS[a][b][c].elements[j].y;
			z3 = nodeS[a][b][c].elements[j].z;
			dx = x3-x1;
			dy = y3-y1;
			dz = z3-z1;
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<=dd_max){
			for (k=0; k<nodeS[u][v][w].len; ++k){
				// 3er punto en N2
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
				if (d12<=dd_max){
				*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=1;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::make_nodos(Node ***nod, Point3D *dat){
	/*
	Función para crear los nodos con los datos y puntos random
	
	Argumentos
	nod: arreglo donde se crean los nodos.
	dat: datos a dividir en nodos.
	
	*/
	int i, row, col, mom, partitions = (int)(ceil(size_box/size_node));
	float p_med = size_node/2;
	
	// Inicializamos los nodos vacíos:
	for (row=0; row<partitions; row++){
		for (col=0; col<partitions; col++){
			for (mom=0; mom<partitions; mom++){
				nod[row][col][mom].nodepos.x = ((float)(row)*(size_node))+p_med;
				nod[row][col][mom].nodepos.y = ((float)(col)*(size_node))+p_med;
				nod[row][col][mom].nodepos.z = ((float)(mom)*(size_node))+p_med;
				nod[row][col][mom].len = 0;
				nod[row][col][mom].elements = new Point3D[0];
			}
		}
	}
	// Llenamos los nodos con los puntos de dat:
	for (i=0; i<n_pts; i++){
		row = (int)(floor(dat[i].x/size_node));
        	col = (int)(floor(dat[i].y/size_node));
        	mom = (int)(floor(dat[i].z/size_node));
		add(nod[row][col][mom].elements, nod[row][col][mom].len, dat[i].x, dat[i].y, dat[i].z);
	}
}
//=================================================================== 
void NODE3P::add(Point3D *&array, int &lon, float _x, float _y, float _z){
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
void NODE3P::add_front(Node *&front, Point3D *&array, int ln, float xN, float yN, float zN, short int sig1, short int sig2, short int sig3, int &n){
	/* 
	Intento de función para crear copias proyecciones de los nodos en la frontera
	*/
	int i, j, k;
	
	Node *aux = new Node[n+1];
	for (i=0; i<n; i++){
		aux[i].nodepos.x = front[i].nodepos.x;
		aux[i].nodepos.y = front[i].nodepos.y;
		aux[i].nodepos.z = front[i].nodepos.z;
		aux[i].len = front[i].len;
		aux[i].elements = new Point3D[front[i].len];
		for (k=0; k<front[i].len; k++){
			aux[i].elements[k].x = front[i].elements[k].x;
			aux[i].elements[k].y = front[i].elements[k].y;
			aux[i].elements[k].z = front[i].elements[k].z;
		}
	}
	delete[] front;
	front = aux;
	front[n].nodepos.x = xN+(size_box*(sig1)); // el sig es para ubicar la proyeccion fuera de la caja 
	front[n].nodepos.y = yN+(size_box*(sig2));
	front[n].nodepos.z = zN+(size_box*(sig3));
	front[n].len = ln;
	front[n].elements = new Point3D[ln];
	for (i=0; i<ln; i++){
		front[n].elements[i].x = array[i].x+(size_box*(sig1));
		front[n].elements[i].y = array[i].y+(size_box*(sig2));
		front[n].elements[i].z = array[i].z+(size_box*(sig3));
	}
	n++;
}
//=================================================================== 
void NODE3P::make_fron(unsigned int ***XXX , Node *&nodeFront, Node ***node, int &n){
	/*
	Función para crear los nodos proyectados en la frontera de la caja.
	
	Argumentos
	nodeFront: arreglo de nodos proyectados.
	node: nodos en la caja.
	
	NOTA: "los nodos de la frontera NO tienen un orden especial".
	*/
	int row, col, mom, i, j, k, u, v, w;
	int partitions = (int)(ceil(size_box/size_node));
	int lim_n1 = (int)(ceil(d_max/size_node))-1, lim_n2 = (int)(ceil(d_max/size_node));
	float x1N, y1N, z1N, xN, yN, zN, xN_, yN_, zN_, xF, yF, zF;
	float len;
	float dx, dy, dz;
	float dis_nod;
	float d_max_pm = d_max + size_node/2, front_pm = front - size_node/2;
	bool con_xup, con_xdown, con_yup, con_ydown, con_zup, con_zdown, con_x, con_y, con_z;

	for (row=0; row<partitions; ++row){
	x1N = node[row][0][0].nodepos.x; 
	con_xdown = x1N<=d_max_pm;
	con_xup = x1N>=front_pm;
	for (col=0; col<partitions; ++col){
	y1N = node[row][col][0].nodepos.y;
	con_ydown = y1N<=d_max_pm;
	con_yup = y1N>=front_pm;
	for (mom=0; mom<partitions; ++mom){
	z1N = node[row][col][mom].nodepos.z;	
	con_zdown = z1N<=d_max_pm;
	con_zup = z1N>=front_pm;
	
	con_x = con_xup||con_xdown;
	con_y = con_yup||con_ydown;
	con_z = con_zup||con_zdown;
	
	if(con_x){ // Nodos en paredes X
		//===================================================================================
		// Boloque 1
		//===================================================================================
		if (con_xup){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, -1, 0, 0, n);
		xF = x1N-size_box; 
		yF = y1N;
		zF = z1N;
		for(i=0; i<=lim_n1; ++i){
		xN = node[i][0][0].nodepos.x; // Nodos en la caja
		for(j=0; j<partitions; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<partitions; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, -1, 0, 0);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1;  w<partitions ; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, 0, 0);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<partitions ; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w = 0; w < partitions ; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, 0, 0);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<=lim_n1; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, 0, 0);
			}
			}
			}}}
		}	
		} } }
	 	}
	 	//===================================================================================
	 	// Boloque 2
	 	//===================================================================================
	 	else if (con_xdown){
	 	len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, +1, 0, 0, n);
		xF = x1N+size_box; 
		yF = y1N;
		zF = z1N;
		for(i=partitions-1; i>=partitions-lim_n2; --i){
		xN = node[i][0][0].nodepos.x; 
		for(j=0; j<partitions; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<partitions; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, +1, 0, 0);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1;  w<partitions ; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, 0, 0);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<partitions ; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions ; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, 0, 0);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i-1; u>=partitions-lim_n2; --u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, 0, 0);
			}
			}
			}}}
		}	
		} } } 
		}
	}
	if(con_y){ // Nodos en paredes Y
		//===================================================================================
		// Boloque 3
		//===================================================================================
		if (con_yup){
		len = node[row][col][mom].len;
		//Haacemos una copia en la pared Y contraria:
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, 0, -1, 0, n);
		xF = x1N; 
		yF = y1N-size_box;
		zF = z1N;
		for(i=0; i<partitions; ++i){
		xN = node[i][0][0].nodepos.x; 
		for(j=0; j<=lim_n1; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<partitions; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, 0, -1, 0);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1;  w<partitions ; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, -1, 0);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions ; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, -1, 0);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<partitions; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, -1, 0);
			}
			}
			}}}
		}
		} } }
	 	}
	 	//===================================================================================
	 	// Boloque 4
		//===================================================================================
	 	else if (con_ydown){
	 	len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, 0, +1, 0, n);
		xF = x1N;
		yF = y1N+size_box;
		zF = z1N;
		for(i=0; i<partitions; ++i){
		xN = node[i][0][0].nodepos.x; 
		for(j=partitions-1; j>=partitions-lim_n2; --j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<partitions; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, 0, +1, 0);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1;  w<partitions ; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, +1, 0);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions ; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, +1, 0);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<partitions; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=partitions-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, +1, 0);
			}
			}
			}}}
		}	
		} } }
	 	}
	}
	if(con_z){ // Nodos en paredes Z
		//===================================================================================
		// Boloque 5
		//===================================================================================
		if (con_zup){
		len = node[row][col][mom].len;
		//Haacemos una copia en la pared Z contraria:
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, 0, 0, -1, n);
		xF = x1N; 
		yF = y1N;
		zF = z1N-size_box;
		for(i=0; i<partitions; ++i){
		xN = node[i][0][0].nodepos.x; 
		for(j=0; j<partitions; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<=lim_n1; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, 0, 0, -1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1;  w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, 0, -1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, 0, -1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<partitions; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, 0, -1);
			}
			}
			}}}
		}	
		} } }
	 	}
	 	//===================================================================================
	 	// Boloque 6
		//===================================================================================
	 	else if (con_zdown){
	 	len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, 0, 0, +1, n);
		xF = x1N;
		yF = y1N;
		zF = z1N+size_box; 
		for(i=0; i<partitions; ++i){
		xN = node[i][0][0].nodepos.x; 
		for(j=0; j<partitions; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=partitions-1; k>=partitions-lim_n2; --k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, 0, 0, +1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, 0, +1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, 0, +1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<partitions; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, 0, +1);
			}
			}
			}}}
		}	
		} } }
	 	}
	}
	if(con_x && con_y){ // Nodos en esquinas XY
		//===================================================================================
		// Boloque 7
		//===================================================================================
		if (con_xup && con_yup){
		len = node[row][col][mom].len;
		//Haacemos una copia en la esquina XY contraria:
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, -1, -1, 0, n);
		xF = x1N-size_box; 
		yF = y1N-size_box;
		zF = z1N;
		for(i=0; i<=lim_n1; ++i){
		xN = node[i][0][0].nodepos.x;
		for(j=0; j<=lim_n1; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<partitions; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, -1, -1, 0);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, -1, 0);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, -1, 0);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<=lim_n1; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, -1, 0);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 8
		//===================================================================================
		else if (con_xup && con_ydown){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, -1, +1, 0, n);
		xF = x1N-size_box;
		yF = y1N+size_box;
		zF = z1N;
		for(i=0; i<=lim_n1; ++i){
		xN = node[i][0][0].nodepos.x; 
		for(j=partitions-1; j>=partitions-lim_n2; --j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<partitions; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, -1, +1, 0);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, +1, 0);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, +1, 0);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<=lim_n1; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=partitions-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, +1, 0);
			}
			}
			}}}
		}	
		} } }
	 	}
		//===================================================================================
		// Boloque 9
		//===================================================================================
		else if (con_xdown && con_yup){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, +1, -1, 0, n);
		xF = x1N+size_box; 
		yF = y1N-size_box;
		zF = z1N;
		for(i=partitions-1; i>=partitions-lim_n2; --i){
		xN = node[i][0][0].nodepos.x; 
		for(j=0; j<=lim_n1; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<partitions; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, +1, -1, 0);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, -1, 0);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, -1, 0);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i-1; u>=partitions-lim_n2; --u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, -1, 0);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 10
		//===================================================================================
		else if (con_xdown && con_ydown){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, +1, +1, 0, n);
		xF = x1N+size_box; 
		yF = y1N+size_box;
		zF = z1N;
		for(i=partitions-1; i>=partitions-lim_n2; --i){
		xN = node[i][0][0].nodepos.x; 
		for(j=partitions-1; j>=partitions-lim_n2; --j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<partitions; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, +1, +1, 0);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, +1, 0);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, +1, 0);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i-1; u>=partitions-lim_n2; --u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=partitions-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<partitions; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, +1, 0);
			}
			}
			}}}
		}	
		} } }
		}
	}
	if(con_x && con_z){ // Nodos en esquinas XZ
		//===================================================================================
		// Boloque 11
		//===================================================================================
		if (con_xup && con_zup){
		len = node[row][col][mom].len;
		//Haacemos una copia en la esquina XZ contraria:
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, -1, 0, -1, n);
		xF = x1N-size_box; 
		yF = y1N;
		zF = z1N-size_box;
		for(i=0; i<=lim_n1; ++i){
		xN = node[i][0][0].nodepos.x;
		for(j=0; j<partitions; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<=lim_n1; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, -1, 0, -1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, 0, -1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, 0, -1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<=lim_n1; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, 0, -1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 12
		//===================================================================================
		else if (con_xup && con_zdown){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, -1, 0, +1, n);
		xF = x1N-size_box; 
		yF = y1N;
		zF = z1N+size_box;
		for(i=0; i<=lim_n1; ++i){
		xN = node[i][0][0].nodepos.x;
		for(j=0; j<partitions; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=partitions-1; k>=partitions-lim_n2; --k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, -1, 0, +1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, 0, +1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, 0, +1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<=lim_n1; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, 0, +1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 13
		//===================================================================================
		else if (con_xdown && con_zup){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, +1, 0, -1, n);
		xF = x1N+size_box; 
		yF = y1N;
		zF = z1N-size_box;
		for(i=partitions-1; i>=partitions-lim_n2; --i){
		xN = node[i][0][0].nodepos.x;
		for(j=0; j<partitions; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<=lim_n1; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, +1, 0, -1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, 0, -1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, 0, -1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i-1; u>=partitions-lim_n2; --u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, 0, -1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 14
		//===================================================================================
		else if (con_xdown && con_zdown){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, +1, 0, +1, n);
		xF = x1N+size_box; 
		yF = y1N;
		zF = z1N+size_box;
		for(i=partitions-1; i>=partitions-lim_n2; --i){
		xN = node[i][0][0].nodepos.x;
		for(j=0; j<partitions; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=partitions-1; k>=partitions-lim_n2; --k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, +1, 0, +1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, 0, +1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, 0, +1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i-1; u>=partitions-lim_n2; --u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<partitions; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, 0, +1);
			}
			}
			}}}
		}	
		} } }
		}
	}
	if(con_y && con_z){ // Nodos en esquinas YZ
		//===================================================================================
		// Boloque 15
		//===================================================================================
		if (con_yup && con_zup){
		len = node[row][col][mom].len;
		//Haacemos una copia en la esquina YZ contraria:
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, 0, -1, -1, n);
		xF = x1N;
		yF = y1N-size_box; 
		zF = z1N-size_box;
		for(i=0; i<partitions; ++i){
		xN = node[i][0][0].nodepos.x;
		for(j=0; j<=lim_n1; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<=lim_n1; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, 0, -1, -1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, -1, -1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, -1, -1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<partitions; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, -1, -1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 16
		//===================================================================================
		else if (con_yup && con_zdown){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, 0, -1, +1, n);
		xF = x1N;
		yF = y1N-size_box; 
		zF = z1N+size_box;
		for(i=0; i<partitions; ++i){
		xN = node[i][0][0].nodepos.x;
		for(j=0; j<=lim_n1; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=partitions-1; k>=partitions-lim_n2; --k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, 0, -1, +1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, -1, +1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, -1, +1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<partitions; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, -1, +1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 17
		//===================================================================================
		else if (con_ydown && con_zup){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, 0, +1, -1, n);
		xF = x1N;
		yF = y1N+size_box; 
		zF = z1N-size_box;
		for(i=0; i<partitions; ++i){
		xN = node[i][0][0].nodepos.x;
		for(j=partitions-1; j>=partitions-lim_n2; --j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<=lim_n1; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, 0, +1, -1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, +1, -1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, +1, -1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<partitions; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=partitions-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, +1, -1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 18
		//===================================================================================
		else if (con_ydown && con_zdown){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, 0, +1, +1, n);
		xF = x1N;
		yF = y1N+size_box; 
		zF = z1N+size_box;
		for(i=0; i<partitions; ++i){
		xN = node[i][0][0].nodepos.x;
		for(j=partitions-1; j>=partitions-lim_n2; --j){
		yN = node[i][j][0].nodepos.y;
		for(k=partitions-1; k>=partitions-lim_n2; --k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, 0, +1, +1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, +1, +1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, +1, +1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<partitions; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=partitions-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, 0, +1, +1);
			}
			}
			}}}
		}	
		} } }
		}
	}
	if(con_x && con_y && con_z){ // Nodos en esquinas XYZ
		//===================================================================================
		// Boloque 19
		//===================================================================================
		if (con_xup && con_yup && con_zup){
		len = node[row][col][mom].len;
		//Haacemos una copia en la esquina XYZ contraria:
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, -1, -1, -1, n);
		xF = x1N-size_box;
		yF = y1N-size_box;
		zF = z1N-size_box;
		for(i=0; i<=lim_n1; ++i){
		xN = node[i][0][0].nodepos.x;
		for(j=0; j<=lim_n1; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<=lim_n1; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, -1, -1, -1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, -1, -1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, -1, -1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<=lim_n1; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, -1, -1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 20
		//===================================================================================
		else if (con_xup && con_yup && con_zdown){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, -1, -1, +1, n);
		xF = x1N-size_box;
		yF = y1N-size_box;
		zF = z1N+size_box;
		for(i=0; i<=lim_n1; ++i){
		xN = node[i][0][0].nodepos.x;
		for(j=0; j<=lim_n1; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=partitions-1; k>=partitions-lim_n2; --k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, -1, -1, +1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, -1, +1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, -1, +1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<=lim_n1; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, -1, +1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 21
		//===================================================================================
		else if(con_xup && con_ydown && con_zup){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, -1, +1, -1, n);
		xF = x1N-size_box;
		yF = y1N+size_box;
		zF = z1N-size_box;
		for(i=0; i<=lim_n1; ++i){
		xN = node[i][0][0].nodepos.x;
		for(j=partitions-1; j>=partitions-lim_n2; --j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<=lim_n1; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, -1, +1, -1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, +1, -1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, +1, -1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<=lim_n1; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=partitions-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, +1, -1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 22
		//===================================================================================
		else if (con_xdown && con_yup && con_zup){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, +1, -1, -1, n);
		xF = x1N+size_box;
		yF = y1N-size_box;
		zF = z1N-size_box;
		for(i=partitions-1; i>=partitions-lim_n2; --i){
		xN = node[i][0][0].nodepos.x;
		for(j=0; j<=lim_n1; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<=lim_n1; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, +1, -1, -1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, -1, -1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, -1, -1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i-1; u>=partitions-lim_n2; --u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, -1, -1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 23
		//===================================================================================
		else if (con_xup && con_ydown && con_zdown){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, -1, +1, +1, n);
		xF = x1N-size_box;
		yF = y1N+size_box;
		zF = z1N+size_box;
		for(i=0; i<=lim_n1; ++i){
		xN = node[i][0][0].nodepos.x;
		for(j=partitions-1; j>=partitions-lim_n2; --j){
		yN = node[i][j][0].nodepos.y;
		for(k=partitions-1; k>=partitions-lim_n2; --k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, -1, +1, +1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, +1, +1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, +1, +1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i+1; u<=lim_n1; ++u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=partitions-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, -1, +1, +1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 24
		//===================================================================================
		else if (con_xdown && con_yup && con_zdown){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, +1, -1, +1, n);
		xF = x1N+size_box;
		yF = y1N-size_box;
		zF = z1N+size_box;
		for(i=partitions-1; i>=partitions-lim_n2; --i){
		xN = node[i][0][0].nodepos.x;
		for(j=0; j<=lim_n1; ++j){
		yN = node[i][j][0].nodepos.y;
		for(k=partitions-1; k>=partitions-lim_n2; --k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, +1, -1, +1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, -1, +1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j+1; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, -1, +1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i-1; u>=partitions-lim_n2; --u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=0; v<=lim_n1; ++v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, -1, +1);
			}
			}
			}}}
		}	
		} } }
		}
		//===================================================================================
		// Boloque 25
		//===================================================================================
		else if (con_xdown && con_ydown && con_zup){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, +1, +1, -1, n);
		xF = x1N+size_box;
		yF = y1N+size_box;
		zF = z1N-size_box;
		for(i=partitions-1; i>=partitions-lim_n2; --i){
		xN = node[i][0][0].nodepos.x;
		for(j=partitions-1; j>=partitions-lim_n2; --j){
		yN = node[i][j][0].nodepos.y;
		for(k=0; k<=lim_n1; ++k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, +1, +1, -1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k+1; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, +1, -1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, +1, -1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i-1; u>=partitions-lim_n2; --u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=partitions-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=0; w<=lim_n1; ++w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, +1, -1);
			}
			}
			}}}
		}	
		} } }
	 	}
		//===================================================================================
		// Boloque 26
		//===================================================================================
		else if (con_xdown && con_ydown && con_zdown){
		len = node[row][col][mom].len;
		add_front(nodeFront, node[row][col][mom].elements, len, x1N,  y1N,  z1N, +1, +1, +1, n);
		xF = x1N+size_box;
		yF = y1N+size_box;
		zF = z1N+size_box;
		for(i=partitions-1; i>=partitions-lim_n2; --i){
		xN = node[i][0][0].nodepos.x;
		for(j=partitions-1; j>=partitions-lim_n2; --j){
		yN = node[i][j][0].nodepos.y;
		for(k=partitions-1; k>=partitions-lim_n2; --k){
		zN = node[i][j][k].nodepos.z;	
		dx = xF - xN;
		dy = yF - yN;
		dz = zF - zN;
		dis_nod = dx*dx + dy*dy + dz*dz;
		if (dis_nod <= ddmax_nod){
			//======= caja/frontera ===============
			count_3_N112_front(row, col, mom, i, j, k, XXX, node, +1, +1, +1);
			u = i;
			v = j;
			//======= caja/caja/frontera ===============
			for (w=k-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, +1, +1);
			}
			}
			}
			u = i;
			//======= caja/caja/frontera ===============
			for (v=j-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){		
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, +1, +1);
			}
			}
			}}
			//======= caja/caja/frontera ===============
			for (u=i-1; u>=partitions-lim_n2; --u){
			xN_ = node[u][0][0].nodepos.x;
			dx = xN-xN_;
			for (v=partitions-1; v>=partitions-lim_n2; --v){
			yN_ = node[u][v][0].nodepos.y;
			dy = yN-yN_;
			for (w=partitions-1; w>=partitions-lim_n2; --w){
			zN_ = node[u][v][w].nodepos.z;
			dz = zN-zN_;
			dis_nod = dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
			dx = xF-xN_;
			dy = yF-yN_;
			dz = zF-zN_;
			dis_nod = dx*dx + dy*dy + dz*dz;
			if (dis_nod <= ddmax_nod){
				count_3_N123_front(row, col, mom, i, j, k, u, v, w, XXX, node, +1, +1, +1);
			}
			}
			}}}
		}	
		} } }
		}
	}
	
	}
	}
	}

}

//=================================================================== 
void NODE3P::count_3_N112_front(int row, int col, int mom, int u, int v, int w, unsigned int ***XXX, Node ***nodeS, short int sig1, short int sig2, short int sig3){
	/*
	Funcion para contar los triangulos en dos 
	nodos con dos puntos en N1 y un punto en N2.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	
	*/
	int i,j,k;
	float dx,dy,dz;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	int a, b, c;
	for (i=0; i<nodeS[u][v][w].len; ++i){
		// 1er punto en N2 (caja)
		x1 = nodeS[u][v][w].elements[i].x;
		y1 = nodeS[u][v][w].elements[i].y;
		z1 = nodeS[u][v][w].elements[i].z;
		for (j=0; j<nodeS[row][col][mom].len; ++j){
			// 2do punto en N1 (frontera)
			x2 = nodeS[row][col][mom].elements[j].x+(size_box*(sig1));
			y2 = nodeS[row][col][mom].elements[j].y+(size_box*(sig2));
			z2 = nodeS[row][col][mom].elements[j].z+(size_box*(sig3));
			dx = x2-x1;
			dy = y2-y1;
			dz = z2-z1;
			d12 = dx*dx+dy*dy+dz*dz;
			if (d12<=dd_max){
			for (k=j+1; k<nodeS[row][col][mom].len; ++k){
				// 3er punto en N1 (frontera)
				x3 = nodeS[row][col][mom].elements[k].x+(size_box*(sig1));
				y3 = nodeS[row][col][mom].elements[k].y+(size_box*(sig2));
				z3 = nodeS[row][col][mom].elements[k].z+(size_box*(sig3));
				dx = x3-x1;
				dy = y3-y1;
				dz = z3-z1;
				d13 = dx*dx+dy*dy+dz*dz;
				if (d13<=dd_max){
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23<=dd_max){
				a = (int)(sqrt(d12)*ds);
				b = (int)(sqrt(d13)*ds);
				c = (int)(sqrt(d23)*ds);
				*(*(*(XXX+a)+b)+c)+=1;
				*(*(*(XXX+b)+a)+c)+=1;
				}
				}
			}
			for (k=i+1; k<nodeS[u][v][w].len; ++k){
				// 3er punto en N2 
				x3 = nodeS[u][v][w].elements[k].x;
				y3 = nodeS[u][v][w].elements[k].y;
				z3 = nodeS[u][v][w].elements[k].z;
				dx = x3-x1;
				dy = y3-y1;
				dz = z3-z1;
				d13 = dx*dx+dy*dy+dz*dz;
				if (d13<=dd_max){
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23<=dd_max){
				a = (int)(sqrt(d12)*ds);
				b = (int)(sqrt(d13)*ds);
				c = (int)(sqrt(d23)*ds);
				*(*(*(XXX+a)+b)+c)+=1;
				*(*(*(XXX+c)+b)+a)+=1;
				*(*(*(XXX+b)+a)+c)+=1;
				*(*(*(XXX+b)+c)+a)+=1;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N123_front(int row, int col, int mom, int a, int b, int c, int u, int v, int w, unsigned int ***XXX, Node ***nodeS, short int sig1, short int sig2, short int sig3){
	/*
	Funcion para contar los triangulos en tres 
	nodos con un puntos en N1, un punto en N2
	y un punto en N3.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	a, b, c => posición de N3.
	
	*/
	int i,j,k, n, m, l;
	float dx,dy,dz;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	for (i=0; i<nodeS[row][col][mom].len; ++i){
		// 1er punto en N1 (frontera)
		x1 = nodeS[row][col][mom].elements[i].x+(size_box*(sig1));
		y1 = nodeS[row][col][mom].elements[i].y+(size_box*(sig2));
		z1 = nodeS[row][col][mom].elements[i].z+(size_box*(sig3));
		for (j=0; j<nodeS[a][b][c].len; ++j){
			// 2do punto en N3 (caja)
			x2 = nodeS[a][b][c].elements[j].x;
			y2 = nodeS[a][b][c].elements[j].y;
			z2 = nodeS[a][b][c].elements[j].z;
			dx = x2-x1;
			dy = y2-y1;
			dz = z2-z1;
			d12 = dx*dx+dy*dy+dz*dz;
			if (d12<=dd_max){
			for (k=0; k<nodeS[u][v][w].len; ++k){
				// 3er punto en N2 (caja)
				x3 = nodeS[u][v][w].elements[k].x;
				y3 = nodeS[u][v][w].elements[k].y;
				z3 = nodeS[u][v][w].elements[k].z;
				dx = x3-x1;
				dy = y3-y1;
				dz = z3-z1;
				d13 = dx*dx+dy*dy+dz*dz;
				if (d13<=dd_max){
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23<=dd_max){
				n = (int)(sqrt(d12)*ds);
				m = (int)(sqrt(d13)*ds);
				l = (int)(sqrt(d23)*ds);
				*(*(*(XXX+n)+m)+l)+=1;
				*(*(*(XXX+m)+n)+l)+=1;
				*(*(*(XXX+l)+m)+n)+=1;
				*(*(*(XXX+l)+n)+m)+=1;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N123_front_FFN(int n, int m, int u, int v, int w, unsigned int ***XXX, Node ***nodeS, Node *nodeF){
	/*
	Funcion para contar los triangulos en tres 
	nodos con un puntos en N1, un punto en N2
	y un punto en N3.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	a, b, c => posición de N3.
	
	*/
	int i, j, k, a, b, c;
	float dx,dy,dz;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	//std::cout << nodeF[n].len << std::endl;
	for (i=0; i<nodeF[n].len; ++i){
		// 1er punto en N1 (frontera)
		x1 = nodeF[n].elements[i].x;
		y1 = nodeF[n].elements[i].y;
		z1 = nodeF[n].elements[i].z;
		for (j=0; j<nodeF[m].len; ++j){
			// 2do punto en N3 (frontera)
			x2 = nodeF[m].elements[j].x;
			y2 = nodeF[m].elements[j].y;
			z2 = nodeF[m].elements[j].z;
			dx = x2-x1;
			dy = y2-y1;
			dz = z2-z1;
			d12 = dx*dx+dy*dy+dz*dz;
			if (d12<=dd_max){
			for (k=0; k<nodeS[u][v][w].len; ++k){
				// 3er punto en N2 (caja)
				x3 = nodeS[u][v][w].elements[k].x;
				y3 = nodeS[u][v][w].elements[k].y;
				z3 = nodeS[u][v][w].elements[k].z;
				dx = x3-x1;
				dy = y3-y1;
				dz = z3-z1;
				d13 = dx*dx+dy*dy+dz*dz;
				if (d13<=dd_max){
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23<=dd_max){
				a = (int)(sqrt(d12)*ds);
				b = (int)(sqrt(d13)*ds);
				c = (int)(sqrt(d23)*ds);
				//std::cout << "-..." << std::endl;
				*(*(*(XXX+a)+b)+c)+=1;
				*(*(*(XXX+b)+a)+c)+=1;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
NODE3P::~NODE3P(){	
}

