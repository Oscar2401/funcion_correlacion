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

struct conditions{
	float disx1;
	float disx2;
	float disy1;
	float disy2;
	float disz1;
	float disz2;
	float dx1;
	float dx2;
	float dy1;
	float dy2;
	float dz1;
	float dz2; 
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
		void count_3_N111(int, int, int, unsigned int ***, Node ***);
		void count_3_N112(int, int, int, int, int, int, unsigned int ***, Node ***);
		void count_3_N122(int, int, int, int, int, int, unsigned int ***, Node ***);
		void count_3_N123(int, int, int, int, int, int, int, int, int, unsigned int ***, Node ***);
		void symmetrize(unsigned int ***);
		//void histo_front_XX(unsigned int *, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int);
		//void histo_front_XY(unsigned int *, Node ***, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int);
		~NODE3P();
};

//=================================================================== 
//==================== Funciones ==================================== 
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
	
	float dis_f;
	bool con_x, con_y, con_z, con_x2, con_y2, con_z2;
	conditions cond;
	
	// x1N, y1N, z1N => Nodo pivote
	for (row=0; row<partitions; ++row){
	x1N = nodeX[row][0][0].nodepos.x;
	for (col=0; col<partitions; ++col){
	y1N = nodeX[row][col][0].nodepos.y;
	for (mom=0; mom<partitions; ++mom){
	z1N = nodeX[row][col][mom].nodepos.z;
	
	u = row;
	v = col;
	//=======================
	// Nodo 2 movil en Z:
	//=======================
	if(z1N<=d_max_pm){.
	x2N = nodeX[u][0][0].nodepos.x;
	y2N = nodeX[u][v][0].nodepos.y;
	for (w=mom+1;  w<partitions; ++w){	
	z2N = nodeX[u][v][w].nodepos.z;
		
		con_z = ((z1N<=d_max_pm)&&(z2N>=front_pm))||((z2N<=d_max_pm)&&(z1N>=front_pm));
		if(con_z){
			dz_nod = fabs(z2N-z1N);
			front_NN2(XXX, nodeX, 0.0, 0.0, dn_z, false, false, con_z, row, col, mom, u, v, w);
		}
		a = u;
		b = v;
		//=======================
		// Nodo 3 movil en Z:
		//=======================
		for (c=w+1;  c<partitions; ++c){
		z3N = nodeX[a][b][c].nodepos.z;
		
		con_z21 = ((z1N<=d_max_pm)&&(z2N<=d_max_pm)&&(z3N>=front_pm))||((z1N>=front_pm)&&(z2N>=front_pm)&&(z3N<=d_max_pm));
		con_z22 = ((z1N<=d_max_pm)&&(z2N>=front_pm)&&(z3N>=front_pm))||((z1N>=front_pm)&&(z2N<=d_max_pm)&&(z3N<=d_max_pm));
		con_z23 = ((z1N<=d_max_pm)&&(z2N>=front_pm)&&(z3N<=d_max_pm))||((z1N>=front_pm)&&(z2N<=d_max_pm)&&(z3N>=front_pm));
		
		con_z2 = con_z21||con_z22||con_z23;
		
		if (con_z2){
			dz31 = z3N - z1N;
			disN31 = dz31*dz31;
			dz32 = z3N - z2N;
			disN32 = dz32*dz32;
			dz21 = z2N - z1N;
			disN21 = dz21*dz21;
			if (con_z21){
				cond.dz1 =  fabs(dz31);
				cond.dz2 =  fabs(dz32);
				cond.disz1 = disN31;
				cond.disz2 = disN32;
				
			}
			else if (con_z22){
				cond.dz1 =  fabs(dz21);
				cond.dz2 =  fabs(dz31);
				cond.disz1 = disN21;
				cond.disz2 = disN31;
			}
			else if (con_z23){
				cond.dz1 =  fabs(dz21);
				cond.dz2 =  fabs(dz32);
				cond.disz1 = disN21;
				cond.disz2 = disN32;
			}
		}
		
		}
		//=======================
		// Nodo 3 movil en ZY:
		//=======================
		for (b=v+1; b<partitions; ++b){
		y3N = nodeX[a][b][0].nodepos.y;
		for (c=0;  c<partitions; ++c){
		z3N = nodeX[a][b][c].nodepos.z;
			
		con_y21 = ((y1N<=d_max_pm)&&(y2N<=d_max_pm)&&(y3N>=front_pm))||((y1N>=front_pm)&&(y2N>=front_pm)&&(y3N<=d_max_pm));
		con_y22 = ((y1N<=d_max_pm)&&(y2N>=front_pm)&&(y3N>=front_pm))||((y1N>=front_pm)&&(y2N<=d_max_pm)&&(y3N<=d_max_pm));
		con_y23 = ((y1N<=d_max_pm)&&(y2N>=front_pm)&&(y3N<=d_max_pm))||((y1N>=front_pm)&&(y2N<=d_max_pm)&&(y3N>=front_pm));
		
		con_z21 = ((z1N<=d_max_pm)&&(z2N<=d_max_pm)&&(z3N>=front_pm))||((z1N>=front_pm)&&(z2N>=front_pm)&&(z3N<=d_max_pm));
		con_z22 = ((z1N<=d_max_pm)&&(z2N>=front_pm)&&(z3N>=front_pm))||((z1N>=front_pm)&&(z2N<=d_max_pm)&&(z3N<=d_max_pm));
		con_z23 = ((z1N<=d_max_pm)&&(z2N>=front_pm)&&(z3N<=d_max_pm))||((z1N>=front_pm)&&(z2N<=d_max_pm)&&(z3N>=front_pm));
		
		con_y2 = con_y21||con_y22||con_y23;
		con_z2 = con_z21||con_z22||con_z23;
		
		if (con_y2 || con_z2){
			dy31 = y3N - y1N;
			dz31 = z3N - z1N;
			disN31 = dy31*dy31 + dz31*dz31;
			dy32 = y3N - y2N;
			dz32 = z3N - z2N;
			disN32 = dy32*dy32 + dz32*dz32;
			dy21 = y2N - y1N;
			dz21 = z2N - z1N;
			disN21 = dy21*dy21 + dz21*dz21;
			if (con_y21){
				if (con_z21){
					if (disN21 <= ddmax_nod){
						cond.dy1 =  fabs(dy31);
						cond.dy2 =  fabs(dy32);
						cond.disy1 = disN31;
						cond.disy2 = disN32;
						cond.dz1 =  fabs(dz31);
						cond.dz2 =  fabs(dz32);
						cond.disz1 = disN31;
						cond.disz2 = disN32;
					
					}
				}
				else if (con_z22){
					cond.dy1 =  fabs(dy31);
					cond.dy2 =  fabs(dy32);
					cond.disy1 = disN31;
					cond.disy2 = disN32;
					cond.dz1 =  fabs(dz21);
					cond.dz2 =  fabs(dz31);
					cond.disz1 = disN21;
					cond.disz2 = disN31;
					
				}
				else if (con_z23){
					cond.dy1 =  fabs(dy31);
					cond.dy2 =  fabs(dy32);
					cond.disy1 = disN31;
					cond.disy2 = disN32;
					cond.dz1 =  fabs(dz21);
					cond.dz2 =  fabs(dz32);
					cond.disz1 = disN21;
					cond.disz2 = disN32;
					
				}
				else{
					
				}
			}
			else if (con_y22){
				if (con_z21){
					cond.dy1 =  fabs(dy21);
					cond.dy2 =  fabs(dy31);
					cond.disy1 = disN21;
					cond.disy2 = disN31;
					cond.dz1 =  fabs(dz31);
					cond.dz2 =  fabs(dz32);
					cond.disz1 = disN31;
					cond.disz2 = disN32;
					
				}
				else if (con_z22){
					cond.dy1 =  fabs(dy21);
					cond.dy2 =  fabs(dy31);
					cond.disy1 = disN21;
					cond.disy2 = disN31;
					cond.dz1 =  fabs(dz21);
					cond.dz2 =  fabs(dz31);
					cond.disz1 = disN21;
					cond.disz2 = disN31;
					
				}
				else if (con_z23){
					cond.dy1 =  fabs(dy21);
					cond.dy2 =  fabs(dy31);
					cond.disy1 = disN21;
					cond.disy2 = disN31;
					cond.dz1 =  fabs(dz21);
					cond.dz2 =  fabs(dz32);
					cond.disz1 = disN21;
					cond.disz2 = disN32;
					
				}
				else{
					
				}
			}
			else if (con_y23){
				if (con_z21){
					cond.dy1 =  fabs(dy21);
					cond.dy2 =  fabs(dy32);
					cond.disy1 = disN21;
					cond.disy2 = disN32;
					cond.dz1 =  fabs(dz31);
					cond.dz2 =  fabs(dz32);
					cond.disz1 = disN31;
					cond.disz2 = disN32;
					
				}
				else if (con_z22){
					cond.dy1 =  fabs(dy21);
					cond.dy2 =  fabs(dy32);
					cond.disy1 = disN21;
					cond.disy2 = disN32;
					cond.dz1 =  fabs(dz21);
					cond.dz2 =  fabs(dz31);
					cond.disz1 = disN21;
					cond.disz2 = disN31;
					
				}
				else if (con_z23){
					cond.dy1 =  fabs(dy21);
					cond.dy2 =  fabs(dy32);
					cond.disy1 = disN21;
					cond.disy2 = disN32;
					cond.dz1 =  fabs(dz21);
					cond.dz2 =  fabs(dz32);
					cond.disz1 = disN21;
					cond.disz2 = disN32;
					
				}
				else{
					
				}
			}
			else{
				if (con_z21){
				
				}
				else if (con_z22){
			
				}
				else if (con_z23){
			
				}
				else{
					
				}
			}
		}
		
		}
		}
		//=======================
		// Nodo 3 movil en ZYX:
		//=======================
		for (a=u+1; a<partitions; ++a){
		x3N = nodeX[a][0][0].nodepos.x;
		for (b=0; b<partitions; ++b){
		y3N = nodeX[a][b][0].nodepos.y;
		for (c=0;  c<partitions; ++c){
		z3N = nodeX[a][b][c].nodepos.z;
		
		con_x2 = ((x1N<=d_max_pm)&&(x2N<=d_max_pm)&&(x3N>=front_pm))||((x1N>=front_pm)&&(x2N>=front_pm)&&(x3N<=d_max_pm));
		con_y2 = ((y1N<=d_max_pm)&&(y2N<=d_max_pm)&&(y3N>=front_pm))||((y1N>=front_pm)&&(y2N>=front_pm)&&(y3N<=d_max_pm));
		con_z2 = ((z1N<=d_max_pm)&&(z2N<=d_max_pm)&&(z3N>=front_pm))||((z1N>=front_pm)&&(z2N>=front_pm)&&(z3N<=d_max_pm));
		if(con_x2||con_y2||con_z2){
			
		}
		
		con_x2 = ((x1N<=d_max_pm)&&(x2N>=front_pm)&&(x3N>=front_pm))||((x1N>=front_pm)&&(x2N<=d_max_pm)&&(x3N<=d_max_pm));
		con_y2 = ((y1N<=d_max_pm)&&(y2N>=front_pm)&&(y3N>=front_pm))||((y1N>=front_pm)&&(y2N<=d_max_pm)&&(y3N<=d_max_pm));
		con_z2 = ((z1N<=d_max_pm)&&(z2N>=front_pm)&&(z3N>=front_pm))||((z1N>=front_pm)&&(z2N<=d_max_pm)&&(z3N<=d_max_pm));
		if(con_x2||con_y2||con_z2){
			
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
	for (w=0; w<partitions ; ++w){		
	z2N = nodeX[u][v][w].nodepos.z;
		
		con_y = ((y1N<=d_max_pm)&&(y2N>=front_pm))||((y2N<=d_max_pm)&&(y1N>=front_pm));
		con_z = ((z1N<=d_max_pm)&&(z2N>=front_pm))||((z2N<=d_max_pm)&&(z1N>=front_pm));
		if(con_y||con_z){
		dy_nod = fabs(y2N-y1N);
		dz_nod = fabs(z2N-z1N);
		front_NN2(XXX, nodeX, 0.0, dn_y, dn_z, false, con_y, con_z, row, col, mom, u, v, w);
		}
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
	for (v=0; v<partitions; ++v){
	y2N = nodeX[u][v][0].nodepos.y;
	for (w=0; w<partitions; ++w){
	z2N = nodeX[u][v][w].nodepos.z;
		con_x = ((x1N<=d_max_pm)&&(x2N>=front_pm))||((x2N<=d_max_pm)&&(x1N>=front_pm));
		con_y = ((y1N<=d_max_pm)&&(y2N>=front_pm))||((y2N<=d_max_pm)&&(y1N>=front_pm));
		con_z = ((z1N<=d_max_pm)&&(z2N>=front_pm))||((z2N<=d_max_pm)&&(z1N>=front_pm));
		if(con_x||con_y||con_z){
			dx_nod = fabs(x2N-x1N);
			dy_nod = fabs(y2N-y1N);
			dz_nod = fabs(z2N-z1N);
			front_NN2(XXX, nodeX, dx_nod, dn_y, dn_z, con_x, con_y, con_z, row, col, mom, u, v, w);
		}	
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
void NODE3P::front_NN2(unsigned int *XXX, Node ***nodeX, float disn, float dn_x, float dn_y, float dn_z, bool con_x, bool con_y, bool con_z, int row, int col, int mom, int u, int v, int w){
	int i, j;
	float dis_f;
	if( con_x ){
		dis_f = disn + ll - 2*dn_x*size_box;
		if (dis_f <= ddmax_nod){
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, 1, 0, 0);
		}
	}
	if( con_y ){
		dis_f = disn + ll - 2*dn_y*size_box;
		if (dis_f <= ddmax_nod){
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, 0, 1, 0);
		}
	}
	if( con_z ){
		dis_f = disn + ll - 2*dn_z*size_box;
		if (dis_f <= ddmax_nod){
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, 0, 0, 1);
		}
	}		
	if( con_x && con_y ){
		dis_f = disn + 2*ll - 2*(dn_x+dn_y)*size_box;
		if (dis_f < ddmax_nod){
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, 1, 1, 0);
		}
	}				
	if( con_in_x && con_in_z ){
		dis_f = disn + 2*ll - 2*(dn_x+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, 1, 0, 1);
		}
	}			
	if( con_in_y && con_in_z ){
		dis_f = disn + 2*ll - 2*(dn_y+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, 0, 1, 1);
		}
	}		
	if( con_in_x && con_in_y && con_in_z ){
		dis_f = disn + 3*ll - 2*(dn_x+dn_y+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, 1, 1, 1);
		}
	}
}
//=================================================================== 
void NODE3P::front_NNN3(unsigned int *XXX, Node ***nodeX, conditions cond, ){
	int i, j;
	float dis_f;
	if( con_x ){
		dis_f = cond.disx1 + ll - 2*cond.dx1*size_box;
		if (dis_f <= ddmax_nod){
			dis_f = cond.disx2 + ll - 2*cond.dx2*size_box;
			if (dis_f <= ddmax_nod){
				count_3_N123_front(row, col, mom, u, v, w, XXX, nodeX, 1, 0, 0);
			}
		}
	}
	if( con_y ){
		dis_f = cond.disy1 + ll - 2*cond.dy1*size_box;
		if (dis_f <= ddmax_nod){
			dis_f = cond.disy2 + ll - 2*cond.dy2*size_box;
			if (dis_f <= ddmax_nod){
				count_3_N123_front(row, col, mom, u, v, w, XXX, nodeX, 1, 0, 0);
			}
		}
	}
	if( con_z ){
		dis_f = disn1 + ll - 2*dn_z1*size_box;
		if (dis_f <= ddmax_nod){
			dis_f = disn2 + ll - 2*dn_z2*size_box;
			if (dis_f <= ddmax_nod){
				count_3_N123_front(row, col, mom, u, v, w, XXX, nodeX, 0, 0, 1);
			}
		}
	}		
	if( con_x && con_y ){
		dis_f = disn + 2*ll - 2*(dn_x+dn_y)*size_box;
		if (dis_f < ddmax_nod){
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, 1, 1, 0);
		}
	}				
	if( con_in_x && con_in_z ){
		dis_f = disn + 2*ll - 2*(dn_x+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, 1, 0, 1);
		}
	}			
	if( con_in_y && con_in_z ){
		dis_f = disn + 2*ll - 2*(dn_y+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, 0, 1, 1);
		}
	}		
	if( con_in_x && con_in_y && con_in_z ){
		dis_f = disn + 3*ll - 2*(dn_x+dn_y+dn_z)*size_box;
		if (dis_f <= ddmax_nod){
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, 1, 1, 1);
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N112_front(int row, int col, int mom, int u, int v, int w, unsigned int ***XXX, Node ***nodeS, short int Mx, short int My, short int Mz){
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
			dx = fabs(x2-x1)-(size_box*Mx);
			dy = fabs(y2-y1)-(size_box*My);
			dz = fabs(z2-z1)-(size_box*Mz);
			d12 = dx*dx+dy*dy+dz*dz;
			if (d12<=dd_max){
			for (k=j+1; k<nodeS[row][col][mom].len; ++k){
				// 3er punto en N1
				x3 = nodeS[row][col][mom].elements[k].x;
				y3 = nodeS[row][col][mom].elements[k].y;
				z3 = nodeS[row][col][mom].elements[k].z;
				dx = fabs(x3-x1)-(size_box*Mx);
				dy = fabs(y3-y1)-(size_box*My);
				dz = fabs(z3-z1)-(size_box*Mz);
				d13 = dx*dx+dy*dy+dz*dz;
				if (d13<=dd_max){
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23<=dd_max){
				*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=2;
				}
				}
			}
			for (k=i+1; k<nodeS[u][v][w].len; ++k){
				// 3er punto en N2
				x3 = nodeS[u][v][w].elements[k].x;
				y3 = nodeS[u][v][w].elements[k].y;
				z3 = nodeS[u][v][w].elements[k].z;
				dx = fabs(x3-x1)-(size_box*Mx);
				dy = fabs(y3-y1)-(size_box*My);
				dz = fabs(z3-z1)-(size_box*Mz);
				d13 = dx*dx+dy*dy+dz*dz;
				if (d13<=dd_max){
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23<=dd_max){
				*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=2;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N123_front(int row, int col, int mom, int u, int v, int w, int a, int b, int c, unsigned int ***XXX, Node ***nodeS, short int Mx, short int My, short int Mz, short int n1, short int n2){
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
			dx = fabs(x3-x1)-(size_box*Mx);
			dy = fabs(y3-y1)-(size_box*My);
			dz = fabs(z3-z1)-(size_box*Mz);
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<=dd_max){
			for (k=0; k<nodeS[u][v][w].len; ++k){
				// 3er punto en N2
				x2 = nodeS[u][v][w].elements[k].x;
				y2 = nodeS[u][v][w].elements[k].y;
				z2 = nodeS[u][v][w].elements[k].z;
				dx = fabs(x3-x2)-(size_box*Mx)*n1; //cuando N3 y N2 están en proyectados n1 = 0 y n2 = 1 
				dy = fabs(y3-y2)-(size_box*My)*n1;
				dz = fabs(z3-z2)-(size_box*Mz)*n1;
				d23 = dx*dx+dy*dy+dz*dz;
				if (d23<=dd_max){
				dx = fabs(x2-x1)-(size_box*Mx)*n2; //cuando solo N3 están en proyectados n1 = 1 y n2 = 0 
				dy = fabs(y2-y1)-(size_box*My)*n2;
				dz = fabs(z2-z1)-(size_box*Mz)*n2;
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
NODE3P::~NODE3P(){	
}

