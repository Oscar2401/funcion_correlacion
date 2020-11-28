#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <omp.h>

struct Point3D{
	float x;
	float y; 
	float z;
};


struct PointW3D{
	float x;
	float y; 
	float z;
	float w;
};

struct Node{
	Point3D nodepos;	// Coordenadas del nodo (posición del nodo).
	int len;		// Cantidad de elementos en el nodo.
	PointW3D *elements;	// Elementos del nodo.
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
		int partitions;
		float size_box;
		float size_node;
		float d_max;
		Node ***nodeD;
		PointW3D *dataD;
		Node ***nodeR;
		PointW3D *dataR;
		// Derivados
		float ll;
		float dd_max;
		float corr;
		float front;
		float ds;
		float ds_th;
		float ddmax_nod;
		
	private: 
		void make_nodos(Node ***, PointW3D *);
		void add(PointW3D *&, int&, float, float, float, float);
	
	// Métodos de Clase:
	public:
		//Constructor de clase:
		NODE3P(int _bn, int _n_pts, float _size_box, float _size_node, float _d_max, PointW3D *_dataD, Node ***_nodeD, PointW3D *_dataR, Node ***_nodeR){
			bn = _bn;
			n_pts = _n_pts;
			size_box = _size_box;
			size_node = _size_node;
			d_max = _d_max;
			dataD = _dataD;
			nodeD = _nodeD;
			dataR = _dataR;
			nodeR = _nodeR;
			dd_max = d_max*d_max;
			front = size_box - d_max;
			corr = size_node*sqrt(3);
			ds = ((float)(bn))/d_max;
			ds_th = (float)(bn)/2;
			ddmax_nod = (d_max+corr)*(d_max+corr);
			partitions = (int)(ceil(size_box/size_node));
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
		void make_histoXXX(double *****, Node ***);
		void count_3_N111(int, int, int, double *****, Node ***);
		void count_3_N112(int, int, int, int, int, int, double *****, Node ***);
		void count_3_N123(int, int, int, int, int, int, int, int, int, double *****, Node ***);
		
		void make_histoXXY(double *****, Node ***, Node ***);
		void count_3_N112_xxy(int, int, int, int, int, int, double *****, Node ***, Node ***);
		void count_3_N123_xxy(int, int, int, int, int, int, int, int, int, double *****, Node ***, Node ***);
		~NODE3P();
};

//=================================================================== 
//===================== Funciones =================================== 
//=================================================================== 

void NODE3P::make_nodos(Node ***nod, PointW3D *dat){
	/*
	Función para crear los nodos con los datos y puntos random
	
	Argumentos
	nod: arreglo donde se crean los nodos.
	dat: datos a dividir en nodos.
	
	*/
	int i, row, col, mom;
	float p_med = size_node/2;
	
	// Inicializamos los nodos vacíos:
	for (row=0; row<partitions; row++){
		for (col=0; col<partitions; col++){
			for (mom=0; mom<partitions; mom++){
				nod[row][col][mom].nodepos.x = ((float)(row)*(size_node))+p_med;
				nod[row][col][mom].nodepos.y = ((float)(col)*(size_node))+p_med;
				nod[row][col][mom].nodepos.z = ((float)(mom)*(size_node))+p_med;
				nod[row][col][mom].len = 0;
				nod[row][col][mom].elements = new PointW3D[0];
			}
		}
	}
	// Llenamos los nodos con los puntos de dat:
	for (i=0; i<n_pts; i++){
		row = (int)(dat[i].x/size_node);
        	col = (int)(dat[i].y/size_node);
        	mom = (int)(dat[i].z/size_node);
		add(nod[row][col][mom].elements, nod[row][col][mom].len, dat[i].x, dat[i].y, dat[i].z, dat[i].w);
	}
}
//=================================================================== 
void NODE3P::add(PointW3D *&array, int &lon, float _x, float _y, float _z, float _w){
	lon++;
	PointW3D *array_aux = new PointW3D[lon];
	for (int i=0; i<lon-1; i++){
		array_aux[i].x = array[i].x;
		array_aux[i].y = array[i].y;
		array_aux[i].z = array[i].z;
		array_aux[i].w = array[i].w;
	}
	delete[] array;
	array = array_aux;
	array[lon-1].x = _x;
	array[lon-1].y = _y; 
	array[lon-1].z = _z;
	array[lon-1].w = _w; 
}
//=================================================================== 
void NODE3P::make_histoXXX(double *****XXX, Node ***nodeX){
	/*
	Función para crear los histogramas DDD y RRR.
	
	Argumentos
	DDD: arreglo donde se creará el histograma DDD.
	*/ 
	int i, j, k, row, col, mom, u, v, w, a ,b, c;
	float dis, dis_nod, dis_nod2, dis_nod3;
	float x1N, y1N, z1N, x2N, y2N, z2N, x3N, y3N, z3N;
	float x, y, z;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod, dx_nod2, dy_nod2, dz_nod2, dx_nod3, dy_nod3, dz_nod3;
	bool con_x, con_y, con_z;
	
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
			dy_nod2 *= dy_nod2;
			for (c=0;  c<partitions; ++c){
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod2 = z3N-z1N;
				dz_nod2 *= dz_nod2;
				dis_nod2 = dy_nod2 + dz_nod2;
				if (dis_nod2 <= ddmax_nod){
				dy_nod3 = y3N-y2N;
				dy_nod3 *= dy_nod3;
				dz_nod3 = z3N-z2N;
				dz_nod3 *= dz_nod3;
				dis_nod3 = dy_nod3 + dz_nod3;
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
			dx_nod2 *= dx_nod2;
			for (b=0; b<partitions; ++b){
				y3N = nodeX[a][b][0].nodepos.y;
				dy_nod2 = y3N-y1N;
				dy_nod2 *= dy_nod2;
				for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					dz_nod2 *= dz_nod2;
					dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
					if (dis_nod2 <= ddmax_nod){
					dx_nod3 = x3N-x2N;
					dx_nod3 *= dx_nod3;
					dy_nod3 = y3N-y2N;
					dy_nod3 *= dy_nod3;
					dz_nod3 = z3N-z2N;
					dz_nod3 *= dz_nod3;
					dis_nod3 = dx_nod3 + dy_nod3 + dz_nod3;
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
	for (v=col+1; v<partitions ; ++v){
		y2N = nodeX[u][v][0].nodepos.y;
		dy_nod = y2N-y1N;
		dy_nod *= dy_nod;
		for (w=0; w<partitions ; ++w){		
			z2N = nodeX[u][v][w].nodepos.z;
			dz_nod = z2N-z1N;
			dz_nod *= dz_nod;
			dis_nod = dy_nod + dz_nod;
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
			dy_nod2 *= dy_nod2;
			for (c=w+1;  c<partitions; ++c){
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod2 = z3N-z1N;
				dz_nod2 *= dz_nod2;
				dis_nod2 = dy_nod2 + dz_nod2;
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
				dy_nod2 *= dy_nod2;
				for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					dz_nod2 *= dz_nod2;
					dis_nod2 = dy_nod2 + dz_nod2;
					if (dis_nod2 <= ddmax_nod){
					dy_nod3 = y3N-y2N;
					dy_nod3 *= dy_nod3;
					dz_nod3 = z3N-z2N;
					dz_nod3 *= dz_nod3;
					dis_nod3 = dy_nod3 + dz_nod3;
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
				dx_nod2 *= dx_nod2;
				for (b=0; b<partitions; ++b){
					y3N = nodeX[a][b][0].nodepos.y;
					dy_nod2 = y3N-y1N;
					dy_nod2 *= dy_nod2;
					for (c=0;  c<partitions; ++c){
						z3N = nodeX[a][b][c].nodepos.z;
						dz_nod2 = z3N-z1N;
						dz_nod2 *= dz_nod2;
						dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
						if (dis_nod2 <= ddmax_nod){
						dx_nod3 = x3N-x2N;
						dx_nod3 *= dx_nod3;
						dy_nod3 = y3N-y2N;
						dy_nod3 *= dy_nod3;
						dz_nod3 = z3N-z2N;
						dz_nod3 *= dz_nod3;
						dis_nod3 = dx_nod3 + dy_nod3 + dz_nod3;
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
		dx_nod *= dx_nod;
		for (v=0; v<partitions; ++v){
			y2N = nodeX[u][v][0].nodepos.y;
			dy_nod = y2N-y1N;
			dy_nod *= dy_nod;
			for (w=0; w<partitions; ++w){
				z2N = nodeX[u][v][w].nodepos.z;
				dz_nod = z2N-z1N;
				dz_nod *= dz_nod;
				dis_nod = dx_nod + dy_nod + dz_nod;
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
				dx_nod2 *= dx_nod2;
				dy_nod2 = y3N-y1N;
				dy_nod2 *= dy_nod2;
				for (c=w+1;  c<partitions; ++c){	
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					dz_nod2 *= dz_nod2;
					dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
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
					dy_nod2 *= dy_nod2;
					for (c=0;  c<partitions; ++c){
						z3N = nodeX[a][b][c].nodepos.z;
						dz_nod2 = z3N-z1N;
						dz_nod2 *= dz_nod2;
						dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
						if (dis_nod2 <= ddmax_nod){
						dy_nod3 = y3N-y2N;
						dy_nod3 *= dy_nod3;
						dz_nod3 = z3N-z2N;
						dz_nod3 *= dz_nod3;
						dis_nod3 = dy_nod3 + dz_nod3;
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
					dx_nod2 *= dx_nod2;
					for (b=0; b<partitions; ++b){
						y3N = nodeX[a][b][0].nodepos.y;
						dy_nod2 = y3N-y1N;
						dy_nod2 *= dy_nod2;
						for (c=0;  c<partitions; ++c){
							z3N = nodeX[a][b][c].nodepos.z;
							dz_nod2 = z3N-z1N;
							dz_nod2 *= dz_nod2;
							dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
							if (dis_nod2 <= ddmax_nod){
							dx_nod3 = x3N-x2N;
							dx_nod3 *= dx_nod3;
							dy_nod3 = y3N-y2N;
							dy_nod3 *= dy_nod3;
							dz_nod3 = z3N-z2N;
							dz_nod3 *= dz_nod3;
							dis_nod3 = dx_nod3 + dy_nod3 + dz_nod3;
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
}
//=================================================================== 
void NODE3P::count_3_N111(int row, int col, int mom, double *****XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en un mismo Nodo.
	
	row, col, mom => posición del Nodo. Esto define al Nodo.
	
	*/
	int i,j,k;
	int a_,b_,c_;
	int t,p,q, t_,p_,q_;
	float dx,dy,dz1,dz2,dz3;
	float d12,d13,d23;
	float cth1,cth2,cth3, cth1_,cth2_,cth3_;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3,W;

	for (i=0; i<nodeS[row][col][mom].len-2; ++i){
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=i+1; j<nodeS[row][col][mom].len-1; ++j){
		x2 = nodeS[row][col][mom].elements[j].x;
		y2 = nodeS[row][col][mom].elements[j].y;
		z2 = nodeS[row][col][mom].elements[j].z;
		w2 = nodeS[row][col][mom].elements[j].w;
		dx = x2-x1;
		dy = y2-y1;
		dz1 = z2-z1;
		d12 = dx*dx+dy*dy+dz1*dz1;
		if (d12<dd_max){
		d12 = sqrt(d12);
			for (k=j+1; k<nodeS[row][col][mom].len; ++k){ 
			x3 = nodeS[row][col][mom].elements[k].x;
			y3 = nodeS[row][col][mom].elements[k].y;
			z3 = nodeS[row][col][mom].elements[k].z;
			w3 = nodeS[row][col][mom].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz2 = z3-z1;
			d13 = dx*dx+dy*dy+dz2*dz2;
			if (d13<dd_max){
			d13 = sqrt(d13);
			dx = x3-x2;
			dy = y3-y2;
			dz3 = z3-z2;
			d23 = dx*dx+dy*dy+dz3*dz3;
			if (d23<dd_max){
			d23 = sqrt(d23);
				
				// ángulo entre r y ẑ
				
				cth1 = dz1/d12;
				cth2 = dz2/d13;
				cth3 = dz3/d23;
				
				cth1_ = 1 - cth1;
				cth2_ = 1 - cth2;
				cth3_ = 1 - cth3;
				
				cth1 += 1;
				cth2 += 1;
				cth3 += 1;
				
				// Indices 
				a_ = (int) (d12*ds);
				b_ = (int) (d13*ds);
				c_ = (int) (d23*ds);
				
				t = (int) (cth1*ds_th);
				p = (int) (cth2*ds_th);
				q = (int) (cth3*ds_th);
				
				t_ = (int) (cth1_*ds_th);
				p_ = (int) (cth2_*ds_th);
				q_ = (int) (cth3_*ds_th);
				
				// Guardamos
				W = w1*w2*w3;
				*(*(*(*(*(XXX+a_)+b_)+c_)+t)+p)+=W;
				*(*(*(*(*(XXX+c_)+a_)+b_)+p_)+q_)+=W;
				*(*(*(*(*(XXX+b_)+c_)+a_)+q)+t_)+=W;
				*(*(*(*(*(XXX+b_)+a_)+c_)+p)+t)+=W;
				*(*(*(*(*(XXX+c_)+b_)+a_)+q_)+p_)+=W;
				*(*(*(*(*(XXX+a_)+c_)+b_)+t_)+q)+=W;
			}
			}
			}
		}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N112(int row, int col, int mom, int u, int v, int w, double *****XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en dos 
	nodos con dos puntos en N1 y un punto en N2.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	
	*/
	int i,j,k;
	int a_,b_,c_;
	int t,p,q, t_,p_,q_;
	float dx,dy,dz1,dz2,dz3;
	float d12,d13,d23;
	float cth1,cth2,cth3, cth1_,cth2_,cth3_;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3,W;

	for (i=0; i<nodeS[row][col][mom].len; ++i){
	// 1er punto en N1
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeS[u][v][w].len; ++j){
		// 2do punto en N2
		x2 = nodeS[u][v][w].elements[j].x;
		y2 = nodeS[u][v][w].elements[j].y;
		z2 = nodeS[u][v][w].elements[j].z;
		w2 = nodeS[u][v][w].elements[j].w;
		dx = x2-x1;
		dy = y2-y1;
		dz1 = z2-z1;
		d12 = dx*dx+dy*dy+dz1*dz1;
		if (d12<dd_max){
		d12 = sqrt(d12);
			for (k=i+1; k<nodeS[row][col][mom].len; ++k){
			// 3er punto en N1
			x3 = nodeS[row][col][mom].elements[k].x;
			y3 = nodeS[row][col][mom].elements[k].y;
			z3 = nodeS[row][col][mom].elements[k].z;
			w3 = nodeS[row][col][mom].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz2 = z3-z1;
			d13 = dx*dx+dy*dy+dz2*dz2;
			if (d13<dd_max){
			d13 = sqrt(d13);
			dx = x3-x2;
			dy = y3-y2;
			dz3 = z3-z2;
			d23 = dx*dx+dy*dy+dz3*dz3;
			if (d23<dd_max){
			d23 = sqrt(d23);
				
				// ángulo entre r y ẑ
				
				cth1 = dz1/d12;
				cth2 = dz2/d13;
				cth3 = dz3/d23;
				
				cth1_ = 1 - cth1;
				cth2_ = 1 - cth1;
				cth3_ = 1 - cth1;
				
				cth1 += 1;
				cth2 += 1;
				cth3 += 1;
				
				// Indices 
				a_ = (int) (d12*ds);
				b_ = (int) (d13*ds);
				c_ = (int) (d23*ds);
				
				t = (int) (cth1*ds_th);
				p = (int) (cth2*ds_th);
				q = (int) (cth3*ds_th);
				
				t_ = (int) (cth1_*ds_th);
				p_ = (int) (cth2_*ds_th);
				q_ = (int) (cth3_*ds_th);
				
				// Guardamos
				W = w1*w2*w3;
				*(*(*(*(*(XXX+a_)+b_)+c_)+t)+p)+=W;
				*(*(*(*(*(XXX+c_)+a_)+b_)+p_)+q_)+=W;
				*(*(*(*(*(XXX+b_)+c_)+a_)+q)+t_)+=W;
				*(*(*(*(*(XXX+b_)+a_)+c_)+p)+t)+=W;
				*(*(*(*(*(XXX+c_)+b_)+a_)+q_)+p_)+=W;
				*(*(*(*(*(XXX+a_)+c_)+b_)+t_)+q)+=W;
			}
			}
			}
			for (k=j+1; k<nodeS[u][v][w].len; ++k){
			// 3er punto en N2
			x3 = nodeS[u][v][w].elements[k].x;
			y3 = nodeS[u][v][w].elements[k].y;
			z3 = nodeS[u][v][w].elements[k].z;
			w3 = nodeS[u][v][w].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz2 = z3-z1;
			d13 = dx*dx+dy*dy+dz2*dz2;
			if (d13<dd_max){
			d13 = sqrt(d13);
			dx = x3-x2;
			dy = y3-y2;
			dz3 = z3-z2;
			d23 = dx*dx+dy*dy+dz3*dz3;
			if (d23<dd_max){
			d23 = sqrt(d23);
				
				// ángulo entre r y ẑ
				
				cth1 = dz1/d12;
				cth2 = dz2/d13;
				cth3 = dz3/d23;
				
				cth1_ = 1 - cth1;
				cth2_ = 1 - cth2;
				cth3_ = 1 - cth3;
				
				cth1 += 1;
				cth2 += 1;
				cth3 += 1;
				
				// Indices 
				a_ = (int) (d12*ds);
				b_ = (int) (d13*ds);
				c_ = (int) (d23*ds);
				
				t = (int) (cth1*ds_th);
				p = (int) (cth2*ds_th);
				q = (int) (cth3*ds_th);
				
				t_ = (int) (cth1_*ds_th);
				p_ = (int) (cth2_*ds_th);
				q_ = (int) (cth3_*ds_th);
				
				// Guardamos
				W = w1*w2*w3;
				*(*(*(*(*(XXX+a_)+b_)+c_)+t)+p)+=W;
				*(*(*(*(*(XXX+c_)+a_)+b_)+p_)+q_)+=W;
				*(*(*(*(*(XXX+b_)+c_)+a_)+q)+t_)+=W;
				*(*(*(*(*(XXX+b_)+a_)+c_)+p)+t)+=W;
				*(*(*(*(*(XXX+c_)+b_)+a_)+q_)+p_)+=W;
				*(*(*(*(*(XXX+a_)+c_)+b_)+t_)+q)+=W;
			}
			}
			}
		}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N123(int row, int col, int mom, int u, int v, int w, int a, int b, int c, double *****XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en tres 
	nodos con un puntos en N1, un punto en N2
	y un punto en N3.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	a, b, c => posición de N3.
	
	*/
	int i,j,k;
	int a_,b_,c_;
	int t,p,q, t_,p_,q_;
	float dx,dy,dz1,dz2,dz3;
	float d12,d13,d23;
	float cth1,cth2,cth3, cth1_,cth2_,cth3_;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3,W;
	
	for (i=0; i<nodeS[row][col][mom].len; ++i){
	// 1er punto en N1
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeS[a][b][c].len; ++j){
		// 2do punto en N3
		x2 = nodeS[a][b][c].elements[j].x;
		y2 = nodeS[a][b][c].elements[j].y;
		z2 = nodeS[a][b][c].elements[j].z;
		w2 = nodeS[a][b][c].elements[j].w;
		dx = x2-x1;
		dy = y2-y1;
		dz1 = z2-z1;
		d12 = dx*dx+dy*dy+dz1*dz1;
		if (d12<dd_max){
		d12 = sqrt(d12);
			for (k=0; k<nodeS[u][v][w].len; ++k){
			// 3er punto en N2
			x3 = nodeS[u][v][w].elements[k].x;
			y3 = nodeS[u][v][w].elements[k].y;
			z3 = nodeS[u][v][w].elements[k].z;
			w3 = nodeS[u][v][w].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz2 = z3-z1;
			d13 = dx*dx+dy*dy+dz2*dz2;
			if (d13<dd_max){
			d13 = sqrt(d13);
			dx = x3-x2;
			dy = y3-y2;
			dz3 = z3-z2;
			d23 = dx*dx+dy*dy+dz3*dz3;
			if (d23<dd_max){
			d23 = sqrt(d23);
				
				// ángulo entre r y ẑ
				
				cth1 = dz1/d12;
				cth2 = dz2/d13;
				cth3 = dz3/d23;
				
				cth1_ = 1 - cth1;
				cth2_ = 1 - cth2;
				cth3_ = 1 - cth3;
				
				cth1 += 1;
				cth2 += 1;
				cth3 += 1;
				
				// Indices 
				a_ = (int) (d12*ds);
				b_ = (int) (d13*ds);
				c_ = (int) (d23*ds);
				
				t = (int) (cth1*ds_th);
				p = (int) (cth2*ds_th);
				q = (int) (cth3*ds_th);
				
				t_ = (int) (cth1_*ds_th);
				p_ = (int) (cth2_*ds_th);
				q_ = (int) (cth3_*ds_th);
				
				// Guardamos
				W = w1*w2*w3;
				*(*(*(*(*(XXX+a_)+b_)+c_)+t)+p)+=W;
				*(*(*(*(*(XXX+c_)+a_)+b_)+p_)+q_)+=W;
				*(*(*(*(*(XXX+b_)+c_)+a_)+q)+t_)+=W;
				*(*(*(*(*(XXX+b_)+a_)+c_)+p)+t)+=W;
				*(*(*(*(*(XXX+c_)+b_)+a_)+q_)+p_)+=W;
				*(*(*(*(*(XXX+a_)+c_)+b_)+t_)+q)+=W;
				
			}
			}
			}
		}
		}
	}
}
//=================================================================== 
void NODE3P::make_histoXXY(double *****XXY, Node ***nodeX, Node ***nodeY){
	/*
	Función para crear los histogramas DDR y DRR.
	
	Argumentos
	XXY: arreglo donde se creará el histograma.
	nodeX: malla 1.
	nodeY: malla 2.
	
	*/ 
	int i, j, k, row, col, mom, u, v, w, a, b, c;
	float dis, dis_nod, dis_nod2, dis_nod3;
	float x1N, y1N, z1N, x2N, y2N, z2N, x3N, y3N, z3N;
	float x, y, z;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod, dx_nod2, dy_nod2, dz_nod2, dx_nod3, dy_nod3, dz_nod3;
	bool con_x, con_y, con_z;
	
	// x1N, y1N, z1N => Nodo pivote
	for (row=0; row<partitions; ++row){
	x1N = nodeX[row][0][0].nodepos.x;
	for (col=0; col<partitions; ++col){
	y1N = nodeX[row][col][0].nodepos.y;
	for (mom=0; mom<partitions; ++mom){
	z1N = nodeX[row][col][mom].nodepos.z;			
				
	//=======================
	// Nodo 2 movil en ZYX:
	//=======================
	for (u=0; u<partitions; ++u){
	x2N = nodeY[u][0][0].nodepos.x;
	dx_nod = x2N-x1N;
	dx_nod *= dx_nod;
		for (v=0; v<partitions; ++v){
		y2N = nodeY[u][v][0].nodepos.y;
		dy_nod = y2N-y1N;
		dy_nod *= dy_nod;
			for (w=0; w<partitions; ++w){
			z2N = nodeY[u][v][w].nodepos.z;
			dz_nod = z2N-z1N;
			dz_nod *= dz_nod;
			dis_nod = dx_nod + dy_nod + dz_nod;
			if (dis_nod <= ddmax_nod){
			//==============================================
			// 2 puntos en N y 1 punto en N'
			//==============================================
			count_3_N112_xxy(row, col, mom, u, v, w, XXY, nodeX, nodeY);
			//==============================================
			// 1 punto en N1, 1 punto en N2 y 1 punto en N3
			//==============================================
			a = row;
			b = col;
			//=======================
			// Nodo 3 movil en Z:
			//=======================
			x3N = nodeX[a][0][0].nodepos.x;
			y3N = nodeX[a][b][0].nodepos.y;
			dx_nod2 = x3N-x1N;
			dx_nod2 *= dx_nod2;
			dy_nod2 = y3N-y1N;
			dy_nod2 *= dy_nod2;
				for (c=mom+1;  c<partitions; ++c){	
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod2 = z3N-z1N;
				dz_nod2 *= dz_nod2;
				dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
				if (dis_nod2 <= ddmax_nod){
				dz_nod3 = z3N-z2N;
				dis_nod3 = dz_nod3*dz_nod3;
				if (dis_nod3 <= ddmax_nod){
					count_3_N123_xxy(row, col, mom, u, v, w, a, b, c, XXY, nodeX, nodeY);
				}
				}
				}
			//=======================
			// Nodo 3 movil en ZY:
			//=======================
				for (b=col+1; b<partitions; ++b){
				y3N = nodeX[a][b][0].nodepos.y;
				dy_nod2 = y3N-y1N;
				dy_nod2 *= dy_nod2;
					for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					dz_nod2 *= dz_nod2;
					dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
					if (dis_nod2 <= ddmax_nod){
					dy_nod3 = y3N-y2N;
					dy_nod3 *= dy_nod3;
					dz_nod3 = z3N-z2N;
					dz_nod3 *= dz_nod3;
					dis_nod3 = dy_nod3 + dz_nod3;
					if (dis_nod3 <= ddmax_nod){
						count_3_N123_xxy(row, col, mom, u, v, w, a, b, c, XXY, nodeX, nodeY);
					}
					}
				}
				}
			//=======================
			// Nodo 3 movil en ZYX:
			//=======================		
				for (a=row+1; a<partitions; ++a){
				x3N = nodeX[a][0][0].nodepos.x;
				dx_nod2 = x3N-x1N;
				dx_nod2 *= dx_nod2;
					for (b=0; b<partitions; ++b){
					y3N = nodeX[a][b][0].nodepos.y;
					dy_nod2 = y3N-y1N;
					dy_nod2 *= dy_nod2;
						for (c=0;  c<partitions; ++c){
						z3N = nodeX[a][b][c].nodepos.z;
						dz_nod2 = z3N-z1N;
						dz_nod2 *= dz_nod2;
						dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
						if (dis_nod2 <= ddmax_nod){
						dx_nod3 = x3N-x2N;
						dx_nod3 *= dx_nod3;
						dy_nod3 = y3N-y2N;
						dy_nod3 *= dy_nod3;
						dz_nod3 = z3N-z2N;
						dz_nod3 *= dz_nod3;
						dis_nod3 = dx_nod3 + dy_nod3 + dz_nod3;
						if (dis_nod3 <= ddmax_nod){
							count_3_N123_xxy(row, col, mom, u, v, w, a, b, c, XXY, nodeX, nodeY);
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
}
//=================================================================== 
void NODE3P::count_3_N112_xxy(int row, int col, int mom, int u, int v, int w, double *****XXY, Node ***nodeS, Node ***nodeT){
	/*
	Funcion para contar los triangulos en dos 
	nodos con dos puntos en N1 y un punto en N2.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	*/
	int i,j,k;
	int a_,b_,c_;
	int t,p, p_,q_;
	float dx,dy,dz1,dz2,dz3;
	float d12,d13,d23;
	float cth1,cth2, cth2_,cth3_;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3,W;

	for (i=0; i<nodeS[row][col][mom].len; ++i){
	// 1er punto en N1
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeT[u][v][w].len; ++j){
		// 2do punto en N2
		x2 = nodeT[u][v][w].elements[j].x;
		y2 = nodeT[u][v][w].elements[j].y;
		z2 = nodeT[u][v][w].elements[j].z;
		w2 = nodeT[u][v][w].elements[j].w;
		dx = x2-x1;
		dy = y2-y1;
		dz1 = z2-z1;
		d12 = dx*dx+dy*dy+dz1*dz1;
		if (d12<dd_max){
		d12 = sqrt(d12);
			for (k=i+1; k<nodeS[row][col][mom].len; ++k){
			// 3er punto en N1
			x3 = nodeS[row][col][mom].elements[k].x;
			y3 = nodeS[row][col][mom].elements[k].y;
			z3 = nodeS[row][col][mom].elements[k].z;
			w3 = nodeS[row][col][mom].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz2 = z3-z1;
			d13 = dx*dx+dy*dy+dz2*dz2;
			if (d13<dd_max){
			d13 = sqrt(d13);
			dx = x3-x2;
			dy = y3-y2;
			dz3 = z3-z2;
			d23 = dx*dx+dy*dy+dz3*dz3;
			if (d23<dd_max){
			d23 = sqrt(d23);
				
				// ángulo entre r y ẑ
				
				cth1 = dz1/d12;
				cth2 = dz2/d13;
				
				cth2_ = 1 - cth2;
				cth3_ = 1 - (dz3/d23);
				
				cth1 += 1;
				cth2 += 1;
				
				// Indices 
				a_ = (int) (d12*ds);
				b_ = (int) (d13*ds);
				c_ = (int) (d23*ds);
				
				t = (int) (cth1*ds_th);
				p = (int) (cth2*ds_th);
				
				p_ = (int) (cth2_*ds_th);
				q_ = (int) (cth3_*ds_th);
					
				// Guardamos
				W = w1*w2*w3;
				*(*(*(*(*(XXY+a_)+b_)+c_)+t)+p)+=W;
				*(*(*(*(*(XXY+b_)+a_)+c_)+p)+t)+=W;
				*(*(*(*(*(XXY+c_)+b_)+a_)+q_)+p_)+=W;
				*(*(*(*(*(XXY+c_)+a_)+b_)+p_)+q_)+=W;
			
			}
			}
			}
		}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N123_xxy(int row, int col, int mom, int u, int v, int w, int a, int b, int c, double *****XXY, Node ***nodeS, Node ***nodeT){
	/*
	Funcion para contar los triangulos en tres 
	nodos con un puntos en N1, un punto en N2
	y un punto en N3.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	a, b, c => posición de N3.
	*/
	int i,j,k;
	int a_,b_,c_;
	int t,p, p_,q_;
	float dx,dy,dz1,dz2,dz3;
	float d12,d13,d23;
	float cth1,cth2, cth2_,cth3_;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3,W;

	for (i=0; i<nodeS[row][col][mom].len; ++i){
		// 1er punto en N1
		x1 = nodeS[row][col][mom].elements[i].x;
		y1 = nodeS[row][col][mom].elements[i].y;
		z1 = nodeS[row][col][mom].elements[i].z;
		w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeT[u][v][w].len; ++j){
			// 2do punto en N3
			x2 = nodeT[u][v][w].elements[j].x;
			y2 = nodeT[u][v][w].elements[j].y;
			z2 = nodeT[u][v][w].elements[j].z;
			w2 = nodeT[u][v][w].elements[j].w;
			dx = x2-x1;
			dy = y2-y1;
			dz1 = z2-z1;
			d12 = dx*dx+dy*dy+dz1*dz1;
			if (d12<dd_max){
			d12 = sqrt(d12);
			for (k=0; k<nodeS[a][b][c].len; ++k){
				// 3er punto en N2
				x3 = nodeS[a][b][c].elements[k].x;
				y3 = nodeS[a][b][c].elements[k].y;
				z3 = nodeS[a][b][c].elements[k].z;
				w3 = nodeS[a][b][c].elements[k].w;
				dx = x3-x1;
				dy = y3-y1;
				dz2 = z3-z1;
				d13 = dx*dx+dy*dy+dz2*dz2;
				if (d13<dd_max){
				d13 = sqrt(d13);
				dx = x3-x2;
				dy = y3-y2;
				dz3 = z3-z2;
				d23 = dx*dx+dy*dy+dz3*dz3;
				if (d23<dd_max){
				d23 = sqrt(d23);
				
					// ángulo entre r y ẑ
					
					cth1 = dz1/d12;
					cth2 = dz2/d13;
					
					cth2_ = 1 - cth2;
					cth3_ = 1 - (dz3/d23);
					
					cth1 += 1;
					cth2 += 1;
					
					// Indices 
					a_ = (int) (d12*ds);
					b_ = (int) (d13*ds);
					c_ = (int) (d23*ds);
					
					t = (int) (cth1*ds_th);
					p = (int) (cth2*ds_th);
					
					p_ = (int) (cth2_*ds_th);
					q_ = (int) (cth3_*ds_th);
						
					// Guardamos
					W = w1*w2*w3;
					*(*(*(*(*(XXY+a_)+b_)+c_)+t)+p)+=W;
					*(*(*(*(*(XXY+b_)+a_)+c_)+p)+t)+=W;
					*(*(*(*(*(XXY+c_)+b_)+a_)+q_)+p_)+=W;
					*(*(*(*(*(XXY+c_)+a_)+b_)+p_)+q_)+=W;
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


