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
		void count_3_N112_front(int, int, int, int, int, int, unsigned int ***, Node ***, short int, short int, short int);
		void count_3_N123_front(int, int, int, int, int, int, int, int, int, unsigned int ***, Node ***, short int, short int, short int, short int, short int, short int);
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
	
	short int sx1, sy1, sz1, sx2, sy2, sz2;
	float d_max_nod = size_box - (d_max + size_node);
	float _d_max_nod = - d_max_nod;
	float dd_max_nod = d_max_nod*d_max_nod;
	
	//=======================
	// Nodo 2 movil en ZYX:
	//=======================
	
	// x1N, y1N, z1N => Nodo pivote
	for (row=0; row<partitions; ++row){
	x1N = nodeX[row][0][0].nodepos.x;
	for (col=0; col<partitions; ++col){
	y1N = nodeX[row][col][0].nodepos.y;
	for (mom=0; mom<partitions; ++mom){
	z1N = nodeX[row][col][mom].nodepos.z;	
	
	u = row;
	v = col;
	
	x2N = nodeX[u][0][0].nodepos.x;
	sx1 = 0;

	y2N = nodeX[u][v][0].nodepos.y;
	sy1 = 0;

	//=======================
	// Nodo 2 movil en Z:
	//=======================
	for (w=mom+1;  w<partitions; ++w){
		z2N = nodeX[u][v][w].nodepos.z;
		dz_nod = z2N-z1N;
		sz1 = 0;
		if ( d_max_nod <= dz_nod ){
			dz_nod = size_box - dz_nod;
			sz1 = 1;
		}
		else if (_d_max_nod >= dz_nod ){
			dz_nod = size_box + dz_nod;
			sz1 = 1;
		}
		if (sx1==1 || sy1==1 || sz1==1){
		dis_nod = dz_nod*dz_nod;
		if (dis_nod <= ddmax_nod){
			//==============================================
			// 2 puntos en N y 1 punto en N'
			//==============================================
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, sx1, sy1, sz1);
			//==============================================
			// 1 punto en N1, 1 punto en N2 y un punto en N3
			//==============================================
			a = u;
			b = v;
			//=======================
			// Nodo 3 movil en Z:
			//=======================
			sx2 = 0;
			sy2 = 0;
			for (c=w+1;  c<partitions; ++c){
				z3N = nodeX[a][b][c].nodepos.z; 
				dz_nod2 = z3N-z1N;
				sz2 = 0;
				if (d_max_nod <= dz_nod2){
					dz_nod2 = size_box - dz_nod2;
					sz2 = 1;
				}
				else if (_d_max_nod >= dz_nod2){
					dz_nod2 = size_box + dz_nod2;
					sz2 = 1;
				}
				if (sz2==1){
				dis_nod2 = dz_nod2*dz_nod2;
				if (dis_nod2 <= ddmax_nod){
				
				dz_nod3 = z3N-z2N;
				dz_nod3 *= dz_nod3;
				if ( dd_max_nod <= dz_nod3 ){
					dz_nod3 = size_box - sqrt(dz_nod3);
					dz_nod3 *= dz_nod3;
				}
				dis_nod3 = dz_nod3;
				if (dis_nod3 <= ddmax_nod){
				
				count_3_N123_front(row, col, mom, u, v, w, a, b, c, XXX, nodeX, sx1, sy1, sz1, sx2, sy2, sz2);
				
				}
				}
				}
			}
			//=======================
			// Nodo 3 movil en ZY:
			//=======================
			sx2 = 0;
			for (b=v+1; b<partitions; ++b){
				y3N = nodeX[a][b][0].nodepos.y;
				dy_nod2 = y3N-y1N;
				sy2 = 0;
				if (d_max_nod <= dy_nod2){
					dy_nod2 = size_box - dy_nod2;
					sy2 = 1;
				}
				else if (_d_max_nod >= dy_nod2){
					dy_nod2 = size_box + dy_nod2;
					sy2 = 1;
				}
				for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					sz2 = 0;
					if (d_max_nod <= dz_nod2){
						dz_nod2 = size_box - dz_nod2;
						sz2 = 1;
					}
					else if (_d_max_nod >= dz_nod2){
						dz_nod2 = size_box + dz_nod2;
						sz2 = 1;
					}
					if (sy2==1 || sz2==1){
					dis_nod2 = dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
					if (dis_nod2 <= ddmax_nod){
					//==================================================================
					
					dy_nod3 = y3N-y2N;
					dy_nod3 *= dy_nod3;
					if ( dd_max_nod <= dy_nod3 ){
						dy_nod3 = size_box - sqrt(dy_nod3);
						dy_nod3 *= dy_nod3;
					}
					
					dz_nod3 = z3N-z2N;
					dz_nod3 *= dz_nod3;
					if ( dd_max_nod <= dz_nod3 ){
						dz_nod3 = size_box - sqrt(dz_nod3);
						dz_nod3 *= dz_nod3;
					}
					
					dis_nod3 = dy_nod3 + dz_nod3;
					if (dis_nod3 <= ddmax_nod){
					count_3_N123_front(row, col, mom, u, v, w, a, b, c, XXX, nodeX, sx1, sy1, sz1, sx2, sy2, sz2);
					}
					//==================================================================
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
				sx2=0;
				if (d_max_nod <= dx_nod2){
					dx_nod2 = size_box - dx_nod2;
					sx2=1;
				}
				else if (_d_max_nod >= dx_nod2){
					dx_nod2 = size_box + dx_nod2;
					sx2=1;
				}
				for (b=0; b<partitions; ++b){
					y3N = nodeX[a][b][0].nodepos.y;
					dy_nod2 = y3N-y1N;
					sy2=0;
					if (d_max_nod <= dy_nod2){
						dy_nod2 = size_box - dy_nod2;
						sy2=1;
					}
					else if (_d_max_nod >= dy_nod2){
						dy_nod2 = size_box + dy_nod2;
						sy2=1;
					}
					for (c=0;  c<partitions; ++c){
						z3N = nodeX[a][b][c].nodepos.z;
						dz_nod2 = z3N-z1N;
						sz2=0;
						if (d_max_nod <= dz_nod2){
							dz_nod2 = size_box - dz_nod2;
							sz2=1;
						}
						else if (_d_max_nod >= dz_nod2){
							dz_nod2 = size_box + dz_nod2;
							sz2=1;
						}
						if (sx2==1 || sy2==1 || sz2==1){
						dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
						if (dis_nod2 <= ddmax_nod){
						//==================================================================
						dx_nod3 = x3N-x2N;
						dx_nod3 *= dx_nod3;
						if ( dd_max_nod <= dx_nod3 ){
							dx_nod3 = size_box - sqrt(dx_nod3);
							dx_nod3 *= dx_nod3;
						}
						
						dy_nod3 = y3N-y2N;
						dy_nod3 *= dy_nod3;
						if ( dd_max_nod <= dy_nod3 ){
							dy_nod3 = size_box - sqrt(dy_nod3);
							dy_nod3 *= dy_nod3;
						}
						
						dz_nod3 = z3N-z2N;
						dz_nod3 *= dz_nod3;
						if ( dd_max_nod <= dz_nod3 ){
							dz_nod3 = size_box - sqrt(dz_nod3);
							dz_nod3 *= dz_nod3;
						}
						
						dis_nod3 = dx_nod3 + dy_nod3 + dz_nod3;
						if (dis_nod3 <= ddmax_nod){
						count_3_N123_front(row, col, mom, u, v, w, a, b, c, XXX, nodeX, sx1, sy1, sz1, sx2, sy2, sz2);
						}
						//==================================================================
						}
						}
					}
				}
			}		
		}
		}
	}
	
	u = row;
	x2N = nodeX[u][0][0].nodepos.x;
	sx1=0;
	//=======================
	// Nodo 2 movil en ZY:
	//=======================
	for (v=col+1; v<partitions ; ++v){
		y2N = nodeX[u][v][0].nodepos.y;
		dy_nod = y2N-y1N;
		sy1=0;
		if (d_max_nod <= dy_nod){
			dy_nod = size_box - dy_nod;
			sy1=1;
		}
		else if (_d_max_nod >= dy_nod){
			dy_nod = size_box + dy_nod;
			sy1=1;
		}
		for (w=0; w<partitions ; ++w){		
			z2N = nodeX[u][v][w].nodepos.z;
			dz_nod = z2N-z1N;
			sz1=0;
			if (d_max_nod <= dz_nod){
				dz_nod = size_box - dz_nod;
				sz1=1;
			}
			else if (_d_max_nod >= dz_nod){
				dz_nod = size_box + dz_nod;
				sz1=1;
			}
			if (sx1==1 || sy1==1 || sz1==1){
			dis_nod = dy_nod*dy_nod + dz_nod*dz_nod;
			if (dis_nod <= ddmax_nod){
			//==============================================
			// 2 puntos en N y 1 punto en N'
			//==============================================
			count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, sx1, sy1, sz1);
			//==============================================
			// 1 punto en N1, 1 punto en N2 y un punto en N3
			//==============================================
			a = u;
			b = v;
			//=======================
			// Nodo 3 movil en Z:
			//=======================
			
			sx2 = 0;
			
			y3N = nodeX[a][b][0].nodepos.y;
			dy_nod2 = y3N-y1N;
			sy2 = 0;
			if (d_max_nod <= dy_nod2){
				dy_nod2 = size_box - dy_nod2;
				sy2 = 1;
			}
			else if (_d_max_nod >= dy_nod2){
				dy_nod2 = size_box + dy_nod2;
				sy2 = 1;
			}
			for (c=w+1;  c<partitions; ++c){
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod2 = z3N-z1N;
				sz2 = 0;
				if (d_max_nod <= dz_nod2){
					dz_nod2 = size_box - dz_nod2;
					sz2 = 1;
				}
				else if (_d_max_nod >= dz_nod2){
					dz_nod2 = size_box + dz_nod2;
					sz2 = 1;
				}
				if (sy2==1 || sz2==1){
				dis_nod2 = dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
				if (dis_nod2 <= ddmax_nod){
				//==================================================================
				
				dz_nod3 = z3N-z2N;
				dz_nod3 *= dz_nod3;
				if ( dd_max_nod <= dz_nod3 ){
					dz_nod3 = size_box - sqrt(dz_nod3);
					dz_nod3 *= dz_nod3;
				}
				
				dis_nod3 = dz_nod3;
				if (dis_nod3 <= ddmax_nod){
				count_3_N123_front(row, col, mom, u, v, w, a, b, c, XXX, nodeX, sx1, sy1, sz1, sx2, sy2, sz2);
				}
				//==================================================================
				}
				}
			}
			//=======================
			// Nodo 3 movil en ZY:
			//=======================
			sx2 = 0;
			for (b=v+1; b<partitions; ++b){
				y3N = nodeX[a][b][0].nodepos.y;
				dy_nod2 = y3N-y1N;
				sy2 = 0;
				if (d_max_nod <= dy_nod2){
					dy_nod2 = size_box - dy_nod2;
					sy2 = 1;
				}
				else if (_d_max_nod >= dy_nod2){
					dy_nod2 = size_box + dy_nod2;
					sy2 = 1;
				}
				for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					sz2 = 0;
					if (d_max_nod <= dz_nod2){
						dz_nod2 = size_box - dz_nod2;
						sz2 = 1;
					}
					else if (_d_max_nod >= dz_nod2){
						dz_nod2 = size_box + dz_nod2;
						sz2 = 1;
					}
					if (sy2==1 || sz2==1){
					dis_nod2 = dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
					if (dis_nod2 <= ddmax_nod){
					//==================================================================
					dy_nod3 = y3N-y2N;
					dy_nod3 *= dy_nod3;
					if ( dd_max_nod <= dy_nod3 ){
						dy_nod3 = size_box - sqrt(dy_nod3);
						dy_nod3 *= dy_nod3;
					}
					
					dz_nod3 = z3N-z2N;
					dz_nod3 *= dz_nod3;
					if ( dd_max_nod <= dz_nod3 ){
						dz_nod3 = size_box - sqrt(dz_nod3);
						dz_nod3 *= dz_nod3;
					}
					
					dis_nod3 = dy_nod3 + dz_nod3;
					if (dis_nod3 <= ddmax_nod){
					count_3_N123_front(row, col, mom, u, v, w, a, b, c, XXX, nodeX, sx1, sy1, sz1, sx2, sy2, sz2);
					}
					//==================================================================
					}
					}
				}
			}
			//=======================
			// Nodo 3 movil en ZYX:
			//=======================
			for (a=u+1; a<partitions; ++a){
				x3N = nodeX[a][0][0].nodepos.x;
				dx_nod2 = x2N-x1N;
				sx2 = 0;
				if (d_max_nod <= dx_nod2){
					dx_nod2 = size_box - dx_nod2;
					sx2 = 1;
				}
				else if (_d_max_nod >= dx_nod2){
					dx_nod2 = size_box + dx_nod2;
					sx2 = 1;
				}
				for (b=0; b<partitions; ++b){
					y3N = nodeX[a][b][0].nodepos.y;
					dy_nod2 = y3N-y1N;
					sy2 = 0;
					if (d_max_nod <= dy_nod2){
						dy_nod2 = size_box - dy_nod2;
						sy2 = 1;
					}
					else if (_d_max_nod >= dy_nod2){
						dy_nod2 = size_box + dy_nod2;
						sy2 = 1;
					}
					for (c=0;  c<partitions; ++c){
						z3N = nodeX[a][b][c].nodepos.z;
						dz_nod2 = z3N-z1N;
						sz2 = 0;
						if (d_max_nod <= dz_nod2){
							dz_nod2 = size_box - dz_nod2;
							sz2 = 1;
						}
						else if (_d_max_nod >= dz_nod2){
							dz_nod2 = size_box + dz_nod2;
							sz2 = 1;
						}
						if (sx2==1 || sy2==1 || sz2==1){
						dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
						if (dis_nod2 <= ddmax_nod){
						//==================================================================
						dx_nod3 = x3N-x2N;
						dx_nod3 *= dx_nod3;
						if ( dd_max_nod <= dx_nod3 ){
							dx_nod3 = size_box - sqrt(dx_nod3);
							dx_nod3 *= dx_nod3;
						}
						
						dy_nod3 = y3N-y2N;
						dy_nod3 *= dy_nod3;
						if ( dd_max_nod <= dy_nod3 ){
							dy_nod3 = size_box - sqrt(dy_nod3);
							dy_nod3 *= dy_nod3;
						}
						
						dz_nod3 = z3N-z2N;
						dz_nod3 *= dz_nod3;
						if ( dd_max_nod <= dz_nod3 ){
							dz_nod3 = size_box - sqrt(dz_nod3);
							dz_nod3 *= dz_nod3;
						}
						
						dis_nod3 = dx_nod3 + dy_nod3 + dz_nod3;
						if (dis_nod3 <= ddmax_nod){
						count_3_N123_front(row, col, mom, u, v, w, a, b, c, XXX, nodeX, sx1, sy1, sz1, sx2, sy2, sz2);
						}
						//==================================================================
						}
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
		sx1 = 0;
		if (d_max_nod <= dx_nod){
			dx_nod = size_box - dx_nod;
			sx1 = 1;
		}
		else if (_d_max_nod >= dx_nod){
			dx_nod = size_box + dx_nod;
			sx1 = 1;
		}
		for (v=0; v<partitions; ++v){
			y2N = nodeX[u][v][0].nodepos.y;
			dy_nod = y2N-y1N;
			sy1 = 0;
			if (d_max_nod <= dy_nod){
				dy_nod = size_box - dy_nod;
				sy1 = 1;
			}
			else if (_d_max_nod >= dy_nod){
				dy_nod = size_box + dy_nod;
				sy1 = 1;
			}
			for (w=0; w<partitions; ++w){
				z2N = nodeX[u][v][w].nodepos.z;
				dz_nod = z2N-z1N;
				sz1 = 0;
				if (d_max_nod <= dz_nod){
					dz_nod = size_box - dz_nod;
					sz1 = 1;
				}
				else if (_d_max_nod >= dz_nod){
					dz_nod = size_box + dz_nod;
					sz1 = 1;
				}
				if ( sx1==1 || sy1==1 || sz1==1){ //tenemos frontera
				dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod; // distancia entre nodo caja y nodo frontera
				if (dis_nod <= ddmax_nod){
				//==============================================
				// 2 puntos en N y 1 punto en N'
				//==============================================
				count_3_N112_front(row, col, mom, u, v, w, XXX, nodeX, sx1, sy1, sz1);
				//==============================================
				// 1 punto en N1, 1 punto en N2 y un punto en N3
				//==============================================
				a = u;
				b = v;
				//=======================
				// Nodo 3 movil en Z:
				//=======================
				x3N = nodeX[a][0][0].nodepos.x;
				dx_nod2 = x3N-x1N;
				sx2 = 0;
				if (d_max_nod <= dx_nod2){
					dx_nod2 = size_box - dx_nod2;
					sx2 = 1;
				}
				else if (_d_max_nod >= dx_nod2){
					dx_nod2 = size_box + dx_nod2;
					sx2 = 1;
				}
				y3N = nodeX[a][b][0].nodepos.y;
				dy_nod2 = y3N-y1N;
				sy2 = 0;
				if (d_max_nod <= dy_nod2){
					dy_nod2 = size_box - dy_nod2;
					sy2 = 1;
				}
				else if (_d_max_nod >= dy_nod2){
					dy_nod2 = size_box + dy_nod2;
					sy2 = 1;
				}
				for (c=w+1;  c<partitions; ++c){	
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					sz2 = 0;
					if (d_max_nod <= dz_nod2){
						dz_nod2 = size_box - dz_nod2;
						sz2 = 1;
					}
					else if (_d_max_nod >= dz_nod2){
						dz_nod2 = size_box + dz_nod2;
						sz2 = 1;
					}
					if (sx2==1 || sy2==1 || sz2==1){
					dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
					if (dis_nod2 <= ddmax_nod){
					//==================================================================
					
					dz_nod3 = z3N-z2N;
					dz_nod3 *= dz_nod3;
					if ( dd_max_nod <= dz_nod3 ){
						dz_nod3 = size_box - sqrt(dz_nod3);
						dz_nod3 *= dz_nod3;
					}
					
					dis_nod3 = dz_nod3;
					if (dis_nod3 <= ddmax_nod){
					count_3_N123_front(row, col, mom, u, v, w, a, b, c, XXX, nodeX, sx1, sy1, sz1, sx2, sy2, sz2);
					}
					//==================================================================
					}
					}
				}
				//=======================
				// Nodo 3 movil en ZY:
				//=======================
				x3N = nodeX[a][0][0].nodepos.x;
				dx_nod2 = x3N-x1N;
				sx2 = 0;
				if (d_max_nod <= dx_nod2){
					dx_nod2 = size_box - dx_nod2;
					sx2 = 1;
				}
				else if (_d_max_nod >= dx_nod2){
					dx_nod2 = size_box + dx_nod2;
					sx2 = 1;
				}
				for (b=v+1; b<partitions; ++b){
					y3N = nodeX[a][b][0].nodepos.y;
					dy_nod2 = y3N-y1N;
					sy2 = 0;
					if (d_max_nod <= dy_nod2){
						dy_nod2 = size_box - dy_nod2;
						sy2 = 1;
					}
					else if (_d_max_nod >= dy_nod2){
						dy_nod2 = size_box + dy_nod2;
						sy2 = 1;
					}
					for (c=0;  c<partitions; ++c){
						z3N = nodeX[a][b][c].nodepos.z;
						dz_nod2 = z3N-z1N;
						sz2 = 0;
						if (d_max_nod <= dz_nod2){
							dz_nod2 = size_box - dz_nod2;
							sz2 = 1;
						}
						else if (_d_max_nod >= dz_nod2){
							dz_nod2 = size_box + dz_nod2;
							sz2 = 1;
						}
						if (sx2==1 || sy2==1 || sz2==1){
						dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
						if (dis_nod2 <= ddmax_nod){
						//==================================================================
						
						dy_nod3 = y3N-y2N;
						dy_nod3 *= dy_nod3;
						if ( dd_max_nod <= dy_nod3 ){
							dy_nod3 = size_box - sqrt(dy_nod3);
							dy_nod3 *= dy_nod3;
						}
					
						dz_nod3 = z3N-z2N;
						dz_nod3 *= dz_nod3;
						if ( dd_max_nod <= dz_nod3 ){
							dz_nod3 = size_box - sqrt(dz_nod3);
							dz_nod3 *= dz_nod3;
						}
						
						dis_nod3 = dy_nod3 + dz_nod3;
						if (dis_nod3 <= ddmax_nod){
						count_3_N123_front(row, col, mom, u, v, w, a, b, c, XXX, nodeX, sx1, sy1, sz1, sx2, sy2, sz2);
						}
						//==================================================================
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
					sx2 = 0;
					if (d_max_nod <= dx_nod2){
						dx_nod2 = size_box - dx_nod2;
						sx2 = 1;
					}
					else if (_d_max_nod >= dx_nod2){
						dx_nod2 = size_box + dx_nod2;
						sx2 = 1;
					}
					for (b=0; b<partitions; ++b){
						y3N = nodeX[a][b][0].nodepos.y;
						dy_nod2 = y3N-y1N;
						sy2 = 0;
						if (d_max_nod <= dy_nod2){
							dy_nod2 = size_box - dy_nod2;
							sy2 = 1;
						}
						else if (_d_max_nod >= dy_nod2){
							dy_nod2 = size_box + dy_nod2;
							sy2 = 1;
						}
						for (c=0;  c<partitions; ++c){
							z3N = nodeX[a][b][c].nodepos.z;
							dz_nod2 = z3N-z1N;
							sz2 = 0;
							if (d_max_nod <= dz_nod2){
								dz_nod2 = size_box - dz_nod2;
								sz2 = 1;
							}
							else if (_d_max_nod >= dz_nod2){
								dz_nod2 = size_box + dz_nod2;
								sz2 = 1;
							}
							
							if (sx2==1 || sy2==1 || sz2==1){
							
							dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
							if (dis_nod2 <= ddmax_nod){
							//==================================================================
							
							dx_nod3 = x3N-x2N;
							dx_nod3 *= dx_nod3;
							if ( dd_max_nod <= dx_nod3 ){
								dx_nod3 = size_box - sqrt(dx_nod3);
								dx_nod3 *= dx_nod3;
							}
							
							dy_nod3 = y3N-y2N;
							dy_nod3 *= dy_nod3;
							if ( dd_max_nod <= dy_nod3 ){
								dy_nod3 = size_box - sqrt(dy_nod3);
								dy_nod3 *= dy_nod3;
							}
							
							dz_nod3 = z3N-z2N;
							dz_nod3 *= dz_nod3;
							if ( dd_max_nod <= dz_nod3 ){
								dz_nod3 = size_box - sqrt(dz_nod3);
								dz_nod3 *= dz_nod3;
							}
							
							dis_nod3 = dx_nod3 + dy_nod3 + dz_nod3;
							if (dis_nod3 <= ddmax_nod){
							count_3_N123_front(row, col, mom, u, v, w, a, b, c, XXX, nodeX, sx1, sy1, sz1, sx2, sy2, sz2);
							}
							//==================================================================
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
void NODE3P::count_3_N112_front(int row, int col, int mom, int u, int v, int w, unsigned int ***XXX, Node ***nodeS, short int sx1, short int sy1, short int sz1){
	/*
	Funcion para contar los triangulos en dos 
	nodos con dos puntos en N1 y un punto en N2.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	
	*/
	int i,j,k;
	float dx12,dy12,dz12,dx13,dy13,dz13,dx23,dy23,dz23;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	int a, b, c;
	for (i=0; i<nodeS[row][col][mom].len; ++i){
		// 1er punto en N1 (caja)
		x1 = nodeS[row][col][mom].elements[i].x;
		y1 = nodeS[row][col][mom].elements[i].y;
		z1 = nodeS[row][col][mom].elements[i].z;
		for (j=0; j<nodeS[u][v][w].len; ++j){
			// 2do punto en N2 (frontera)
			x2 = nodeS[u][v][w].elements[j].x;
			y2 = nodeS[u][v][w].elements[j].y;
			z2 = nodeS[u][v][w].elements[j].z;
			dx12 = (size_box*(sx1)) - fabs(x2-x1);
			dy12 = (size_box*(sy1)) - fabs(y2-y1);
			dz12 = (size_box*(sz1)) - fabs(z2-z1);
			d12 = dx12*dx12+dy12*dy12+dz12*dz12;
			if (d12<=dd_max){
			for (k=j+1; k<nodeS[u][v][w].len; ++k){
				// 3er punto en N2 (frontera)
				x3 = nodeS[u][v][w].elements[k].x;
				y3 = nodeS[u][v][w].elements[k].y;
				z3 = nodeS[u][v][w].elements[k].z;
				dx13 = (size_box*(sx1)) - fabs(x3-x1);
				dy13 = (size_box*(sy1)) - fabs(y3-y1);
				dz13 = (size_box*(sz1)) - fabs(z3-z1);
				d13 = dx13*dx13+dy13*dy13+dz13*dz13;
				if (d13<=dd_max){
				dx23 = x3-x2;
				dy23 = y3-y2;
				dz23 = z3-z2;
				d23 = dx23*dx23 + dy23*dy23 + dz23*dz23;
				if (d23<=dd_max){
				a = (int)(sqrt(d12)*ds);
				b = (int)(sqrt(d13)*ds);
				c = (int)(sqrt(d23)*ds);
				*(*(*(XXX+a)+b)+c)+=1;
				*(*(*(XXX+c)+a)+b)+=1;
				*(*(*(XXX+b)+c)+a)+=1;
				*(*(*(XXX+b)+a)+c)+=1;
				*(*(*(XXX+c)+b)+a)+=1;
				*(*(*(XXX+a)+c)+b)+=1;
				}
				}
			}
			for (k=i+1; k<nodeS[row][col][mom].len; ++k){
				// 3er punto en N1 
				x3 = nodeS[row][col][mom].elements[k].x;
				y3 = nodeS[row][col][mom].elements[k].y;
				z3 = nodeS[row][col][mom].elements[k].z;
				dx13 = x3-x1;
				dy13 = y3-y1;
				dz13 = z3-z1;
				d13 = dx13*dx13+dy13*dy13+dz13*dz13;
				if (d13<=dd_max){
				dx23 = (size_box*(sx1)) - fabs(x3-x2);
				dy23 = (size_box*(sy1)) - fabs(y3-y2);
				dz23 = (size_box*(sz1)) - fabs(z3-z2);
				d23 = dx23*dx23+dy23*dy23+dz23*dz23;
				if (d23<=dd_max){
				a = (int)(sqrt(d12)*ds);
				b = (int)(sqrt(d13)*ds);
				c = (int)(sqrt(d23)*ds);
				*(*(*(XXX+a)+b)+c)+=1;
				*(*(*(XXX+c)+a)+b)+=1;
				*(*(*(XXX+b)+c)+a)+=1;
				*(*(*(XXX+b)+a)+c)+=1;
				*(*(*(XXX+c)+b)+a)+=1;
				*(*(*(XXX+a)+c)+b)+=1;
				}
				}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N123_front(int row, int col, int mom, int u, int v, int w, int a, int b, int c, unsigned int ***XXX, Node ***nodeS, short int sx1, short int sy1, short int sz1, short int sx2, short int sy2, short int sz2){
	/*
	Funcion para contar los triangulos en tres 
	nodos con un puntos en N1, un punto en N2
	y un punto en N3.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	a, b, c => posición de N3.
	
	*/
	int i,j,k, n, m, l;
	float dx12,dy12,dz12,dx13,dy13,dz13,dx23,dy23,dz23;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	float d_max_c = size_box - d_max ;
	d_max_c *= d_max_c;
	
	for (i=0; i<nodeS[row][col][mom].len; ++i){
		// 1er punto en N1
		x1 = nodeS[row][col][mom].elements[i].x;
		y1 = nodeS[row][col][mom].elements[i].y;
		z1 = nodeS[row][col][mom].elements[i].z;
		for (j=0; j<nodeS[u][v][w].len; ++j){
			// 2do punto en N2
			x2 = nodeS[u][v][w].elements[j].x;
			y2 = nodeS[u][v][w].elements[j].y;
			z2 = nodeS[u][v][w].elements[j].z;
			dx12 = (size_box*(sx1)) - fabs(x2-x1);
			dy12 = (size_box*(sy1)) - fabs(y2-y1);
			dz12 = (size_box*(sz1)) - fabs(z2-z1);
			d12 = dx12*dx12+dy12*dy12+dz12*dz12;
			if (d12<=dd_max){
			for (k=0; k<nodeS[a][b][c].len; ++k){
				// 3er punto en N3 
				x3 = nodeS[a][b][c].elements[k].x;
				y3 = nodeS[a][b][c].elements[k].y;
				z3 = nodeS[a][b][c].elements[k].z;
				dx13 = (size_box*(sx2)) - fabs(x3-x1);
				dy13 = (size_box*(sy2)) - fabs(y3-y1);
				dz13 = (size_box*(sz2)) - fabs(z3-z1);
				d13 = dx13*dx13+dy13*dy13+dz13*dz13;
				if (d13<=dd_max){
				
				dx23 = x3-x2;
				dx23 *= dx23;
				if ( d_max_c <= dx23 ){
					dx23 = size_box - sqrt(dx23);
					dx23 *= dx23;
				}
				
				dy23 = y3-y2;
				dy23 *= dy23;
				if ( d_max_c <= dy23 ){
					dy23 = size_box - sqrt(dy23);
					dy23 *= dy23;
				}
				
				dz23 = z3-z2;
				dz23 *= dz23;
				if ( d_max_c <= dz23 ){
					dz23 = size_box - sqrt(dz23);
					dz23 *= dz23;
				}
				
				d23 = dx23+dy23+dz23;
				if (d23<=dd_max){
				n = (int)(sqrt(d12)*ds);
				m = (int)(sqrt(d13)*ds);
				l = (int)(sqrt(d23)*ds);
				*(*(*(XXX+n)+m)+l)+=1;
				*(*(*(XXX+l)+n)+m)+=1;
				*(*(*(XXX+m)+l)+n)+=1;
				*(*(*(XXX+m)+n)+l)+=1;
				*(*(*(XXX+l)+m)+n)+=1;
				*(*(*(XXX+n)+l)+m)+=1;
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


