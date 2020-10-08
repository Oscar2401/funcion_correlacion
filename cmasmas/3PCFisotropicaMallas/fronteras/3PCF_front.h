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
		int partitions;
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
			dd_max = d_max*d_max;
			front = size_box - d_max;
			corr = size_node*sqrt(3);
			ds = ((float)(bn))/d_max;
			ddmax_nod = (d_max+corr)*(d_max+corr);
			partitions = (int)(ceil(size_box/size_node));
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
		void count_3_N123(int, int, int, int, int, int, int, int, int, unsigned int ***, Node ***);
		void symmetrize(unsigned int ***);
		void symmetrize_analitic(float ***);
		void BPC(unsigned int ***, Point3D *);
		
		void make_histo_RRR(float ***, float ***, Node ***);
		void make_histo_DDR(float ***, float ***, Node ***);
		void make_histoXX(unsigned int *, float *, Node ***, int);
		void histo_front_XX(unsigned int *, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int, float);
		
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
	int i, j, k, row, col, mom, u, v, w, a ,b, c;
	float dis, dis_nod, dis_nod2, dis_nod3;
	float x1N, y1N, z1N, x2N, y2N, z2N, x3N, y3N, z3N;
	float x, y, z;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod, dx_nod2, dy_nod2, dz_nod2, dx_nod3, dy_nod3, dz_nod3;
	bool con_x, con_y, con_z;
	
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
	
	//===============================================================================
	// Condiciones Preriódicas de frontera
	//===============================================================================

	BPC(XXX,dataD);
	
	//================================
	// Simetrización:
	//================================
	
	symmetrize(XXX); 	
}
//=================================================================== 
void NODE3P::make_histo_RRR(float ***XXY, float ***XXX, Node ***nodeX){
	
	//======================================
	// Histograma RRR (ANALITICA)
	//======================================
	
	int i, j, k, u, v, w, a, b, c;
	int v_in;
	
	float dr = d_max/(float)bn;
	float dr2 = dr/2, r1, r2, r3;
	float ri, rj, rk;
	
	float V = size_box*size_box*size_box;
	float beta1 = n_pts*n_pts/V;
	float beta2 = n_pts/V;
	float alph = 8*(3.14159265359)*(3.14159265359)*(beta1*beta2)*dr*dr*dr;
	
	float bn_ref = 200; 
	float dr_ref = dr/bn_ref;
	float dr_ref2 = dr_ref/2;
	float ru, rv, rw;
	
	float alph_ref = 8*(3.14159265359)*(3.14159265359)*(beta1*beta2)*dr_ref*dr_ref*dr_ref;
	
	float f_av, con;
	
	for(i=0; i<bn; ++i) {
	ri = (float)(i)*dr;
	for(j=i; j<bn; ++j) {
	rj = (float)(j)*dr;
	for(k=j; k<bn; ++k) {
	rk = (float)(k)*dr;
		// Checamos vertices del cubo para hacer 
		// o no refinamiento  
		
		v_in = 0;
		
		for (a = 0; a < 2; ++a){
		r1 = ri + (float)(a)*dr;
		for (b = 0; b < 2; ++b){
		r2 = rj + (float)(b)*dr;
		for (c = 0; c < 2; ++c){
		r3 = rk + (float)(c)*dr;	
			if (r1 + r2 >= r3 && r1 + r3 >= r2 && r2 + r3 >= r1) v_in++; 
		}
		}
		}
		if (v_in==8) *(*(*(XXX+i)+j)+k) += alph*(ri+dr2)*(rj+dr2)*(rk+dr2);
		else if (v_in != 8 && v_in != 0){
			
			f_av = 0;
			con = 0;
			
			for(u=0; u<bn_ref; ++u) {
			ru = ri + (float)(u)*dr_ref;
			for(v=0; v<bn_ref; ++v) {
			rv = rj + (float)(v)*dr_ref;
			for(w=0; w<bn_ref; ++w) {
			rw = rk + (float)(w)*dr_ref;
			
				v_in = 0;
		
				for (a = 0; a < 2; ++a){
				r1 = ru + (float)(a)*dr_ref;
				for (b = 0; b < 2; ++b){
				r2 = rv + (float)(b)*dr_ref;
				for (c = 0; c < 2; ++c){
				r3 = rw + (float)(c)*dr_ref;	
					if (r1 + r2 >= r3 && r1 + r3 >= r2 && r2 + r3 >= r1) v_in++;
				}
				}
				}
				if (v_in==8){
					f_av += (ru + dr_ref2)*(rv + dr_ref2)*(rw + dr_ref2);
					con++;
				}
			}
			}
			}
			if (con!=0) *(*(*(XXX+i)+j)+k) += alph_ref*f_av;
		}
	}
	}
	}
	std::cout << " Termine HITOGRAMA RRR:" << std::endl;
	
	make_histo_DDR(XXY, XXX, nodeX);
	
	symmetrize_analitic(XXX);
	symmetrize_analitic(XXY); 
	
}
//=================================================================== 
void NODE3P::make_histo_DDR(float ***XXY, float ***XXX, Node ***nodeX){
	
	//======================================
	// Histograma DDR (ANALITICA)
	//======================================
	
	int i, j, k;
	
	float dr = d_max/(float)bn;
	float dr2 = dr/2, r1, r2, r3;
	float ri, rj, rk;
	
	float f_av;
	
	// Hacemos la funcion de correlación con 
	// ptt-veces el número de bins.
	
	int ptt = 100;
	int bins = ptt*bn;
	float dr_ptt = d_max/(float)bins;
	unsigned int  *DD;
	float *RR;
	DD = new unsigned int[bins];
	RR = new float[bins];
	for (i = 0; i < bins; i++){
		*(DD+i) = 0; 
		*(RR+i) = 0.0;
	}
	
	// Hacemos los histogramas DD y RR
	make_histoXX(DD, RR, nodeX, bins);
	
	// iniciamos un arreglo para la funcion f_averrage
	float *ff_av;
	ff_av = new float[bn];
	for (i = 0; i < bn; i++){
		*(ff_av+i) = 0;
	}
	
	// Calculamos las f_averrage
	for(i=0; i<bn; ++i){
		f_av = 0;
		float ri = (float)(i)*dr;
		for(j = 0; j<ptt; ++j){
			float rj = ((float)(j)+0.5)*dr_ptt;
			f_av += (ri + rj)*(((float)(*(DD+(i*ptt)+j)) /(*(RR+(i*ptt)+j))) - 1);
		}
		*(ff_av+i) = f_av/ptt;
	}
	
	// Calculamos el histograma DDR 
	for(i=0; i<bn; ++i) {
	ri = ((float)(i)+0.5)*dr;
	for(j=i; j<bn; ++j) {
	ri = ((float)(j)+0.5)*dr;
	for(k=j; k<bn; ++k) {
	ri = ((float)(k)+0.5)*dr;
		
		*(*(*(XXY+i)+j)+k) = *(*(*(XXX+i)+j)+k)*(1 + (*(ff_av+i)/(3*ri)) + (*(ff_av+j)/(3*rj)) + (*(ff_av+k)/(3*rk)));
		
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
				nod[row][col][mom].elements = new Point3D[0];
			}
		}
	}
	// Llenamos los nodos con los puntos de dat:
	for (i=0; i<n_pts; i++){
		row = (int)(dat[i].x/size_node);
        	col = (int)(dat[i].y/size_node);
        	mom = (int)(dat[i].z/size_node);
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
void NODE3P::symmetrize_analitic(float ***XXX){
	int i,j,k;
	float elem;
	for (i=0; i<bn; i++){
	for (j=0; j<bn; j++){
	for (k=0; k<bn; k++){
		elem = XXX[i][j][k];
		//XXX[i][j][k] = elem;
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
	int i,j,k;
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
void NODE3P::BPC(unsigned int ***XXX, Point3D *data){
	
	int a ,b, c;
	float d12, d13, d23;
	float x1, y1, z1, x2, y2, z2, x3, y3, z3;
	float dx, dy, dz;
	bool sx1, sy1, sz1, sx2, sy2, sz2;
	
	float d_max_f = size_box - d_max;
	float _d_max_f =  -1*d_max_f;

	for (a=0; a<n_pts-2; ++a){
		x1 = data[a].x;
		y1 = data[a].y;
		z1 = data[a].z;
		for (b=a+1; b<n_pts-1; ++b){
			x2 = data[b].x;
			y2 = data[b].y;
			z2 = data[b].z;
			
			dx = x2 - x1;
			sx1 = false;
			if ( d_max_f < dx){
				dx = size_box - dx;
				x2 = x1 - dx;
				sx1 = true;
			}
			else if (_d_max_f > dx){
				dx = size_box + dx;
				x2 = x1 + dx;
				sx1 = true;
			}
			
			if(dx < d_max){
			
			dy = y2 - y1;
			sy1 = false;
			if ( d_max_f < dy){
				dy = size_box - dy;
				y2 = y1 - dy;
				sy1 = true;
			}
			else if (_d_max_f > dy){
				dy = size_box + dy;
				y2 = y1 + dy;
				sy1 = true;
			}
			
			if(dy < d_max){
			
			dz = z2 - z1;
			sz1 = false;
			if ( d_max_f < dz){
				dz = size_box - dz;
				z2 = z1 - dz;
				sz1 = true;
			}
			else if ( _d_max_f > dz){
				dz = size_box + dz;
				z2 = z1 + dz;
				sz1 = true;
			}
			
			if(dz < d_max){
			
			d12 = dx*dx + dy*dy + dz*dz; 
			
			if (d12 < dd_max){
			
			for (c=b+1; c<n_pts; ++c){
				x3 = data[c].x;
				y3 = data[c].y;
				z3 = data[c].z;
				
				dx = x3 - x1;
				sx2 = false;
				if ( d_max_f < dx){
					dx = size_box - dx;
					x3 = x1 - dx;
					sx2 = true;
				}
				else if ( _d_max_f > dx){
					dx = size_box + dx;
					x3 = x1 + dx;
					sx2 = true;
				}
				
				dy = y3 - y1;
				sy2 = false;
				if ( d_max_f < dy){
					dy = size_box - dy;
					y3 = y1 - dy;
					sy2 = true;
				}
				else if ( _d_max_f > dy){
					dy = size_box + dy;
					y3 = y1 + dy;
					sy2 = true;
				}
				
				dz = z3 - z1;
				sz2 = false;
				if ( d_max_f < dz){
					dz = size_box - dz;
					z3 = z1 - dz;
					sz2 = true;
				}
				else if ( _d_max_f > dz){
					dz = size_box + dz;
					z3 = z1 + dz;
					sz2 = true;
				}
				
				d13 = dx*dx + dy*dy + dz*dz; 
				
				if (d13 < dd_max){
				
				if(sx2 || sy2 || sz2 || sx1 || sy1 || sz1){
				
					dx = x3 - x2;
					dy = y3 - y2;
					dz = z3 - z2;
					
					d23 = dx*dx + dy*dy + dz*dz; 
					
					if (d23 < dd_max){
						*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=1;
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

void NODE3P::make_histoXX(unsigned int *XX, float *YY, Node ***nodeX, int bins){
	/*
	Función para crear los histogramas DD y RR.
	
	Argumentos
	DD: arreglo donde se creará el histograma DD.
	RR: arreglo donde se creará el histograma RR.
	
	*/
	//Variables compartidas en hilos: 
	int i, j, row, col, mom, u, v, w;
	float dis, dis_nod;
	float x1D, y1D, z1D, x2D, y2D, z2D;
	float x, y, z;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod;
	bool con_x, con_y, con_z;
	float d_max_pm = d_max + size_node/2, front_pm = front - size_node/2;
	
	float ds_new = ((float)(bins))/d_max;
	
	for (row = 0; row < partitions; ++row){
	x1D = nodeX[row][0][0].nodepos.x;
	for (col = 0; col < partitions; ++col){
	y1D = nodeX[row][col][0].nodepos.y;
	for (mom = 0; mom < partitions; ++mom){
	z1D = nodeX[row][col][mom].nodepos.z;			
		//==================================================
		// Distancias entre puntos del mismo nodo:
		//==================================================
		for ( i= 0; i <nodeX[row][col][mom].len - 1; ++i){
			x = nodeX[row][col][mom].elements[i].x;
			y = nodeX[row][col][mom].elements[i].y;
			z = nodeX[row][col][mom].elements[i].z;
			for ( j = i+1; j < nodeX[row][col][mom].len; ++j){
				dx = x-nodeX[row][col][mom].elements[j].x;
				dy = y-nodeX[row][col][mom].elements[j].y;
				dz = z-nodeX[row][col][mom].elements[j].z;
				dis = dx*dx+dy*dy+dz*dz;
				if (dis <= dd_max){
				*(XX + (int)(sqrt(dis)*ds_new)) += 2;
				}
			}
		}
		//==================================================
		// Distancias entre puntos del diferente nodo:
		//==================================================
		u = row;
		v = col;
		//=========================
		// N2 movil en Z
		//=========================
		for (w=mom+1;  w<partitions ; ++w){	
			z2D = nodeX[u][v][w].nodepos.z;
			dz_nod = z1D-z2D;
			dis_nod = dz_nod*dz_nod;
			if (dis_nod <= ddmax_nod){
			for ( i = 0; i < nodeX[row][col][mom].len; ++i){
				x = nodeX[row][col][mom].elements[i].x;
				y = nodeX[row][col][mom].elements[i].y;
				z = nodeX[row][col][mom].elements[i].z;
				for ( j = 0; j < nodeX[u][v][w].len; ++j){
					dx = x-nodeX[u][v][w].elements[j].x;
					dy = y-nodeX[u][v][w].elements[j].y;
					dz = z-nodeX[u][v][w].elements[j].z;
					dis = dx*dx+dy*dy+dz*dz;
					if (dis <= dd_max){
					*(XX + (int)(sqrt(dis)*ds_new)) += 2;
					}
					}
				}
			}
			// Distacia de los puntos frontera XX
			//=======================================
			//Condiciones de nodos en frontera:
			con_z = ((z1D<=d_max_pm)&&(z2D>=front_pm))||((z2D<=d_max_pm)&&(z1D>=front_pm));
			if(con_z){
			histo_front_XX(XX,nodeX,dis_nod,0.0,0.0,fabs(dz_nod),false,false,con_z,row,col,mom,u,v,w,ds_new);
			}
		}
		//=========================
		// N2 movil en ZY
		//=========================
		for (v = col + 1; v < partitions ; ++v){
			y2D = nodeX[u][v][0].nodepos.y;
			dy_nod = y1D-y2D;
			for (w = 0; w < partitions ; ++w){		
				z2D = nodeX[u][v][w].nodepos.z;
				dz_nod = z1D-z2D;
				dis_nod = dy_nod*dy_nod + dz_nod*dz_nod;
				if (dis_nod <= ddmax_nod){
				for ( i = 0; i < nodeX[row][col][mom].len; ++i){
					x = nodeX[row][col][mom].elements[i].x;
						y = nodeX[row][col][mom].elements[i].y;
						z = nodeX[row][col][mom].elements[i].z;
						for ( j = 0; j < nodeX[u][v][w].len; ++j){	
							dx =  x-nodeX[u][v][w].elements[j].x;
							dy =  y-nodeX[u][v][w].elements[j].y;
							dz =  z-nodeX[u][v][w].elements[j].z;
							dis = dx*dx+dy*dy+dz*dz;
							if (dis <= dd_max){
							*(XX + (int)(sqrt(dis)*ds_new)) += 2;
							}
						}
				}
				}
				// Distacia de los puntos frontera
				//=======================================
				
				//Condiciones de nodos en frontera:
				con_y = ((y1D<=d_max_pm)&&(y2D>=front_pm))||((y2D<=d_max_pm)&&(y1D>=front_pm));
				con_z = ((z1D<=d_max_pm)&&(z2D>=front_pm))||((z2D<=d_max_pm)&&(z1D>=front_pm));
				if(con_y || con_z){ 
				histo_front_XX(XX,nodeX,dis_nod,0.0,fabs(dy_nod),fabs(dz_nod),false,con_y,con_z,row,col,mom,u,v,w,ds_new);
				}
			}
		}
		//=========================
		// N2 movil en ZYX
		//=========================
		for ( u = row + 1; u < partitions; ++u){
			x2D = nodeX[u][0][0].nodepos.x;
			dx_nod = x1D-x2D;
			for ( v = 0; v < partitions; ++v){
				y2D = nodeX[u][v][0].nodepos.y;
				dy_nod = y1D-y2D;
				for ( w = 0; w < partitions; ++w){
					z2D = nodeX[u][v][w].nodepos.z;
					dz_nod = z1D-z2D;
					dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
					if (dis_nod <= ddmax_nod){
					for ( i = 0; i < nodeX[row][col][mom].len; ++i){
						x = nodeX[row][col][mom].elements[i].x;
						y = nodeX[row][col][mom].elements[i].y;
						z = nodeX[row][col][mom].elements[i].z;
						for ( j = 0; j < nodeX[u][v][w].len; ++j){	
							dx = x-nodeX[u][v][w].elements[j].x;
							dy = y-nodeX[u][v][w].elements[j].y;
							dz = z-nodeX[u][v][w].elements[j].z;
							dis = dx*dx + dy*dy + dz*dz;
							if (dis <= dd_max){
							*(XX + (int)(sqrt(dis)*ds_new)) += 2;
							}
						}
					}
					}
					// Distacia de los puntos frontera
					//=======================================
			
					//Condiciones de nodos en frontera:
					con_x = ((x1D<=d_max_pm)&&(x2D>=front_pm))||((x2D<=d_max_pm)&&(x1D>=front_pm));
					con_y = ((y1D<=d_max_pm)&&(y2D>=front_pm))||((y2D<=d_max_pm)&&(y1D>=front_pm));
					con_z = ((z1D<=d_max_pm)&&(z2D>=front_pm))||((z2D<=d_max_pm)&&(z1D>=front_pm));
			if(con_x || con_y || con_z){
			histo_front_XX(XX,nodeX,dis_nod,fabs(dx_nod),fabs(dy_nod),fabs(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w,ds_new);
				}	
				}	
			}
		}
	}
	}
	}
	
	//======================================
	// Histograma RR (ANALITICA)
	//======================================
	
	float dr = (d_max/bins);
	float V = size_box*size_box*size_box;
	float beta1 = n_pts*n_pts/V;
	float alph = 4*(3.14159265359)*(beta1)*dr*dr*dr/3;
	float r1, r2;
	for(int a=0; a<bins; ++a) {
		r2 = (float) a;
		r1 = r2+1;
        	*(YY+a) += alph*((r1*r1*r1)-(r2*r2*r2));
	}
}

//=================================================================== 
void NODE3P::histo_front_XX(unsigned int *PP, Node ***dat, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int _row, int _col, int _mom, int _u, int _v, int _w, float ds_new){
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
						*(PP + (int)(sqrt(_dis)*ds_new)) += 2;
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
						*(PP + (int)(sqrt(_dis)*ds_new)) += 2;
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
						*(PP + (int)(sqrt(_dis)*ds_new)) += 2;
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
						*(PP + (int)(sqrt(_dis)*ds_new)) += 2;
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
						*(PP + (int)(sqrt(_dis)*ds_new)) += 2;
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
						*(PP + (int)(sqrt(_dis)*ds_new)) += 2;
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
						*(PP + (int)(sqrt(_dis)*ds_new)) += 2;
					}
				}
			}
		}
	}
}

//=================================================================== 
NODE3P::~NODE3P(){	
}


