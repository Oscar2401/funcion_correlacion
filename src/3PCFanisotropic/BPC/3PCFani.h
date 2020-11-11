#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <omp.h>

struct Point3D{
	float x;
	float y; 
	float z;
	short int fx;
	short int fy;
	short int fz;
};

struct PointW3D{
	float x;
	float y; 
	float z;
	float w;
	short int fx;
	short int fy;
	short int fz;
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
		float size_box_2;
		
	private: 
		void make_nodos(Node ***, PointW3D *);
		void add(PointW3D *&, int&, float, float, float, float, short int, short int, short int);
	
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
			size_box_2 =  size_box/2;
			
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
		void make_histoXXX(float *****, Node ***);
		void count_3_N111(int, int, int, float *****, Node ***);
		void count_3_N112(int, int, int, int, int, int, float *****, Node ***);
		void count_3_N123(int, int, int, int, int, int, int, int, int, float *****, Node ***);
		
		void front_node_112(int, int, int, int, int, int, bool, bool, bool, float *****, Node ***);
		void front_112(int,int,int,int,int,int,short int,short int,short int,float *****, Node ***);
		void front_node_123(int,int,int,int,int,int,int,int,int,bool,bool,bool,float *****, Node ***);
		void front_123(int,int,int,int,int,int,int,int,int,short int,short int,short int,short int,short int,short int,short int,short int,short int,float *****, Node ***);
		
		void make_histoXXY(float *****, Node ***, Node ***, PointW3D *, PointW3D *);
		void count_3_N112_xxy(int, int, int, int, int, int, float *****, Node ***, Node ***);
		void count_3_N123_xxy(int, int, int, int, int, int, int, int, int, float *****, Node ***, Node ***);
		
		void front_node_112_xxy(int, int, int, int, int, int, bool, bool, bool, float *****, Node ***, Node ***);
		void front_112_xxy(int,int,int,int,int,int,short int,short int,short int,float *****, Node ***, Node ***);
		void front_node_123_xxy(int,int,int,int,int,int,int,int,int,bool,bool,bool,float *****, Node ***, Node ***);
		void front_123_xxy(int,int,int,int,int,int,int,int,int,short int,short int,short int,short int,short int,short int,short int,short int,short int,float *****, Node ***, Node ***);
		
		~NODE3P();
};

//=================================================================== 
//==================== Funciones ==================================== 
//=================================================================== 
void NODE3P::make_nodos(Node ***nod, PointW3D *dat){
	/*
	Función para crear los nodos con los datos y puntos random
	
	Argumentos
	nod: arreglo donde se crean los nodos.
	dat: datos a dividir en nodos.
	
	*/
	int i, row, col, mom;
	float d_max_pm = d_max + size_node/2;
	float front_pm = front - size_node/2;
	float p_med = size_node/2;
	float posx, posy, posz;
	
	// Inicializamos los nodos vacíos:
	for (row=0; row<partitions; ++row){
	for (col=0; col<partitions; ++col){
	for (mom=0; mom<partitions; ++mom){
		posx = ((float)(row)*(size_node))+p_med;
		posy = ((float)(col)*(size_node))+p_med;
		posz = ((float)(mom)*(size_node))+p_med;
		
		nod[row][col][mom].nodepos.x = posx;
		nod[row][col][mom].nodepos.y = posy;
		nod[row][col][mom].nodepos.z = posz;
		
		// Vemos si el nodo esta en la frontera
		// frontera x:
		if(posx<=d_max_pm) nod[row][col][mom].nodepos.fx = -1;
		else if(posx>=front_pm) nod[row][col][mom].nodepos.fx = 1;
		else nod[row][col][mom].nodepos.fx = 0;
		// frontera y:
		if(posy<=d_max_pm) nod[row][col][mom].nodepos.fy = -1;
		else if(posy>=front_pm) nod[row][col][mom].nodepos.fy = 1;
		else nod[row][col][mom].nodepos.fy = 0;
		// frontera z:
		if(posz<=d_max_pm) nod[row][col][mom].nodepos.fz = -1;
		else if(posz>=front_pm) nod[row][col][mom].nodepos.fz = 1;
		else nod[row][col][mom].nodepos.fz = 0;
		
		nod[row][col][mom].len = 0;
		nod[row][col][mom].elements = new PointW3D[0];
	}
	}
	}
	// Llenamos los nodos con los puntos de dat:
	for (i=0; i<n_pts; ++i){
		row = (int)(dat[i].x/size_node);
        	col = (int)(dat[i].y/size_node);
        	mom = (int)(dat[i].z/size_node);
		add(nod[row][col][mom].elements,nod[row][col][mom].len,dat[i].x,dat[i].y,dat[i].z,dat[i].w,dat[i].fx,dat[i].fy,dat[i].fz);
	}
}
//=================================================================== 
void NODE3P::add(PointW3D *&array, int &lon, float _x, float _y, float _z, float _w, short int _fx, short int _fy, short int _fz){
	lon++;
	PointW3D *array_aux = new PointW3D[lon];
	for (int i=0; i<lon-1; ++i){
		array_aux[i].x = array[i].x;
		array_aux[i].y = array[i].y;
		array_aux[i].z = array[i].z;
		array_aux[i].w = array[i].w;
		array_aux[i].fx = array[i].fx;
		array_aux[i].fy = array[i].fy;
		array_aux[i].fz = array[i].fz;
	}
	delete[] array;
	array = array_aux;
	array[lon-1].x = _x;
	array[lon-1].y = _y; 
	array[lon-1].z = _z;
	array[lon-1].w = _w; 
	array[lon-1].fx = _fx;
	array[lon-1].fy = _fy; 
	array[lon-1].fz = _fz;
}
//=================================================================== 
void NODE3P::make_histoXXX(float *****XXX, Node ***nodeX){
	/*
	Función para crear los histogramas DDD y RRR.
	
	Argumentos
	XXX: arreglo donde se creará el histograma DDD y RRR.
	nodeX: malla de datos.
	
	*/ 
	int i, j, k, row, col, mom, u, v, w, a ,b, c;
	float dis, dis_nod, dis_nod2, dis_nod3;
	float x1N, y1N, z1N, x2N, y2N, z2N, x3N, y3N, z3N;
	float x, y, z;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod, dx_nod2, dy_nod2, dz_nod2, dx_nod3, dy_nod3, dz_nod3;
	
	float fx1N, fy1N, fz1N;
	float fx2N, fy2N, fz2N;
	float fx3N, fy3N, fz3N;
	
	bool con1_x, con1_y, con1_z;
	bool con2_x, con2_y, con2_z;
	bool con3_x, con3_y, con3_z;
	bool con12_x, con12_y, con12_z;
	bool con13_x, con13_y, con13_z;
	bool conx, cony, conz, conx_, cony_, conz_;
	
	// x1N, y1N, z1N => Nodo pivote
	for (row=0; row<partitions; ++row){
	x1N = nodeX[row][0][0].nodepos.x;
	fx1N = nodeX[row][0][0].nodepos.fx;
	con1_x = fx1N != 0;
	for (col=0; col<partitions; ++col){
	y1N = nodeX[row][col][0].nodepos.y;
	fy1N = nodeX[row][col][0].nodepos.fy;
	con1_y = fy1N != 0;
	for (mom=0; mom<partitions; ++mom){
	fz1N = nodeX[row][col][mom].nodepos.fz;
	con1_z = fz1N != 0;
		
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
	
	//  =========================  Frontera  =========================
	
	if (con1_x||con1_y||con1_z){  // si todo es falso, entonces N1 no está en ninguna frontera
		
		u = row;
		v = col;
		
		fx2N = nodeX[u][0][0].nodepos.fx;
		con2_x = fx2N != 0;
		con12_x = fx1N == fx2N;
		conx_ = !con12_x && con1_x && con2_x;
		
		fy2N = nodeX[u][v][0].nodepos.fy;
		con2_y = fy2N != 0;
		con12_y = fy1N == fy2N;
		cony_ = !con12_y && con1_y && con2_y;
		
		//=======================
		// Nodo 2 movil en Z:
		//=======================
		for (w=mom+1;  w<partitions; ++w){	
		fz2N = nodeX[u][v][w].nodepos.fz;
		con2_z = fz2N != 0;
		con12_z = fz1N == fz2N;
		conz_ = !con12_z && con1_z && con2_z;
		
		if(con2_x || con2_y || con2_z){ // si todo es falso, entonces N2 no está en ninguna frontera
			
			// =======================================================================
			if (conx_||cony_||conz_) front_node_112(row,col,mom,u,v,w,conx_,cony_,conz_,XXX,nodeX);
			// =======================================================================
			
			a = u;
			b = v;
			
			fx3N = nodeX[a][0][0].nodepos.fx;
			con3_x = fx3N != 0;
			con13_x = fx1N == fx3N;
			conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
			
			fy3N = nodeX[a][b][0].nodepos.fy;
			con3_y = fy3N != 0;
			con13_y = fy1N == fy3N;
			cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
			
			for (c=w+1; c<partitions; ++c){
			fz3N = nodeX[a][b][c].nodepos.fz;
			con3_z = fz3N != 0;
			con13_z = fz1N == fz3N;
			conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
			if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXX,nodeX);
			
			}
			
			for (b=v+1; b<partitions; ++b){
			fy3N = nodeX[a][b][0].nodepos.fy;
			con3_y = fy3N != 0;
			con13_y = fy1N == fy3N;
			cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
				for (c=0; c<partitions; ++c){
				fz3N = nodeX[a][b][c].nodepos.fz;
				con3_z = fz3N != 0;
				con13_z = fz1N == fz3N;
				conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
				if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXX,nodeX);
				
				}
			}
			
			for (a=u+1; a<partitions; ++a){
			fx3N = nodeX[a][0][0].nodepos.fx;
			con3_x = fx3N != 0;
			con13_x = fx1N == fx3N;
			conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
				for (b=0; b<partitions; ++b){
				fy3N = nodeX[a][b][0].nodepos.fy;
				con3_y = fy3N != 0;
				con13_y = fy1N == fy3N;
				cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
					for (c=0;  c<partitions; ++c){
					fz3N = nodeX[a][b][c].nodepos.fz;
					con3_z = fz3N != 0;
					con13_z = fz1N == fz3N;
					conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
						
					if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXX,nodeX);
					
					}
				}
			}
		}
		}
		//=======================
		// Nodo 2 movil en ZY:
		//=======================
		for (v=col+1; v<partitions ; ++v){
		fy2N = nodeX[u][v][0].nodepos.fy;
		con2_y = fy2N != 0;
		con12_y = fy1N == fy2N;
		cony_ = !con12_y && con1_y && con2_y;
			for (w=0; w<partitions ; ++w){
			fz2N = nodeX[u][v][w].nodepos.fz;
			con2_z = fz2N != 0;
			con12_z = fz1N == fz2N;
			conz_ = !con12_z && con1_z && con2_z;
			
			if(con2_x || con2_y || con2_z){
			
				// =======================================================================
				if (conx_||cony_||conz_) front_node_112(row,col,mom,u,v,w,conx_,cony_,conz_,XXX,nodeX);
				// =======================================================================
				
				a = u;
				b = v;
				
				fx3N = nodeX[a][0][0].nodepos.fx;
				con3_x = fx3N != 0;
				con13_x = fx1N == fx3N;
				conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
				
				fy3N = nodeX[a][b][0].nodepos.fy;
				con3_y = fy3N != 0;
				con13_y = fy1N == fy3N;
				cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
				
				for (c=w+1; c<partitions; ++c){
				fz3N = nodeX[a][b][c].nodepos.fz;
				con3_z = fz3N != 0;
				con13_z = fz1N == fz3N;
				conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
				if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXX,nodeX);
				
				}
				
				for (b=v+1; b<partitions; ++b){
				fy3N = nodeX[a][b][0].nodepos.fy;
				con3_y = fy3N != 0;
				con13_y = fy1N == fy3N;
				cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
					for (c=0; c<partitions; ++c){
					fz3N = nodeX[a][b][c].nodepos.fz;
					con3_z = fz3N != 0;
					con13_z = fz1N == fz3N;
					conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
					if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXX,nodeX);
					
					}
				}
				
				for (a=u+1; a<partitions; ++a){
				fx3N = nodeX[a][0][0].nodepos.fx;
				con3_x = fx3N != 0;
				con13_x = fx1N == fx3N;
				conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
					for (b=0; b<partitions; ++b){
					fy3N = nodeX[a][b][0].nodepos.fy;
					con3_y = fy3N != 0;
					con13_y = fy1N == fy3N;
					cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
						for (c=0;  c<partitions; ++c){
						fz3N = nodeX[a][b][c].nodepos.fz;
						con3_z = fz3N != 0;
						con13_z = fz1N == fz3N;
						conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
						if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXX,nodeX);
						
						}
					}
				}
			}
			}
		}
		//=======================
		// Nodo 2 movil en ZYX:
		//=======================
		for (u=row+1; u<partitions ; ++u){
		fx2N = nodeX[u][0][0].nodepos.fx;
		con2_x = fx2N != 0;
		con12_x = fx1N == fx2N;
		conx_ = !con12_x && con1_x && con2_x;
			for (v=0; v<partitions ; ++v){
			fy2N = nodeX[u][v][0].nodepos.fy;
			con2_y = fy2N != 0;
			con12_y = fy1N == fy2N;
			cony_ = !con12_y && con1_y && con2_y;
				for (w=0; w<partitions ; ++w){
				fz2N = nodeX[u][v][w].nodepos.fz;
				con2_z = fz2N != 0;
				con12_z = fz1N == fz2N;
				conz_ = !con12_z && con1_z && con2_z;
				if(con2_x || con2_y || con2_z){
					
					// =======================================================================
					if (conx_||cony_||conz_) front_node_112(row,col,mom,u,v,w,conx_,cony_,conz_,XXX,nodeX);
					// =======================================================================
					
					a = u;
					b = v;
					
					fx3N = nodeX[a][0][0].nodepos.fx;
					con3_x = fx3N != 0;
					con13_x = fx1N == fx3N;
					conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
					
					fy3N = nodeX[a][b][0].nodepos.fy;
					con3_y = fy3N != 0;
					con13_y = fy1N == fy3N;
					cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
					
					for (c=w+1; c<partitions; ++c){
					fz3N = nodeX[a][b][c].nodepos.fz;
					con3_z = fz3N != 0;
					con13_z = fz1N == fz3N;
					conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
					
					if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXX,nodeX);
					
					}
					
					for (b=v+1; b<partitions; ++b){
					fy3N = nodeX[a][b][0].nodepos.fy;
					con3_y = fy3N != 0;
					con13_y = fy1N == fy3N;
					cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
						for (c=0; c<partitions; ++c){
						fz3N = nodeX[a][b][c].nodepos.fz;
						con3_z = fz3N != 0;
						con13_z = fz1N == fz3N;
						conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
						
						if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXX,nodeX);
						
						}
					}
					
					for (a=u+1; a<partitions; ++a){
					fx3N = nodeX[a][0][0].nodepos.fx;
					con3_x = fx3N != 0;
					con13_x = fx1N == fx3N;
					conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
						for (b=0; b<partitions; ++b){
						fy3N = nodeX[a][b][0].nodepos.fy;
						con3_y = fy3N != 0;
						con13_y = fy1N == fy3N;
						cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
							for (c=0;  c<partitions; ++c){
							fz3N = nodeX[a][b][c].nodepos.fz;
							con3_z = fz3N != 0;
							con13_z = fz1N == fz3N;
							conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
							if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXX,nodeX);
							
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
void NODE3P::count_3_N111(int row, int col, int mom, float *****XXX, Node ***nodeS){
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
void NODE3P::count_3_N112(int row, int col, int mom, int u, int v, int w, float *****XXX, Node ***nodeS){
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
void NODE3P::count_3_N123(int row, int col, int mom, int u, int v, int w, int a, int b, int c, float *****XXX, Node ***nodeS){
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
void NODE3P::front_node_112(int row, int col, int mom, int u, int v, int w, bool conx, bool cony, bool conz, float *****XXX, Node ***nodeS){
	
	float x1N = nodeS[row][0][0].nodepos.x;
	float y1N = nodeS[row][col][0].nodepos.y;
	float z1N = nodeS[row][col][mom].nodepos.z;
	float fx2N = nodeS[u][0][0].nodepos.fx, x2N = nodeS[u][0][0].nodepos.x;
	float fy2N = nodeS[u][v][0].nodepos.fy, y2N = nodeS[u][v][0].nodepos.y;
	float fz2N = nodeS[u][v][w].nodepos.fz, z2N = nodeS[u][v][w].nodepos.z;
	
	float dx, dy, dz;
	
	if (conx){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,fx2N,0,0,XXX,nodeS);
	}
	if (cony){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,0,fy2N,0,XXX,nodeS);
	}
	if (conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz =(z2N-(fz2N*size_box)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,0,0,fz2N,XXX,nodeS);
	}
	if (conx && cony){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,fx2N,fy2N,0,XXX,nodeS);
	}
	if (conx && conz){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = (z2N-fz2N*size_box) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,fx2N,0,fz2N,XXX,nodeS);
	}
	if (cony && conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = (z2N-fz2N*size_box) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,0,fy2N,fz2N,XXX,nodeS);
	}
	if (conx && cony && conz){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = (z2N-fz2N*size_box) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,fx2N,fy2N,fz2N,XXX,nodeS);
	}
}
//=================================================================== 
void NODE3P::front_112(int row,int col,int mom,int u,int v,int w,short int fx2,short int fy2,short int fz2,float *****XXX,Node ***nodeS){
	/*
	Funcion para contar los triangulos en dos 
	nodos con dos puntos en N1 y un punto en N2.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	
	*/
	int i,j,k;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3,W;
	
	int a_,b_,c_;
	int t,p,q, t_,p_,q_;
	float dx,dy,dz1,dz2,dz3;
	float cth1,cth2,cth3, cth1_,cth2_,cth3_;

	for (i=0; i<nodeS[row][col][mom].len; ++i){
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeS[u][v][w].len; ++j){
		x2 = nodeS[u][v][w].elements[j].x - fx2*size_box;
		y2 = nodeS[u][v][w].elements[j].y - fy2*size_box;
		z2 = nodeS[u][v][w].elements[j].z - fz2*size_box;
		w2 = nodeS[u][v][w].elements[j].w;
		dx = x2-x1;
		dy = y2-y1;
		dz1 = z2-z1;
		d12 = dx*dx+dy*dy+dz1*dz1;
		if (d12<dd_max){
		d12 = sqrt(d12);
			for (k=i+1; k<nodeS[row][col][mom].len; ++k){
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
				a_ = (int)(d12*ds);
				b_ = (int)(d13*ds);
				c_ = (int)(d23*ds);
				
				t = (int)(cth1*ds_th);
				p = (int)(cth2*ds_th);
				q = (int)(cth3*ds_th);
				
				t_ = (int)(cth1_*ds_th);
				p_ = (int)(cth2_*ds_th);
				q_ = (int)(cth3_*ds_th);
				
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
			x3 = nodeS[u][v][w].elements[k].x - fx2*size_box;
			y3 = nodeS[u][v][w].elements[k].y - fy2*size_box;
			z3 = nodeS[u][v][w].elements[k].z - fz2*size_box;
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
			if (d23<dd_max) {
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
				a_ = (int)(d12*ds);
				b_ = (int)(d13*ds);
				c_ = (int)(d23*ds);
				
				t = (int)(cth1*ds_th);
				p = (int)(cth2*ds_th);
				q = (int)(cth3*ds_th);
				
				t_ = (int)(cth1_*ds_th);
				p_ = (int)(cth2_*ds_th);
				q_ = (int)(cth3_*ds_th);
				
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
void NODE3P::front_node_123(int row, int col, int mom, int u, int v, int w, int a, int b, int c, bool conx, bool cony, bool conz, float *****XXX, Node ***nodeS){
	
	float fx1N = nodeS[row][0][0].nodepos.fx, x1N = nodeS[row][0][0].nodepos.x;
	float fy1N = nodeS[row][col][0].nodepos.fy, y1N = nodeS[row][col][0].nodepos.y;
	float fz1N = nodeS[row][col][mom].nodepos.fz, z1N = nodeS[row][col][mom].nodepos.z;
	float fx2N = nodeS[u][0][0].nodepos.fx, x2N = nodeS[u][0][0].nodepos.x;
	float fy2N = nodeS[u][v][0].nodepos.fy, y2N = nodeS[u][v][0].nodepos.y;
	float fz2N = nodeS[u][v][w].nodepos.fz, z2N = nodeS[u][v][w].nodepos.z;
	float fx3N = nodeS[a][0][0].nodepos.fx, x3N = nodeS[a][0][0].nodepos.x;
	float fy3N = nodeS[a][b][0].nodepos.fy, y3N = nodeS[a][b][0].nodepos.y;
	float fz3N = nodeS[a][b][c].nodepos.fz, z3N = nodeS[a][b][c].nodepos.z;
	
	float dx,dy,dz;
	
	if (conx){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = y3N-y1N;
			dy *= dy; 
			dz = z3N-z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = y3N-y2N;
				dy *= dy; 
				dz = z3N-z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,fx1N,0,0,fx2N,0,0,fx3N,0,0,XXX,nodeS);
				}
			}
		}
	
	}
	if (cony){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = x3N-x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = z3N-z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = x3N-x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = z3N-z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,0,fy1N,0,0,fy2N,0,0,fy3N,0,XXX,nodeS);
				}
			}
		}
	}
	if (conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz =(z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = x3N-x1N;
			dx *= dx;
			dy = y3N-y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = x3N-x2N;
				dx *= dx;
				dy = y3N-y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,0,0,fz1N,0,0,fz2N,0,0,fz3N,XXX,nodeS);
				}
			}
		}
	}
	if (conx && cony){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = z3N-z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = z3N-z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,fx1N,fy1N,0,fx2N,fy2N,0,fx3N,fy3N,0,XXX,nodeS);
				}
			}
		}
	}
	if (conx && conz){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = (z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = y3N-y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = y3N-y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,fx1N,0,fz1N,fx2N,0,fz2N,fx3N,0,fz3N,XXX,nodeS);
				}
			}
		}
	}
	if (cony && conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = (z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = x3N-x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = x3N-x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,0,fy1N,fz1N,0,fy2N,fz2N,0,fy3N,fz3N,XXX,nodeS);
				}
			}
		}
	}
	if (conx && cony && conz){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = (z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,fx1N,fy1N,fz1N,fx2N,fy2N,fz2N,fx3N,fy3N,fz3N,XXX,nodeS);
				}
			}
		}
	}	
}
//=================================================================== 
void NODE3P::front_123(int row,int col,int mom,int u,int v,int w,int a,int b,int c,short int fx1N,short int fy1N,short int fz1N,short int fx2N,short int fy2N,short int fz2N,short int fx3N,short int fy3N,short int fz3N,float *****XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en tres 
	nodos con un puntos en N1, un punto en N2
	y un punto en N3.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	a, b, c => posición de N3.
	
	*/
	
	int i,j,k;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3,W;
	
	int a_,b_,c_;
	int t,p,q, t_,p_,q_;
	float dx,dy,dz1,dz2,dz3;
	float cth1,cth2,cth3, cth1_,cth2_,cth3_;
	
	for (i=0; i<nodeS[row][col][mom].len; ++i){
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeS[u][v][w].len; ++j){
		x2 = nodeS[u][v][w].elements[j].x;
		y2 = nodeS[u][v][w].elements[j].y;
		z2 = nodeS[u][v][w].elements[j].z;
		w2 = nodeS[u][v][w].elements[j].w;
		dx = (x2-((fx2N-fx1N)*size_box_2))-x1;
		dy = (y2-((fy2N-fy1N)*size_box_2))-y1;
		dz1 = (z2-((fz2N-fz1N)*size_box_2))-z1;
		d12 = dx*dx+dy*dy+dz1*dz1;
		if (d12<dd_max){
		d12 = sqrt(d12);
			for (k=0; k<nodeS[a][b][c].len; ++k){
			x3 = nodeS[a][b][c].elements[k].x;
			y3 = nodeS[a][b][c].elements[k].y;
			z3 = nodeS[a][b][c].elements[k].z;
			w3 = nodeS[a][b][c].elements[k].w;
			dx = (x3-((fx3N-fx1N)*size_box_2))-x1;
			dy = (y3-((fy3N-fy1N)*size_box_2))-y1;
			dz2 = (z3-((fz3N-fz1N)*size_box_2))-z1;
			d13 = dx*dx+dy*dy+dz2*dz2;
			if (d13<dd_max){
			d13 = sqrt(d13);
			dx = (x3-((fx3N-fx2N)*size_box_2))-x2;
			dy = (y3-((fy3N-fy2N)*size_box_2))-y2;
			dz3 = (z3-((fz3N-fz2N)*size_box_2))-z2;
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
void NODE3P::make_histoXXY(float *****XXY, Node ***nodeX, Node ***nodeY, PointW3D *datx, PointW3D *daty){
	/*
	Función para crear los histogramas DDD y RRR.
	
	Argumentos
	DDD: arreglo donde se creará el histograma DDD.
	
	*/ 
	int i, j, k, row, col, mom, u, v, w, a, b, c;
	float dis, dis_nod, dis_nod2, dis_nod3;
	float x1N, y1N, z1N, x2N, y2N, z2N, x3N, y3N, z3N;
	float x, y, z;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod, dx_nod2, dy_nod2, dz_nod2, dx_nod3, dy_nod3, dz_nod3;
	
	float fx1N, fy1N, fz1N;
	float fx2N, fy2N, fz2N;
	float fx3N, fy3N, fz3N;
	
	bool con1_x, con1_y, con1_z;
	bool con2_x, con2_y, con2_z;
	bool con3_x, con3_y, con3_z;
	bool con12_x, con12_y, con12_z;
	bool con13_x, con13_y, con13_z;
	bool conx, cony, conz, conx_, cony_, conz_;
	
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
		for (v=0; v<partitions; ++v){
		y2N = nodeY[u][v][0].nodepos.y;
		dy_nod = y2N-y1N;
			for (w=0; w<partitions; ++w){
			z2N = nodeY[u][v][w].nodepos.z;
			dz_nod = z2N-z1N;
			dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
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
			dy_nod2 = y3N-y1N;
				for (c=mom+1;  c<partitions; ++c){	
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod2 = z3N-z1N;
				dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
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
					for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					dis_nod2 = dx_nod2*dx_nod2 + dy_nod2*dy_nod2 + dz_nod2*dz_nod2;
					if (dis_nod2 <= ddmax_nod){
					dy_nod3 = y3N-y2N;
					dz_nod3 = z3N-z2N;
					dis_nod3 = dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
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
	
	//  =========================  Frontera  =========================
	
	if (con1_x||con1_y||con1_z){  // si todo es falso, entonces N1 no está en ninguna frontera
		
		for (u=0; u<partitions ; ++u){
		fx2N = nodeY[u][0][0].nodepos.fx;
		con2_x = fx2N != 0;
		con12_x = fx1N == fx2N;
		conx_ = !con12_x && con1_x && con2_x;
			for (v=0; v<partitions ; ++v){
			fy2N = nodeY[u][v][0].nodepos.fy;
			con2_y = fy2N != 0;
			con12_y = fy1N == fy2N;
			cony_ = !con12_y && con1_y && con2_y;
				for (w=0; w<partitions ; ++w){
				fz2N = nodeY[u][v][w].nodepos.fz;
				con2_z = fz2N != 0;
				con12_z = fz1N == fz2N;
				conz_ = !con12_z && con1_z && con2_z;
				if(con2_x || con2_y || con2_z){
					
					// =======================================================================
					if (conx_||cony_||conz_) front_node_112_xxy(row,col,mom,u,v,w,conx_,cony_,conz_,XXY,nodeX,nodeY);
					// =======================================================================
					
					a = row;
					b = col;
					
					fx3N = nodeX[a][0][0].nodepos.fx;
					con3_x = fx3N != 0;
					con13_x = fx1N == fx3N;
					conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
					
					fy3N = nodeX[a][b][0].nodepos.fy;
					con3_y = fy3N != 0;
					con13_y = fy1N == fy3N;
					cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
					
					for (c=mom+1; c<partitions; ++c){
					fz3N = nodeX[a][b][c].nodepos.fz;
					con3_z = fz3N != 0;
					con13_z = fz1N == fz3N;
					conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
					
					if (conx||cony||conz){
						front_node_123_xxy(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXY,nodeX,nodeY);
					}
					}
					
					for (b=col+1; b<partitions; ++b){
					fy3N = nodeX[a][b][0].nodepos.fy;
					con3_y = fy3N != 0;
					con13_y = fy1N == fy3N;
					cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
						for (c=0; c<partitions; ++c){
						fz3N = nodeX[a][b][c].nodepos.fz;
						con3_z = fz3N != 0;
						con13_z = fz1N == fz3N;
						conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
						
						if (conx||cony||conz){
							front_node_123_xxy(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXY,nodeX,nodeY);
						}
						}
					}
					for (a=row+1; a<partitions; ++a){
					fx3N = nodeX[a][0][0].nodepos.fx;
					con3_x = fx3N != 0;
					con13_x = fx1N == fx3N;
					conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
						for (b=0; b<partitions; ++b){
						fy3N = nodeX[a][b][0].nodepos.fy;
						con3_y = fy3N != 0;
						con13_y = fy1N == fy3N;
						cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
							for (c=0;  c<partitions; ++c){
							fz3N = nodeX[a][b][c].nodepos.fz;
							con3_z = fz3N != 0;
							con13_z = fz1N == fz3N;
							conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
							if (conx||cony||conz){
								front_node_123_xxy(row,col,mom,u,v,w,a,b,c,conx,cony,conz,XXY,nodeX,nodeY);
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
void NODE3P::count_3_N112_xxy(int row, int col, int mom, int u, int v, int w, float *****XXY, Node ***nodeS, Node ***nodeT){
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
void NODE3P::count_3_N123_xxy(int row, int col, int mom, int u, int v, int w, int a, int b, int c, float *****XXY, Node ***nodeS, Node ***nodeT){
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
void NODE3P::front_node_112_xxy(int row, int col, int mom, int u, int v, int w, bool conx, bool cony, bool conz, float *****XXY, Node ***nodeS, Node ***nodeT){
	
	float x1N = nodeS[row][0][0].nodepos.x;
	float y1N = nodeS[row][col][0].nodepos.y;
	float z1N = nodeS[row][col][mom].nodepos.z;
	float fx2N = nodeT[u][0][0].nodepos.fx, x2N = nodeT[u][0][0].nodepos.x;
	float fy2N = nodeT[u][v][0].nodepos.fy, y2N = nodeT[u][v][0].nodepos.y;
	float fz2N = nodeT[u][v][w].nodepos.fz, z2N = nodeT[u][v][w].nodepos.z;
	
	float dx, dy, dz;
	
	if (conx){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112_xxy(row,col,mom,u,v,w,fx2N,0,0,XXY,nodeS,nodeT);
	}
	if (cony){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112_xxy(row,col,mom,u,v,w,0,fy2N,0,XXY,nodeS,nodeT);
	}
	if (conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz =(z2N-(fz2N*size_box)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112_xxy(row,col,mom,u,v,w,0,0,fz2N,XXY,nodeS,nodeT);
	}
	if (conx && cony){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112_xxy(row,col,mom,u,v,w,fx2N,fy2N,0,XXY,nodeS,nodeT);
	}
	if (conx && conz){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = (z2N-fz2N*size_box) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112_xxy(row,col,mom,u,v,w,fx2N,0,fz2N,XXY,nodeS,nodeT);
	}
	if (cony && conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = (z2N-fz2N*size_box) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112_xxy(row,col,mom,u,v,w,0,fy2N,fz2N,XXY,nodeS,nodeT);
	}
	if (conx && cony && conz){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = (z2N-fz2N*size_box) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112_xxy(row,col,mom,u,v,w,fx2N,fy2N,fz2N,XXY,nodeS,nodeT);
	}
}
//=================================================================== 
void NODE3P::front_112_xxy(int row,int col,int mom,int u,int v,int w,short int fx2,short int fy2,short int fz2,float *****XXY,Node ***nodeS,Node ***nodeT){
	/*
	Funcion para contar los triangulos en dos 
	nodos con dos puntos en N1 y un punto en N2.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	
	*/
	int i,j,k;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3,W;
	
	int a_,b_,c_;
	int t,p,q, t_,p_,q_;
	float dx,dy,dz1,dz2,dz3;
	float cth1,cth2,cth3, cth1_,cth2_,cth3_;

	for (i=0; i<nodeS[row][col][mom].len; ++i){
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeS[u][v][w].len; ++j){
		x2 = nodeT[u][v][w].elements[j].x - fx2*size_box;
		y2 = nodeT[u][v][w].elements[j].y - fy2*size_box;
		z2 = nodeT[u][v][w].elements[j].z - fz2*size_box;
		w2 = nodeT[u][v][w].elements[j].w;
		dx = x2-x1;
		dy = y2-y1;
		dz1 = z2-z1;
		d12 = dx*dx+dy*dy+dz1*dz1;
		if (d12<dd_max){
		d12 = sqrt(d12);
			for (k=i+1; k<nodeS[row][col][mom].len; ++k){
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
				a_ = (int)(d12*ds);
				b_ = (int)(d13*ds);
				c_ = (int)(d23*ds);
				
				t = (int)(cth1*ds_th);
				p = (int)(cth2*ds_th);
				q = (int)(cth3*ds_th);
				
				t_ = (int)(cth1_*ds_th);
				p_ = (int)(cth2_*ds_th);
				q_ = (int)(cth3_*ds_th);
				
				// Guardamos
				W = w1*w2*w3;
				*(*(*(*(*(XXY+a_)+b_)+c_)+t)+p)+=W;
				*(*(*(*(*(XXY+c_)+a_)+b_)+p_)+q_)+=W;
				*(*(*(*(*(XXY+b_)+c_)+a_)+q)+t_)+=W;
				*(*(*(*(*(XXY+b_)+a_)+c_)+p)+t)+=W;
				*(*(*(*(*(XXY+c_)+b_)+a_)+q_)+p_)+=W;
				*(*(*(*(*(XXY+a_)+c_)+b_)+t_)+q)+=W;
			}
			}
			}
		}
		}
	}
}
//=================================================================== 
void NODE3P::front_node_123_xxy(int row, int col, int mom, int u, int v, int w, int a, int b, int c, bool conx, bool cony, bool conz, float *****XXY, Node ***nodeS, Node ***nodeT){
	
	float fx1N = nodeS[row][0][0].nodepos.fx, x1N = nodeS[row][0][0].nodepos.x;
	float fy1N = nodeS[row][col][0].nodepos.fy, y1N = nodeS[row][col][0].nodepos.y;
	float fz1N = nodeS[row][col][mom].nodepos.fz, z1N = nodeS[row][col][mom].nodepos.z;
	float fx2N = nodeT[u][0][0].nodepos.fx, x2N = nodeT[u][0][0].nodepos.x;
	float fy2N = nodeT[u][v][0].nodepos.fy, y2N = nodeT[u][v][0].nodepos.y;
	float fz2N = nodeT[u][v][w].nodepos.fz, z2N = nodeT[u][v][w].nodepos.z;
	float fx3N = nodeS[a][0][0].nodepos.fx, x3N = nodeS[a][0][0].nodepos.x;
	float fy3N = nodeS[a][b][0].nodepos.fy, y3N = nodeS[a][b][0].nodepos.y;
	float fz3N = nodeS[a][b][c].nodepos.fz, z3N = nodeS[a][b][c].nodepos.z;
	
	float dx,dy,dz;
	
	if (conx){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = y3N-y1N;
			dy *= dy; 
			dz = z3N-z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = y3N-y2N;
				dy *= dy; 
				dz = z3N-z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123_xxy(row,col,mom,u,v,w,a,b,c,fx1N,0,0,fx2N,0,0,fx3N,0,0,XXY,nodeS,nodeT);
				}
			}
		}
	
	}
	if (cony){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = x3N-x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = z3N-z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = x3N-x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = z3N-z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123_xxy(row,col,mom,u,v,w,a,b,c,0,fy1N,0,0,fy2N,0,0,fy3N,0,XXY,nodeS,nodeT);
				}
			}
		}
	}
	if (conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz =(z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = x3N-x1N;
			dx *= dx;
			dy = y3N-y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = x3N-x2N;
				dx *= dx;
				dy = y3N-y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123_xxy(row,col,mom,u,v,w,a,b,c,0,0,fz1N,0,0,fz2N,0,0,fz3N,XXY,nodeS,nodeT);
				}
			}
		}
	}
	if (conx && cony){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = z3N-z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = z3N-z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123_xxy(row,col,mom,u,v,w,a,b,c,fx1N,fy1N,0,fx2N,fy2N,0,fx3N,fy3N,0,XXY,nodeS,nodeT);
				}
			}
		}
	}
	if (conx && conz){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = (z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = y3N-y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = y3N-y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123_xxy(row,col,mom,u,v,w,a,b,c,fx1N,0,fz1N,fx2N,0,fz2N,fx3N,0,fz3N,XXY,nodeS,nodeT);
				}
			}
		}
	}
	if (cony && conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = (z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = x3N-x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = x3N-x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123_xxy(row,col,mom,u,v,w,a,b,c,0,fy1N,fz1N,0,fy2N,fz2N,0,fy3N,fz3N,XXY,nodeS,nodeT);
				}
			}
		}
	}
	if (conx && cony && conz){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = (z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123_xxy(row,col,mom,u,v,w,a,b,c,fx1N,fy1N,fz1N,fx2N,fy2N,fz2N,fx3N,fy3N,fz3N,XXY,nodeS,nodeT);
				}
			}
		}
	}	
}
//=================================================================== 
void NODE3P::front_123_xxy(int row,int col,int mom,int u,int v,int w,int a,int b,int c,short int fx1N,short int fy1N,short int fz1N,short int fx2N,short int fy2N,short int fz2N,short int fx3N,short int fy3N,short int fz3N,float *****XXY, Node ***nodeS, Node ***nodeT){
	/*
	Funcion para contar los triangulos en tres 
	nodos con un puntos en N1, un punto en N2
	y un punto en N3.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	a, b, c => posición de N3.
	
	*/
	
	int i,j,k;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3,W;
	
	int a_,b_,c_;
	int t,p,q, t_,p_,q_;
	float dx,dy,dz1,dz2,dz3;
	float cth1,cth2,cth3, cth1_,cth2_,cth3_;
	
	for (i=0; i<nodeS[row][col][mom].len; ++i){
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeT[u][v][w].len; ++j){
		x2 = nodeT[u][v][w].elements[j].x;
		y2 = nodeT[u][v][w].elements[j].y;
		z2 = nodeT[u][v][w].elements[j].z;
		w2 = nodeT[u][v][w].elements[j].w;
		dx = (x2-((fx2N-fx1N)*size_box_2))-x1;
		dy = (y2-((fy2N-fy1N)*size_box_2))-y1;
		dz1 = (z2-((fz2N-fz1N)*size_box_2))-z1;
		d12 = dx*dx+dy*dy+dz1*dz1;
		if (d12<dd_max){
		d12 = sqrt(d12);
			for (k=0; k<nodeS[a][b][c].len; ++k){
			x3 = nodeS[a][b][c].elements[k].x;
			y3 = nodeS[a][b][c].elements[k].y;
			z3 = nodeS[a][b][c].elements[k].z;
			w3 = nodeS[a][b][c].elements[k].w;
			dx = (x3-((fx3N-fx1N)*size_box_2))-x1;
			dy = (y3-((fy3N-fy1N)*size_box_2))-y1;
			dz2 = (z3-((fz3N-fz1N)*size_box_2))-z1;
			d13 = dx*dx+dy*dy+dz2*dz2;
			if (d13<dd_max){
			d13 = sqrt(d13);
			dx = (x3-((fx3N-fx2N)*size_box_2))-x2;
			dy = (y3-((fy3N-fy2N)*size_box_2))-y2;
			dz3 = (z3-((fz3N-fz2N)*size_box_2))-z2;
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
				*(*(*(*(*(XXY+a_)+b_)+c_)+t)+p)+=W;
				*(*(*(*(*(XXY+c_)+a_)+b_)+p_)+q_)+=W;
				*(*(*(*(*(XXY+b_)+c_)+a_)+q)+t_)+=W;
				*(*(*(*(*(XXY+b_)+a_)+c_)+p)+t)+=W;
				*(*(*(*(*(XXY+c_)+b_)+a_)+q_)+p_)+=W;
				*(*(*(*(*(XXY+a_)+c_)+b_)+t_)+q)+=W;
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


