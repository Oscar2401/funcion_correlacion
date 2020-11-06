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
		// Derivados
		float ll;
		float dd_max;
		float corr;
		float front;
		float ds;
		float ddmax_nod;
		float d_max_pm;
		float front_pm;
		
	private: 
		void make_nodos(Node ***, PointW3D *);
		void add(PointW3D *&, int&, float, float, float, float, short int, short int, short int);
	
	// Métodos de Clase:
	public:
		//Constructor de clase:
		NODE3P(int _bn, int _n_pts, float _size_box, float _size_node, float _d_max, PointW3D *_dataD, Node ***_nodeD){
			bn = _bn;
			n_pts = _n_pts;
			size_box = _size_box;
			size_node = _size_node;
			d_max = _d_max;
			dataD = _dataD;
			nodeD = _nodeD;
			//======================
			dd_max = d_max*d_max;
			front = size_box - d_max;
			corr = size_node*sqrt(3);
			ds = ((float)(bn))/d_max;
			ddmax_nod = (d_max+corr)*(d_max+corr);
			partitions = (int)(ceil(size_box/size_node));
			d_max_pm = d_max + size_node/2;
			front_pm = front - size_node/2;
			
			make_nodos(nodeD,dataD);
			std::cout << "Terminé de construir nodos..." << std::endl;
		}
		
		Node ***meshData(){
			return nodeD;
		};
		
		// Implementamos Método de mallas:
		void make_histoXXX(float ***, Node ***);
		void count_3_N111(int, int, int, float ***, Node ***);
		void count_3_N112(int, int, int, int, int, int, float ***, Node ***);
		void count_3_N123(int, int, int, int, int, int, int, int, int, float ***, Node ***);
		void symmetrize(float ***);
		void symmetrize_analitic(float ***);
		void BPC(float ***, PointW3D *);
		
		void front_node_112(int, int, int, int, int, int, bool, bool, bool, float ***, Node ***);
		void front_112(int,int,int,int,int,int,short int,short int,short int,float ***,Node ***);
		void front_node_123(int,int,int,int,int,int,int,int,int,bool,bool,bool,float ***,Node ***);
		void front_123(int,int,int,int,int,int,int,int,int,short int,short int,short int,short int,short int,short int,short int,short int,short int,float ***, Node ***);
		
		void make_histo_analitic(float ***, float ***, Node ***);
		void make_histoXX(float *, float *, Node ***, int);
		void histo_front_XX(float *, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int, float);
		
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
void NODE3P::make_histoXXX(float ***XXX, Node ***nodeX){
	/*
	Función para crear los histogramas DDD.
	
	Argumentos
	XXX_: arreglo donde se creará el histograma DDD.
	nodeX: malla de datos
	
	*/ 
	int i, j, k, row, col, mom, u, v, w, a ,b, c;
	float dis, dis_nod, dis_nod2, dis_nod3;
	float x1N, y1N, z1N, x2N, y2N, z2N, x3N, y3N, z3N;
	float x, y, z;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod, dx_nod2, dy_nod2, dz_nod2, dx_nod3, dy_nod3, dz_nod3;

	
	std::cout << "-> Estoy haciendo histograma DDD..." << std::endl;
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
	dz_nod *= dz_nod;
	
	// Caja
	
	if (dz_nod <= ddmax_nod){
	//==============================================
	// 2 puntos en N1 y 1 punto en N2
	//==============================================
	count_3_N112(row, col, mom, u, v, w, XXX, nodeX);
	//==============================================
	// 1 punto en N1, 1 punto en N2 y 1 punto en N3
	//==============================================
	a = u;
	b = v;
		x3N = nodeX[a][0][0].nodepos.x;
		y3N = nodeX[a][b][0].nodepos.y;
		//=======================
		// Nodo 3 movil en Z:
		//=======================
		for (c=w+1; c<partitions; ++c){
		z3N = nodeX[a][b][c].nodepos.z; 
		dz_nod2 = z3N-z1N;
		dz_nod2 *= dz_nod2;
		if (dz_nod2 <= ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
		}
		//=======================
		// Nodo 3 movil en ZY:
		//=======================
		for (b=v+1; b<partitions; ++b){
		y3N = nodeX[a][b][0].nodepos.y;
		dy_nod2 = y3N-y1N;
		dy_nod2 *= dy_nod2;
			for (c=0; c<partitions; ++c){
			z3N = nodeX[a][b][c].nodepos.z;
			dz_nod2 = z3N-z1N;
			dz_nod2 *= dz_nod2;
			dis_nod2 = dy_nod2 + dz_nod2;
			if (dis_nod2 <= ddmax_nod){
			dy_nod3 = y3N-y2N;
			dz_nod3 = z3N-z2N;
			dis_nod3 = dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
			if (dis_nod3 <= ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
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
				for (c=0; c<partitions; ++c){
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod2 = z3N-z1N;
				dz_nod2 *= dz_nod2;
				dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
				if (dis_nod2 <= ddmax_nod){
				dx_nod3 = x3N-x2N;
				dy_nod3 = y3N-y2N;
				dz_nod3 = z3N-z2N;
				dis_nod3 = dx_nod3*dx_nod3 + dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
				if (dis_nod3 <= ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
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
			if (dis_nod3 <= ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
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
				dz_nod3 = z3N-z2N;
				dis_nod3 = dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
				if (dis_nod3 <= ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
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
					dy_nod3 = y3N-y2N;
					dz_nod3 = z3N-z2N;
					dis_nod3 = dx_nod3*dx_nod3 + dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
					if (dis_nod3 <= ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
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
				if (dis_nod3 <= ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
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
					dz_nod3 = z3N-z2N;
					dis_nod3 = dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
					if (dis_nod3 <= ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
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
						dy_nod3 = y3N-y2N;
						dz_nod3 = z3N-z2N;
						dis_nod3 = dx_nod3*dx_nod3 + dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
						if (dis_nod3 <= ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, XXX, nodeX);
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
}
//=================================================================== 
void NODE3P::make_front(Node ***nod, PointW3D *dat){
	/*
	Función para crear los nodos con los datos y puntos random
	
	Argumentos
	nod: arreglo donde se crean los nodos.
	dat: datos a dividir en nodos.
	
	*/
	int i, row, col, mom;
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

void NODE3P::BPC(){

	

}

//=================================================================== 
void NODE3P::count_3_N111(int row, int col, int mom, float ***XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en un mismo Nodo.
	
	row, col, mom => posición del Nodo. Esto define al Nodo.
	
	*/
	int i,j,k;
	float dx,dy,dz;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3;

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
		dz = z2-z1;
		d12 = dx*dx+dy*dy+dz*dz;
		if (d12<=dd_max){
		d12 = sqrt(d12);
			for (k=j+1; k<nodeS[row][col][mom].len; ++k){ 
			x3 = nodeS[row][col][mom].elements[k].x;
			y3 = nodeS[row][col][mom].elements[k].y;
			z3 = nodeS[row][col][mom].elements[k].z;
			w3 = nodeS[row][col][mom].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz = z3-z1;
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<=dd_max){
			d13 = sqrt(d13);
			dx = x3-x2;
			dy = y3-y2;
			dz = z3-z2;
			d23 = dx*dx+dy*dy+dz*dz;
			if (d23<=dd_max){
				d23 = sqrt(d23);
				*(*(*(XXX+(int)(d12*ds))+(int)(d13*ds))+(int)(d23*ds))+=w1*w2*w3;
			}
			}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N112(int row, int col, int mom, int u, int v, int w, float ***XXX, Node ***nodeS){
	/*
	Funcion para contar los triangulos en dos 
	nodos con dos puntos en N1 y un punto en N2.
	
	row, col, mom => posición de N1.
	u, v, w => posición de N2.
	
	*/
	int i,j,k;
	float dx,dy,dz;
	float d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3;

	for (i=0; i<nodeS[u][v][w].len; ++i){
	x1 = nodeS[u][v][w].elements[i].x;
	y1 = nodeS[u][v][w].elements[i].y;
	z1 = nodeS[u][v][w].elements[i].z;
	w1 = nodeS[u][v][w].elements[i].w;
		for (j=0; j<nodeS[row][col][mom].len; ++j){
		x2 = nodeS[row][col][mom].elements[j].x;
		y2 = nodeS[row][col][mom].elements[j].y;
		z2 = nodeS[row][col][mom].elements[j].z;
		w2 = nodeS[row][col][mom].elements[j].w;
		dx = x2-x1;
		dy = y2-y1;
		dz = z2-z1;
		d12 = dx*dx+dy*dy+dz*dz;
		if (d12<=dd_max){
		d12 = sqrt(d12);
			for (k=j+1; k<nodeS[row][col][mom].len; ++k){
			x3 = nodeS[row][col][mom].elements[k].x;
			y3 = nodeS[row][col][mom].elements[k].y;
			z3 = nodeS[row][col][mom].elements[k].z;
			w3 = nodeS[row][col][mom].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz = z3-z1;
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<=dd_max){
			d13 = sqrt(d13);
			dx = x3-x2;
			dy = y3-y2;
			dz = z3-z2;
			d23 = dx*dx+dy*dy+dz*dz;
			if (d23<=dd_max){
				d23 = sqrt(d23);
				*(*(*(XXX+(int)(d12*ds))+(int)(d13*ds))+(int)(d23*ds))+=w1*w2*w3;
			}
			}
			}
			for (k=i+1; k<nodeS[u][v][w].len; ++k){
			x3 = nodeS[u][v][w].elements[k].x;
			y3 = nodeS[u][v][w].elements[k].y;
			z3 = nodeS[u][v][w].elements[k].z;
			w3 = nodeS[u][v][w].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz = z3-z1;
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<=dd_max){
			d13 = sqrt(d13);
			dx = x3-x2;
			dy = y3-y2;
			dz = z3-z2;
			d23 = dx*dx+dy*dy+dz*dz;
			if (d23<=dd_max){
				d23 = sqrt(d23);
				*(*(*(XXX+(int)(d12*ds))+(int)(d13*ds))+(int)(d23*ds))+=w1*w2*w3;
			}
			}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N123(int row, int col, int mom, int u, int v, int w, int a, int b, int c, float ***XXX, Node ***nodeS){
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
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3;
	for (i=0; i<nodeS[row][col][mom].len; ++i){
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeS[a][b][c].len; ++j){
		x3 = nodeS[a][b][c].elements[j].x;
		y3 = nodeS[a][b][c].elements[j].y;
		z3 = nodeS[a][b][c].elements[j].z;
		w3 = nodeS[a][b][c].elements[j].w;
		dx = x3-x1;
		dy = y3-y1;
		dz = z3-z1;
		d13 = dx*dx+dy*dy+dz*dz;
		if (d13<=dd_max){
			for (k=0; k<nodeS[u][v][w].len; ++k){
			x2 = nodeS[u][v][w].elements[k].x;
			y2 = nodeS[u][v][w].elements[k].y;
			z2 = nodeS[u][v][w].elements[k].z;
			w2 = nodeS[u][v][w].elements[k].w;
			dx = x3-x2;
			dy = y3-y2;
			dz = z3-z2;
			d23 = dx*dx+dy*dy+dz*dz;
			if (d23<=dd_max){
			dx = x2-x1;
			dy = y2-y1;
			dz = z2-z1;
			d12 = dx*dx+dy*dy+dz*dz;
			if (d12<=dd_max) *(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=w1*w2*w3;
			}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::symmetrize(float ***XXX){
	/*
	Función para simetrizar histograma
	
	Argumentos
	XXX_: arreglo a simetrizar
	
	*/ 
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
//================ Funciones para formulas analíticas =============== 
//=================================================================== 

void NODE3P::make_histo_analitic(float ***XXY, float ***XXX, Node ***nodeX){
	/*
	Función para construir histogramas con funciones analíticas
	
	Argumentos
	XXX_: arreglo random
	XXY: arreglo mixto
	nodeX: malla de datos
	
	*/ 
	//======================================
	// Histograma RRR, DDR y DRR (ANALITICA)
	//======================================
	std::cout << "-> Estoy haciendo histograma RRR y DDR..." << std::endl;
	
	int i, j, k, u, v, w, a, b, c;
	short int v_in;
	bool con;
	
	// Constantes para RRR
	float dr = d_max/(double)bn;
	float dr2 = dr/2, r1, r2, r3;
	float ri, rj, rk;
	
	float V = size_box*size_box*size_box;
	float beta = n_pts/V;
	float gama  = 3*8*(2*acos(0.0))*(2*acos(0.0))*(n_pts*beta*beta);
	float alph = gama*dr*dr*dr;
	
	// Refinamiento para RRR
	int bn_ref = 200; 
	float dr_ref = dr/bn_ref;
	float dr_ref2 = dr_ref/2;
	float ru, rv, rw;
	
	float alph_ref = gama*dr_ref*dr_ref*dr_ref;
	float S_av;
	float f_av;
	
	//=============================================================================================
	
	// Constantes para DDR
	int ptt = 100;
	int bins = ptt*bn;
	float dr_ptt = d_max/(float)bins;
	float dr_ptt2 = dr_ptt/2;
	float  *DD;
	float *RR;
	DD = new float[bins];
	RR = new float[bins];
	for (i = 0; i < bins; ++i){
		*(DD+i) = 0; 
		*(RR+i) = 0.0;
	}
	
	// Hacemos los histogramas DD y RR
	make_histoXX(DD, RR, nodeX, bins);
	
	// Iniciamos un arreglo para la funcion f_averrage
	float *ff_av;
	ff_av = new float[bn];
	for (i=0; i<bn; ++i) *(ff_av+i) = 0;
	
	// Calculamos las f_averrage
	int i_;
	for(i=0; i<bn; ++i){
	f_av = 0;
	ri = i*dr;
	i_ = i*ptt;
		for(j = 0; j<ptt; ++j){
		rj = (j+0.5)*dr_ptt;
		f_av += (ri + rj)*(((float)(*(DD+i_+j)) /(*(RR+i_+j))) - 1);
		}
	*(ff_av+i) = f_av/ptt;
	}
	
	delete[] DD;
	delete[] RR;
	
	//=============================================================================================
	
	// Refinamiento para DDR
	int bins_ref = ptt*bn_ref*bn;
	float dr_ptt_ref = d_max/(float)bins_ref;
	float dr_ptt_ref2 = dr_ptt_ref/2;
	DD = new float[bins_ref];
	RR = new float[bins_ref];
	for (i = 0; i < bins_ref; i++){
		*(DD+i) = 0; 
		*(RR+i) = 0.0;
	}
	
	// Hacemos los histogramas DD y RR
	make_histoXX(DD, RR, nodeX, bins_ref);
	
	// Iniciamos un arreglo para la funcion f_averrage del refinamiento
	float *ff_av_ref;
	ff_av_ref = new float[bn_ref*bn];
	for (i=0; i<bn_ref*bn; ++i) *(ff_av_ref+i) = 0.0;
	
	// Calculamos las f_averrage del refinamiento
	int j_;
	for(i=0; i<bn; ++i){
	ri = i*dr;
	i_ = i*bn_ref;
		for(j=0; j<bn_ref; ++j){
		f_av = 0;
		rj = j*dr_ref;
		j_ = j*ptt;
			for(k = 0; k<ptt; ++k){
				rk = (k+0.5)*dr_ptt_ref;
				f_av += (ri+rj+rk)*(((float)(*(DD+i_+j_+k))/(*(RR+i_+j_+k)))-1);
			}
			*(ff_av_ref+i_+j) = f_av/ptt;
		}
	}
	
	delete[] DD;
	delete[] RR;
	
	//=============================================================================================
	
	float c_RRR;
	
	for(i=0; i<bn; ++i) {
	ri = i*dr;
	i_ = i*bn_ref;
	for(j=i; j<bn; ++j) {
	rj = j*dr;
	j_ = j*bn_ref;
	for(k=j; k<bn; ++k) {
	rk = k*dr;
		// Checamos vertices del cubo para hacer 
		// o no refinamiento  
		
		v_in = 0;
		
		for (a = 0; a < 2; ++a){
		r1 = ri + (a*dr);
		for (b = 0; b < 2; ++b){
		r2 = rj + (b*dr);
		for (c = 0; c < 2; ++c){
		r3 = rk + (c*dr);	
			if (r1 + r2 >= r3 && r1 + r3 >= r2 && r2 + r3 >= r1) ++v_in; 
		}
		}
		}
		
		if (v_in==8){
			*(*(*(XXX+i)+j)+k) += alph*(ri+dr2)*(rj+dr2)*(rk+dr2);
			*(*(*(XXY+i)+j)+k) += *(*(*(XXX+i)+j)+k)*(1+(*(ff_av+i)/(3*(ri+dr_ptt2)))+(*(ff_av+j)/(3*(rj+dr_ptt2)))+(*(ff_av+k)/(3*(rk+dr_ptt2))));
		}
		
		else if (v_in < 8 && v_in > 0){
			
			S_av = 0;
			f_av = 0;
			con = false;
			
			for(u=0; u<bn_ref; ++u) {
			ru = ri + (u*dr_ref);
			for(v=0; v<bn_ref; ++v) {
			rv = rj + (v*dr_ref);
			for(w=0; w<bn_ref; ++w) {
			rw = rk + (w*dr_ref);
			
				v_in = 0;
		
				for (a = 0; a < 2; ++a){
				r1 = ru + (a*dr_ref);
				for (b = 0; b < 2; ++b){
				r2 = rv + (b*dr_ref);
				for (c = 0; c < 2; ++c){
				r3 = rw + (c*dr_ref);	
					if (r1 + r2 >= r3 && r1 + r3 >= r2 && r2 + r3 >= r1) v_in++;
				}
				}
				}
				if (v_in==8){
					c_RRR = (ru+dr_ref2)*(rv+dr_ref2)*(rw+dr_ref2);
					S_av += c_RRR;
					f_av += c_RRR*(1+(*(ff_av_ref+i_+u)/(3*(ru+dr_ref2)))+(*(ff_av_ref+j_+v)/(3*(rv+dr_ref2)))+(*(ff_av_ref+(k*bn_ref)+w)/(3*(rw+dr_ref2))));
					con = true;
				}
			}
			}
			}
			if (con){
				*(*(*(XXX+i)+j)+k) += alph_ref*S_av;
				*(*(*(XXY+i)+j)+k) += alph_ref*f_av;
			}
		}
	}
	}
	}
	
	symmetrize_analitic(XXX);
	symmetrize_analitic(XXY); 
	
}
//=================================================================== 
void NODE3P::make_histoXX(float *XX, float *YY, Node ***nodeX, int bins){
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
	float x, y, z, a;
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
		for (i= 0; i<nodeX[row][col][mom].len-1; ++i){
		x = nodeX[row][col][mom].elements[i].x;
		y = nodeX[row][col][mom].elements[i].y;
		z = nodeX[row][col][mom].elements[i].z;
		a = nodeX[row][col][mom].elements[i].w;
			for (j=i+1; j<nodeX[row][col][mom].len; ++j){
			dx = x-nodeX[row][col][mom].elements[j].x;
			dy = y-nodeX[row][col][mom].elements[j].y;
			dz = z-nodeX[row][col][mom].elements[j].z;
			dis = dx*dx+dy*dy+dz*dz;
			if (dis <= dd_max) *(XX + (int)(sqrt(dis)*ds_new)) += a*nodeX[row][col][mom].elements[j].w;
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
		dz_nod *= dz_nod;
		if (dz_nod <= ddmax_nod){
			for (i=0; i<nodeX[row][col][mom].len; ++i){
			x = nodeX[row][col][mom].elements[i].x;
			y = nodeX[row][col][mom].elements[i].y;
			z = nodeX[row][col][mom].elements[i].z;
			a = nodeX[row][col][mom].elements[i].w;
				for (j=0; j<nodeX[u][v][w].len; ++j){
				dx = x-nodeX[u][v][w].elements[j].x;
				dy = y-nodeX[u][v][w].elements[j].y;
				dz = z-nodeX[u][v][w].elements[j].z;
				dis = dx*dx+dy*dy+dz*dz;
				if (dis <= dd_max) *(XX + (int)(sqrt(dis)*ds_new)) +=  a*nodeX[row][col][mom].elements[j].w;
				}
			}
			}
			//=======================================
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
		dy_nod *= dy_nod;
			for (w = 0; w < partitions ; ++w){		
			z2D = nodeX[u][v][w].nodepos.z;
			dz_nod = z1D-z2D;
			dz_nod *= dz_nod;
			dis_nod = dy_nod + dz_nod;
			if (dis_nod <= ddmax_nod){
				for ( i = 0; i < nodeX[row][col][mom].len; ++i){
				x = nodeX[row][col][mom].elements[i].x;
				y = nodeX[row][col][mom].elements[i].y;
				z = nodeX[row][col][mom].elements[i].z;
				a = nodeX[row][col][mom].elements[i].w;
					for ( j = 0; j < nodeX[u][v][w].len; ++j){	
					dx =  x-nodeX[u][v][w].elements[j].x;
					dy =  y-nodeX[u][v][w].elements[j].y;
					dz =  z-nodeX[u][v][w].elements[j].z;
					dis = dx*dx+dy*dy+dz*dz;
					if (dis <= dd_max) *(XX + (int)(sqrt(dis)*ds_new)) += a*nodeX[row][col][mom].elements[j].w;
					}
				}
			}
			//=======================================
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
		dx_nod *= dx_nod;
			for ( v = 0; v < partitions; ++v){
			y2D = nodeX[u][v][0].nodepos.y;
			dy_nod = y1D-y2D;
			dy_nod *= dy_nod;
				for ( w = 0; w < partitions; ++w){
				z2D = nodeX[u][v][w].nodepos.z;
				dz_nod = z1D-z2D;
				dz_nod *= dz_nod;
				dis_nod = dx_nod + dy_nod + dz_nod;
				if (dis_nod <= ddmax_nod){
					for ( i = 0; i < nodeX[row][col][mom].len; ++i){
					x = nodeX[row][col][mom].elements[i].x;
					y = nodeX[row][col][mom].elements[i].y;
					z = nodeX[row][col][mom].elements[i].z;
					a = nodeX[row][col][mom].elements[i].w;
						for ( j = 0; j < nodeX[u][v][w].len; ++j){	
						dx = x-nodeX[u][v][w].elements[j].x;
						dy = y-nodeX[u][v][w].elements[j].y;
						dz = z-nodeX[u][v][w].elements[j].z;
						dis = dx*dx + dy*dy + dz*dz;
						if (dis <= dd_max) *(XX + (int)(sqrt(dis)*ds_new)) += a*nodeX[row][col][mom].elements[j].w;
						}
					}
				}
				//=======================================
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
	double dr = (d_max/(double)bins);
	double V = size_box*size_box*size_box;
	double beta1 = n_pts*n_pts/V;
	double alph = 4*(2*acos(0.0))*(beta1)*dr*dr*dr/3;
	double r1, r2;
	for(int a=0; a<bins; ++a) {
		r2 = (double)(a);
		r1 = r2+1;
        	*(YY+a) += alph*((r1*r1*r1)-(r2*r2*r2));
	}
}

//=================================================================== 
void NODE3P::histo_front_XX(float *PP, Node ***dat, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int row, int col, int mom, int u, int v, int w, float ds_new){
	int i, j;
	float dis_f, dis, d_x, d_y, d_z;
	float x, y, z, a;
	//======================================================================
	// Si los puentos estás en las paredes laterales de X
	if( con_in_x ){
	dis_f = disn + ll - 2*dn_x*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = y-dat[u][v][w].elements[j].y;
			d_z = z-dat[u][v][w].elements[j].z;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis <= dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += a*dat[u][v][w].elements[j].w;
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las paredes laterales de Y		
	if( con_in_y ){
	dis_f = disn + ll - 2*dn_y*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = x-dat[u][v][w].elements[j].x;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = z-dat[u][v][w].elements[j].z;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis <= dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += a*dat[u][v][w].elements[j].w;
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las paredes laterales de Z
	if( con_in_z ){
	dis_f = disn + ll - 2*dn_z*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = x-dat[u][v][w].elements[j].x;
			d_y = y-dat[u][v][w].elements[j].y;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis <= dd_max)*(PP + (int)(sqrt(dis)*ds_new)) += a*dat[u][v][w].elements[j].w;
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X y Y			
	if( con_in_x && con_in_y ){
	dis_f = disn + 2*ll - 2*(dn_x+dn_y)*size_box;
	if (dis_f < ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = z-dat[u][v][w].elements[j].z;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis <= dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += a*dat[u][v][w].elements[j].w;
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X y Z				
	if( con_in_x && con_in_z ){
	dis_f = disn + 2*ll - 2*(dn_x+dn_z)*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = y-dat[u][v][w].elements[j].y;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis <= dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += a*dat[u][v][w].elements[j].w;
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de Y y Z			
	if( con_in_y && con_in_z ){
	dis_f = disn + 2*ll - 2*(dn_y+dn_z)*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		a = dat[row][col][mom].elements[i].w;
			for ( j = 0; j < dat[u][v][w].len; ++j){
			d_x = x-dat[u][v][w].elements[j].x;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis <= dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += a*dat[u][v][w].elements[j].w;
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X, Y y Z		
	if( con_in_x && con_in_y && con_in_z ){
	dis_f = disn + 3*ll - 2*(dn_x+dn_y+dn_z)*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = d_x*d_x + d_y*d_y + d_z*d_z;
			if (dis <= dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += a*dat[u][v][w].elements[j].w;
			}
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
NODE3P::~NODE3P(){	
}


