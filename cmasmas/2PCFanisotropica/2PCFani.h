#include <iomanip> 
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

class NODE2P{
	//Atributos de clase:
	private:
		// Asignados
		int bn;
		int n_pts;
		float size_box;
		float size_node;
		float d_max;
		Node ***nodeD;
		PointW3D *dataD;
		Node ***nodeR;
		PointW3D *dataR;
		// Derivados
		float dd_max;
		float corr;
		float front;
		float ds;
		float ddmax_nod;
		
	private: 
		void make_nodos(Node ***, PointW3D *);
		void add(PointW3D *&, int&, float, float, float, float);
	
	// Métodos de Clase:
	public:
		//Constructor de clase:
		NODE2P(int _bn, int _n_pts, float _size_box, float _size_node, float _d_max, PointW3D *_dataD, Node ***_nodeD, PointW3D *_dataR, Node ***_nodeR){
			
			// Asignados
			bn = _bn;
			n_pts = _n_pts;
			size_box = _size_box;
			size_node = _size_node;
			d_max = _d_max;
			dataD = _dataD;
			nodeD = _nodeD;
			dataR = _dataR;
			nodeR = _nodeR;
			
			// Derivados
			dd_max = d_max*d_max;
			front = size_box - d_max;
			corr = size_node*sqrt(3);
			ds = ((float)(bn))/d_max;
			ddmax_nod = d_max+corr;
			ddmax_nod *= ddmax_nod; 
			
			make_nodos(nodeD,dataD); 
			make_nodos(nodeR,dataR); 
			std::cout << "Terminé de contruir nodos..." << std::endl;
		}
		
		Node ***meshData(){
			return nodeD;
		};
		Node ***meshRand(){
			return nodeR;
		};
		
		// Implementamos Método de mallas:
		void make_histoXX(float **, Node ***);
		void make_histoXY(float **, Node ***, Node ***);
		void histo_front(float **, Node ***, Node ***, float, float, float, bool, bool, bool, int, int, int, int, int, int, float);
		~NODE2P();
};

//=================================================================== 
//==================== Funciones ==================================== 
//===================================================================  

void NODE2P::make_nodos(Node ***nod, PointW3D *dat){
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
		nod[row][col][mom].elements = new PointW3D[0];
	}
	}
	}
	// Llenamos los nodos con los puntos de dat:
	for (i=0; i<n_pts; ++i){
		row = (int)(dat[i].x/size_node);
        	col = (int)(dat[i].y/size_node);
        	mom = (int)(dat[i].z/size_node);
		add( nod[row][col][mom].elements, nod[row][col][mom].len, dat[i].x, dat[i].y, dat[i].z, dat[i].w);
	}
}

//=================================================================== 

void NODE2P::add(PointW3D *&array, int &lon, float _x, float _y, float _z, float _w){
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

void NODE2P::make_histoXX(float **XX, Node ***nodeX){
	/*
	Función para crear los histogramas DD y RR.
	
	Argumentos
	DD: arreglo donde se creará el histograma DD.
	RR: arreglo donde se creará el histograma RR.
	
	*/
	
	int partitions = (int)((size_box/size_node)+1);
	int i, j, row, col, mom, u, v, w;
	float dis, dis_nod;
	float x1D, y1D, z1D, x2D, y2D, z2D;
	float x, y, z, w1;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod;
	float r_ort,r_ort_nod;
	bool con_x, con_y, con_z;
	float d_max_pm = d_max + size_node/2, front_pm = front - size_node/2;
	
	std::cout << "-> Estoy haciendo histograma XX..." << std::endl;
	
	for (row = 0; row < partitions; ++row){
	x1D = nodeX[row][0][0].nodepos.x;
	for (col = 0; col < partitions; ++col){
	y1D = nodeX[row][col][0].nodepos.y;
	for (mom = 0; mom < partitions; ++mom){
	z1D = nodeX[row][col][mom].nodepos.z;			
		//==================================================
		// Distancias entre puntos del mismo nodo:
		//==================================================
		for (i=0; i<nodeX[row][col][mom].len-1; ++i){
		x = nodeX[row][col][mom].elements[i].x;
		y = nodeX[row][col][mom].elements[i].y;
		z = nodeX[row][col][mom].elements[i].z;
		w1 = nodeX[row][col][mom].elements[i].w;
			for (j=i+1; j<nodeX[row][col][mom].len; ++j){
			dx = x-nodeX[row][col][mom].elements[j].x;
			dy = y-nodeX[row][col][mom].elements[j].y;
			dz = z-nodeX[row][col][mom].elements[j].z;
			dz *= dz;
			r_ort = dx*dx+dy*dy;
			if (dz < dd_max && r_ort < dd_max){
			dz = int(sqrt(dz)*ds);
			r_ort = int(sqrt(r_ort)*ds);
			*(*(XX+int(dz))+int(r_ort)) += 2*w1*nodeX[row][col][mom].elements[j].w;
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
		dz_nod *= dz_nod;
		if (dz_nod <= ddmax_nod){
			for (i=0; i<nodeX[row][col][mom].len; ++i){
			x = nodeX[row][col][mom].elements[i].x;
			y = nodeX[row][col][mom].elements[i].y;
			z = nodeX[row][col][mom].elements[i].z;
			w1 = nodeX[row][col][mom].elements[i].w;
				for (j=0; j<nodeX[u][v][w].len; ++j){
				dx = x-nodeX[u][v][w].elements[j].x;
				dy = y-nodeX[u][v][w].elements[j].y;
				dz = z-nodeX[u][v][w].elements[j].z;
				dz *= dz;
				r_ort = dx*dx+dy*dy;
				if (dz < dd_max && r_ort < dd_max){
				dz = int(sqrt(dz)*ds);
				r_ort = int(sqrt(r_ort)*ds);
				*(*(XX+int(dz))+int(r_ort)) += 2*w1*nodeX[u][v][w].elements[j].w;
				}
				}
			}
		}
		//=======================================
		// Distacia de los puntos frontera XX
		//=======================================
		//Condiciones de nodos en frontera:
		con_z = ((z1D<=d_max_pm)&&(z2D>=front_pm))||((z2D<=d_max_pm)&&(z1D>=front_pm));
		if(con_z){
		histo_front(XX,nodeX,nodeX,0.0,0.0,sqrt(dz_nod),false,false,con_z,row,col,mom,u,v,w,2);
		}
		}
		u = row;
		//=========================
		// N2 movil en ZY
		//=========================
		for (v=col+1; v<partitions; ++v){
		y2D = nodeX[u][v][0].nodepos.y;
		dy_nod = y1D-y2D;
		dy_nod *= dy_nod;
			for (w = 0; w<partitions; ++w){		
			z2D = nodeX[u][v][w].nodepos.z;
			dz_nod = z1D-z2D;
			dz_nod *= dz_nod;
			r_ort_nod = dy_nod;
			if (dz_nod <= ddmax_nod && r_ort_nod <= ddmax_nod){
				for (i=0; i<nodeX[row][col][mom].len; ++i){
				x = nodeX[row][col][mom].elements[i].x;
				y = nodeX[row][col][mom].elements[i].y;
				z = nodeX[row][col][mom].elements[i].z;
				w1 = nodeX[row][col][mom].elements[i].w;
					for (j=0; j<nodeX[u][v][w].len; ++j){	
					dx =  x-nodeX[u][v][w].elements[j].x;
					dy =  y-nodeX[u][v][w].elements[j].y;
					dz =  z-nodeX[u][v][w].elements[j].z;
					dz *= dz;
					r_ort = dx*dx+dy*dy;
					if (dz < dd_max && r_ort < dd_max){
					dz = int(sqrt(dz)*ds);
					r_ort = int(sqrt(r_ort)*ds);
					*(*(XX+int(dz))+int(r_ort)) += 2*w1*nodeX[u][v][w].elements[j].w;
					}
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
			histo_front(XX,nodeX,nodeX,0.0,sqrt(dy_nod),sqrt(dz_nod),false,con_y,con_z,row,col,mom,u,v,w,2);
			}
			}
		}
		//=========================
		// N2 movil en ZYX
		//=========================
		for (u=row+1; u<partitions; ++u){
		x2D = nodeX[u][0][0].nodepos.x;
		dx_nod = x1D-x2D;
		dx_nod *= dx_nod;
			for (v=0; v<partitions; ++v){
			y2D = nodeX[u][v][0].nodepos.y;
			dy_nod = y1D-y2D;
			dy_nod *= dy_nod;
				for (w=0; w<partitions; ++w){
				z2D = nodeX[u][v][w].nodepos.z;
				dz_nod = z1D-z2D;
				dz_nod *= dz_nod; 
				r_ort_nod = dx_nod + dy_nod;
				if (dz_nod <= ddmax_nod && r_ort_nod <= ddmax_nod){
					for (i=0; i<nodeX[row][col][mom].len; ++i){
					x = nodeX[row][col][mom].elements[i].x;
					y = nodeX[row][col][mom].elements[i].y;
					z = nodeX[row][col][mom].elements[i].z;
					w1 = nodeX[row][col][mom].elements[i].w;
						for ( j = 0; j < nodeX[u][v][w].len; ++j){	
						dx = x-nodeX[u][v][w].elements[j].x;
						dy = y-nodeX[u][v][w].elements[j].y;
						dz = z-nodeX[u][v][w].elements[j].z;
						dz *= dz;
						r_ort = dx*dx+dy*dy;
						if (dz < dd_max && r_ort < dd_max){
						dz = int(sqrt(dz)*ds);
						r_ort = int(sqrt(r_ort)*ds);
						*(*(XX+int(dz))+int(r_ort)) += 2*w1*nodeX[u][v][w].elements[j].w;
						}
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
				histo_front(XX,nodeX,nodeX,sqrt(dx_nod),sqrt(dy_nod),sqrt(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w,2);
				}	
				}	
			}
		}
	}
	}
	}
}
//=================================================================== 

void NODE2P::make_histoXY(float **XY, Node ***nodeX, Node ***nodeY){
	/*
	Función para crear los histogramas DR.
	
	Argumentos
	XY: arreglo donde se creará el histograma DR.
	
	*/
	
	int partitions = (int)((size_box/size_node)+1);
	int i, j, row, col, mom, u, v, w;
	float dis, dis_nod;
	float x1D, y1D, z1D, x2R, y2R, z2R;
	float x, y, z, w1;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod;
	float r_ort,r_ort_nod;
	bool con_x, con_y, con_z;
	float d_max_pm = d_max + size_node/2, front_pm = front - size_node/2;
	
	std::cout << "-> Estoy haciendo histograma XY..." << std::endl;
	
	for (row = 0; row < partitions; ++row){
	x1D = nodeX[row][0][0].nodepos.x;
	for (col = 0; col < partitions; ++col){
	y1D = nodeX[row][col][0].nodepos.y;
	for (mom = 0; mom < partitions; ++mom){
	z1D = nodeX[row][col][mom].nodepos.z;			
	//=========================
	// N2 movil en ZYX
	//=========================
	for (u=0; u<partitions; ++u){
	x2R = nodeY[u][0][0].nodepos.x;
	dx_nod = x1D-x2R;
	dx_nod *= dx_nod;
		for (v=0; v<partitions; ++v){
		y2R = nodeY[u][v][0].nodepos.y;
		dy_nod = y1D-y2R;
		dy_nod *= dy_nod;
			for (w=0; w<partitions; ++w){
			z2R = nodeY[u][v][w].nodepos.z;
			dz_nod = z1D-z2R;
			dz_nod *= dz_nod; 
			r_ort_nod = dx_nod + dy_nod;
			if (dz_nod <= ddmax_nod && r_ort_nod <= ddmax_nod){
				for (i=0; i<nodeX[row][col][mom].len; ++i){
				x = nodeX[row][col][mom].elements[i].x;
				y = nodeX[row][col][mom].elements[i].y;
				z = nodeX[row][col][mom].elements[i].z;
				w1 = nodeX[row][col][mom].elements[i].w;
					for (j=0; j<nodeY[u][v][w].len; ++j){	
					dx = x-nodeY[u][v][w].elements[j].x;
					dy = y-nodeY[u][v][w].elements[j].y;
					dz = z-nodeY[u][v][w].elements[j].z;
					dz *= dz;
					r_ort = dx*dx+dy*dy;
					if (dz < dd_max && r_ort < dd_max){
					dz = int(sqrt(dz)*ds);
					r_ort = int(sqrt(r_ort)*ds);
					*(*(XY+int(dz))+int(r_ort)) += w1*nodeY[u][v][w].elements[j].w;
					}
					}
				}
			}
			//=======================================
			// Distacia de los puntos frontera
			//=======================================
			//Condiciones de nodos en frontera:
			con_x = ((x1D<=d_max_pm)&&(x2R>=front_pm))||((x2R<=d_max_pm)&&(x1D>=front_pm));
			con_y = ((y1D<=d_max_pm)&&(y2R>=front_pm))||((y2R<=d_max_pm)&&(y1D>=front_pm));
			con_z = ((z1D<=d_max_pm)&&(z2R>=front_pm))||((z2R<=d_max_pm)&&(z1D>=front_pm));
			if(con_x || con_y || con_z){
			histo_front(XY,nodeX,nodeY,sqrt(dx_nod),sqrt(dy_nod),sqrt(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w,1);
			}	
			}	
		}
	}
	}
	}
	}
}
//=================================================================== 

void NODE2P::histo_front(float **PP, Node ***dat, Node ***ran, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int row, int col, int mom, int u, int v, int w, float q){
	int i, j;
	float dis_f,dis,d_x,d_y,d_z;
	float dx_nod, dy_nod, dz_nod;
	float x,y,z,w1;
	float r_ort,r_ort_nod;
	//======================================================================
	// Si los puentos estás en las paredes laterales de X
	if( con_in_x ){
	dx_nod = size_box-(dn_x+size_node);
	dy_nod = dn_y;
	dz_nod = dn_z;
	dz_nod *= dz_nod;
	r_ort_nod = dx_nod*dx_nod + dy_nod*dy_nod;
	if (dz_nod <= ddmax_nod && r_ort_nod <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = fabs(x-ran[u][v][w].elements[j].x)-size_box;
			d_y = y-ran[u][v][w].elements[j].y;
			d_z = z-ran[u][v][w].elements[j].z;
			d_z *= d_z;
			r_ort = d_x*d_x+d_y*d_y;
			if (d_z < dd_max && r_ort < dd_max){
			d_z = int(sqrt(d_z)*ds);
			r_ort = int(sqrt(r_ort)*ds);
			*(*(PP+int(d_z))+int(r_ort)) += q*w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las paredes laterales de Y		
	if( con_in_y ){
	dx_nod = dn_x;
	dy_nod = size_box-(dn_y+size_node);
	dz_nod = dn_z;
	dz_nod *= dz_nod;
	r_ort_nod = dx_nod*dx_nod + dy_nod*dy_nod;
	if (dz_nod <= ddmax_nod && r_ort_nod <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = x-ran[u][v][w].elements[j].x;
			d_y = fabs(y-ran[u][v][w].elements[j].y)-size_box;
			d_z = z-ran[u][v][w].elements[j].z;
			d_z *= d_z;
			r_ort = d_x*d_x+d_y*d_y;
			if (d_z < dd_max && r_ort < dd_max){
			d_z = int(sqrt(d_z)*ds);
			r_ort = int(sqrt(r_ort)*ds);
			*(*(PP+int(d_z))+int(r_ort)) += q*w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las paredes laterales de Z
	if( con_in_z ){
	dx_nod = dn_x;
	dy_nod = dn_y;
	dz_nod = size_box-(dn_z+size_node);
	dz_nod *= dz_nod;
	r_ort_nod = dx_nod*dx_nod + dy_nod*dy_nod;
	if (dz_nod <= ddmax_nod && r_ort_nod <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = x-ran[u][v][w].elements[j].x;
			d_y = y-ran[u][v][w].elements[j].y;
			d_z = fabs(z-ran[u][v][w].elements[j].z)-size_box;
			d_z *= d_z;
			r_ort = d_x*d_x+d_y*d_y;
			if (d_z < dd_max && r_ort < dd_max){
			d_z = int(sqrt(d_z)*ds);
			r_ort = int(sqrt(r_ort)*ds);
			*(*(PP+int(d_z))+int(r_ort)) += q*w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X y Y			
	if( con_in_x && con_in_y ){
	dx_nod = size_box-(dn_x+size_node);
	dy_nod = size_box-(dn_y+size_node);
	dz_nod = dn_z;
	dz_nod *= dz_nod;
	r_ort_nod = dx_nod*dx_nod + dy_nod*dy_nod;
	if (dz_nod <= ddmax_nod && r_ort_nod <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = fabs(x-ran[u][v][w].elements[j].x)-size_box;
			d_y = fabs(y-ran[u][v][w].elements[j].y)-size_box;
			d_z = z-ran[u][v][w].elements[j].z;
			d_z *= d_z;
			r_ort = d_x*d_x+d_y*d_y;
			if (d_z < dd_max && r_ort < dd_max){
			d_z = int(sqrt(d_z)*ds);
			r_ort = int(sqrt(r_ort)*ds);
			*(*(PP+int(d_z))+int(r_ort)) += q*w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X y Z				
	if( con_in_x && con_in_z ){
	dx_nod = size_box-(dn_x+size_node);
	dy_nod = dn_y;
	dz_nod = size_box-(dn_z+size_node);
	dz_nod *= dz_nod;
	r_ort_nod = dx_nod*dx_nod + dy_nod*dy_nod;
	if (dz_nod <= ddmax_nod && r_ort_nod <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = fabs(x-ran[u][v][w].elements[j].x)-size_box;
			d_y = y-ran[u][v][w].elements[j].y;
			d_z = fabs(z-ran[u][v][w].elements[j].z)-size_box;
			d_z *= d_z;
			r_ort = d_x*d_x+d_y*d_y;
			if (d_z < dd_max && r_ort < dd_max){
			d_z = int(sqrt(d_z)*ds);
			r_ort = int(sqrt(r_ort)*ds);
			*(*(PP+int(d_z))+int(r_ort)) += q*w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de Y y Z			
	if( con_in_y && con_in_z ){
	dx_nod = dn_x;
	dy_nod = size_box-(dn_y+size_node);
	dz_nod = size_box-(dn_z+size_node);
	dz_nod *= dz_nod;
	r_ort_nod = dx_nod*dx_nod + dy_nod*dy_nod;
	if (dz_nod <= ddmax_nod && r_ort_nod <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = x-ran[u][v][w].elements[j].x;
			d_y = fabs(y-ran[u][v][w].elements[j].y)-size_box;
			d_z = fabs(z-ran[u][v][w].elements[j].z)-size_box;
			d_z *= d_z;
			r_ort = d_x*d_x+d_y*d_y;
			if (d_z < dd_max && r_ort < dd_max){
			d_z = int(sqrt(d_z)*ds);
			r_ort = int(sqrt(r_ort)*ds);
			*(*(PP+int(d_z))+int(r_ort)) += q*w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================
	// Si los puentos estás en las esquinas que cruzan las paredes laterales de X, Y y Z		
	if( con_in_x && con_in_y && con_in_z ){
	dx_nod = size_box-(dn_x+size_node);
	dy_nod = size_box-(dn_y+size_node);
	dz_nod = size_box-(dn_z+size_node);
	dz_nod *= dz_nod;
	r_ort_nod = dx_nod*dx_nod + dy_nod*dy_nod;
	if (dz_nod <= ddmax_nod && r_ort_nod <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = fabs(x-ran[u][v][w].elements[j].x)-size_box;
			d_y = fabs(y-ran[u][v][w].elements[j].y)-size_box;
			d_z = fabs(z-ran[u][v][w].elements[j].z)-size_box;
			d_z *= d_z;
			r_ort = d_x*d_x+d_y*d_y;
			if (d_z < dd_max && r_ort < dd_max){
			d_z = int(sqrt(d_z)*ds);
			r_ort = int(sqrt(r_ort)*ds);
			*(*(PP+int(d_z))+int(r_ort)) += q*w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
}

NODE2P::~NODE2P(){
	
}
