
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
			std::cout << "Terminé de contruir nodos..." << std::endl;
		}
		
		Node ***meshData(){
			return nodeD;
		};
		Node ***meshRand(){
			return nodeR;
		};
		
		// Implementamos Método de mallas:
		void make_histoXX(unsigned int *, unsigned int *, Node ***);
		void make_histoXY(unsigned int *, Node ***, Node ***);
		void histo_front_XX(unsigned int *, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int);
		void histo_front_XY(unsigned int *, Node ***, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int);
		~NODE();
};

//=================================================================== 
//==================== Funciones ==================================== 
//=================================================================== 
//=================================================================== 

void NODE::make_nodos(Node *** nod, Point3D *dat){
	/*
	Función para crear los nodos con los datos y puntos random
	
	Argumentos
	nod: arreglo donde se crean los nodos.
	dat: datos a dividir en nodos.
	
	*/
	int i, row, col, mom, partitions = (int)((size_box/size_node)+1);
	float p_med = size_node/2;
	
	// Inicializamos los nodos vacíos:
	for ( row = 0; row < partitions; row++){
		for ( col = 0; col < partitions; col++){
			for ( mom = 0; mom < partitions; mom++){
				nod[row][col][mom].nodepos.z = ((float)(mom)*(size_node))+p_med;
				nod[row][col][mom].nodepos.y = ((float)(row)*(size_node))+p_med;
				nod[row][col][mom].nodepos.x = ((float)(col)*(size_node))+p_med;
				nod[row][col][mom].len = 0;
				nod[row][col][mom].elements = new Point3D[0];
			}
		}
	}
	// Llenamos los nodos con los puntos de dat:
	for ( i = 0; i < n_pts; i++){
		col = (int)(dat[i].x/size_node);
        	row = (int)(dat[i].y/size_node);
        	mom = (int)(dat[i].z/size_node);
		add( nod[row][col][mom].elements, nod[row][col][mom].len, dat[i].x, dat[i].y, dat[i].z);
	}
}

//=================================================================== 

void NODE::add(Point3D *&array, int &lon, float _x, float _y, float _z){
	lon++;
	Point3D *array_aux = new Point3D[lon];
	for (int i = 0; i < lon-1; i++){
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

void NODE::make_histoXX(unsigned int *XX, unsigned int *YY, Node ***nodeX){
	/*
	Función para crear los histogramas DD y RR.
	
	Argumentos
	DD: arreglo donde se creará el histograma DD.
	RR: arreglo donde se creará el histograma RR.
	
	*/
	//Variables compartidas en hilos: 
	int partitions = (int)((size_box/size_node)+1);
	std::cout << "-> Estoy haciendo histograma XX..." << std::endl;
	
	#pragma omp parallel num_threads(4) 
    	{
    	unsigned int *SS;
    	SS = new unsigned int[bn];
    	for (int k = 0; k < bn; k++){
		*(SS+k) = 0;
	}
	//Variables privadas en los hilos:
	int i, j, row, col, mom, u, v, w;
	float dis, dis_nod;
	float x1D, y1D, z1D, x2D, y2D, z2D;
	float x, y, z;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod;
	bool con_x, con_y, con_z;
	float d_max_pm = d_max + size_node/2, front_pm = front - size_node/2;
	
    	#pragma omp for collapse(3)  schedule(dynamic)
	for (row = 0; row < partitions; row++){
	for (col = 0; col < partitions; col++){
	for (mom = 0; mom < partitions; mom++){
				
		// Distancias entre puntos del mismo nodo:
		//==================================================
		// Histograma DD
		for ( i= 0; i < nodeX[row][col][mom].len - 1; i++){
		x = nodeX[row][col][mom].elements[i].x;
		y = nodeX[row][col][mom].elements[i].y;
		z = nodeX[row][col][mom].elements[i].z;
			for ( j = i+1; j < nodeX[row][col][mom].len; j++){
			dx = x-nodeX[row][col][mom].elements[j].x;
			dy = y-nodeX[row][col][mom].elements[j].y;
			dz = z-nodeX[row][col][mom].elements[j].z;
			dis = dx*dx+dy*dy+dz*dz;
			if (dis <= dd_max){
				*(SS + (int)(sqrt(dis)*ds)) += 2;
			}
			}
		}
		// Distancias entre puntos del diferente nodo:
		//==================================================
		
		x1D = nodeX[row][col][mom].nodepos.x;
		y1D = nodeX[row][col][mom].nodepos.y;
		z1D = nodeX[row][col][mom].nodepos.z;
				
		//=======================================
		u = row;
		v = col;
		//=======================
		// Nodo 2 movil en Z:
		//=======================
		for ( w = mom + 1;  w < partitions ; w ++){	
		x2D = nodeX[u][v][w].nodepos.x;
		y2D = nodeX[u][v][w].nodepos.y;
		z2D = nodeX[u][v][w].nodepos.z;
		dx_nod = x1D-x2D;
		dy_nod = y1D-y2D;
		dz_nod = z1D-z2D;
		dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
		if (dis_nod <= ddmax_nod){
			for ( i = 0; i < nodeX[row][col][mom].len; i++){
			x = nodeX[row][col][mom].elements[i].x;
			y = nodeX[row][col][mom].elements[i].y;
			z = nodeX[row][col][mom].elements[i].z;
				for ( j = 0; j < nodeX[u][v][w].len; j++){
				dx =  x-nodeX[u][v][w].elements[j].x;
				dy =  y-nodeX[u][v][w].elements[j].y;
				dz =  z-nodeX[u][v][w].elements[j].z;
				dis = dx*dx+dy*dy+dz*dz;
				if (dis <= dd_max){
					*(SS + (int)(sqrt(dis)*ds)) += 2;
				}
				}
		}
		}
		// Distacia de los puntos frontera DD
		//=======================================
					
		//Condiciones de nodos en frontera:
		con_x = ((x1D<=d_max_pm)&&(x2D>=front_pm))||((x2D<=d_max_pm)&&(x1D>=front_pm));
		con_y = ((y1D<=d_max_pm)&&(y2D>=front_pm))||((y2D<=d_max_pm)&&(y1D>=front_pm));
		con_z = ((z1D<=d_max_pm)&&(z2D>=front_pm))||((z2D<=d_max_pm)&&(z1D>=front_pm));
				
		if(con_x || con_y || con_z){
			histo_front_XX(SS,nodeX,dis_nod,fabs(dx_nod),fabs(dy_nod),fabs(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w);
		}
		}
		
		//=======================
		// Nodo 2 movil en ZY:
		//=======================
		for (v = col + 1; v < partitions ; v ++){
		for (w = 0; w < partitions ; w ++){
		x2D = nodeX[u][v][w].nodepos.x;
		y2D = nodeX[u][v][w].nodepos.y;
		z2D = nodeX[u][v][w].nodepos.z;
		dx_nod = x1D-x2D;
		dy_nod = y1D-y2D;
		dz_nod = z1D-z2D;
		dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
		if (dis_nod <= ddmax_nod){
			for ( i = 0; i < nodeX[row][col][mom].len; i++){
			x = nodeX[row][col][mom].elements[i].x;
			y = nodeX[row][col][mom].elements[i].y;
			z = nodeX[row][col][mom].elements[i].z;
				for ( j = 0; j < nodeX[u][v][w].len; j++){	
				dx =  x-nodeX[u][v][w].elements[j].x;
				dy =  y-nodeX[u][v][w].elements[j].y;
				dz =  z-nodeX[u][v][w].elements[j].z;
				dis = dx*dx+dy*dy+dz*dz;
				if (dis <= dd_max){
					*(SS + (int)(sqrt(dis)*ds)) += 2;
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
		histo_front_XX(SS,nodeX,dis_nod,fabs(dx_nod),fabs(dy_nod),fabs(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w);
		}
		}
		}
		
		//=======================
		// Nodo 2 movil en ZYZ:
		//=======================
				for ( u = row + 1; u < partitions; u++){
					for ( v = 0; v < partitions; v++){
						for ( w = 0; w < partitions; w++){
							// Histograma DD
							x2D = nodeX[u][v][w].nodepos.x;
							y2D = nodeX[u][v][w].nodepos.y;
							z2D = nodeX[u][v][w].nodepos.z;
							dx_nod = x1D-x2D;
							dy_nod = y1D-y2D;
							dz_nod = z1D-z2D;
							dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
							if (dis_nod <= ddmax_nod){
								for ( i = 0; i < nodeX[row][col][mom].len; i++){
									x = nodeX[row][col][mom].elements[i].x;
									y = nodeX[row][col][mom].elements[i].y;
									z = nodeX[row][col][mom].elements[i].z;
									for ( j = 0; j < nodeX[u][v][w].len; j++){	
									dx =  x-nodeX[u][v][w].elements[j].x;
									dy =  y-nodeX[u][v][w].elements[j].y;
									dz =  z-nodeX[u][v][w].elements[j].z;
									dis = dx*dx + dy*dy + dz*dz;
									if (dis <= dd_max){
										*(SS + (int)(sqrt(dis)*ds)) += 2;
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
				histo_front_XX(SS,nodeX,dis_nod,fabs(dx_nod),fabs(dy_nod),fabs(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w);
				}
							
						}	
					}
				}
			}
		}
	}
	#pragma omp critical
        {
        	for(int a=0; a<bn; a++) {
        		*(XX+a)+=*(SS+a);
        	}
        }
	}
	// Histograma RR (ANALITICA)
	//======================================
	float dr = (d_max/bn);
	float alph = 4*(3.14159265359)*(n_pts*n_pts)*dr*dr*dr/(3*size_box*size_box*size_box);
	float r1;
	for(int a=0; a<bn; a++) {
		r1 = (a+1);
        	*(YY+a) += alph*(r1*r1*r1-a*a*a);
	}
}
//=================================================================== 

void NODE::make_histoXY(unsigned int *XY, Node ***nodeX, Node ***nodeY){
	/*
	Función para crear el histograma DR. 
	
	Argumentos
	DR: arreglo donde se creará el histograma DR.
	
	*/
	int i, j, u, v, w, partitions = (int)((size_box/size_node)+1);
	float dis, dis_nod;
	float pt_m = size_node/2;
	std::cout << "-> Estoy haciendo histograma XY..." << std::endl;
	
	#pragma omp parallel num_threads(4) private(i,j,u,v,w,dis,dis_nod)
    	{
    	unsigned int *SS;
    	SS = new unsigned int[bn];
    	for (int k = 0; k < bn; k++){
		*(SS+k) = 0;
	}
	//otras variables privadas:
	float x1D, y1D, z1D, x2R, y2R, z2R;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod;
	bool con_x, con_y, con_z;
	#pragma omp for  collapse(3)  schedule(dynamic) //private(i,j,u,v,w,dis,dis_nod)
	for (int row = 0; row < partitions; row++){
		for (int col = 0; col < partitions; col++){
			for (int mom = 0; mom < partitions; mom++){
				// Distancias entre puntos de diferentes nodos de diferentes datos
				x1D = nodeX[row][col][mom].nodepos.x;
				y1D = nodeX[row][col][mom].nodepos.y;
				z1D = nodeX[row][col][mom].nodepos.z;
				for ( u = 0; u < partitions; u++){
					for ( v = 0; v < partitions; v++){
						for ( w = 0; w < partitions; w++){
							// Histograma DR
							x2R = nodeY[u][v][w].nodepos.x;
							y2R = nodeY[u][v][w].nodepos.y;
							z2R = nodeY[u][v][w].nodepos.z;
							
							dx_nod = x1D-x2R;
							dy_nod = y1D-y2R;
							dz_nod = z1D-z2R;
							
							dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
							if (dis_nod <= ddmax_nod){
							for ( i = 0; i < nodeD[row][col][mom].len; i++){
							for ( j = 0; j < nodeR[u][v][w].len; j++){	
								dx =  nodeX[row][col][mom].elements[i].x-nodeY[u][v][w].elements[j].x;
								dy =  nodeX[row][col][mom].elements[i].y-nodeY[u][v][w].elements[j].y;
								dz =  nodeX[row][col][mom].elements[i].z-nodeY[u][v][w].elements[j].z;
								dis = dx*dx + dy*dy + dz*dz;
								if (dis < dd_max){
									*(SS+ (int)(sqrt(dis)*ds)) += 1;
								}
							}
							}	
							}
							// Distacia de los puntos frontera
							//=======================================
					
							//Condiciones de nodos en frontera:
					con_x = ((x1D-pt_m)<d_max && (x2R+pt_m)>front)||((x2R-pt_m)<d_max && x1D+pt_m>front);
					con_y = ((y1D-pt_m)<d_max && (y2R+pt_m)>front)||((y2R-pt_m)<d_max && y1D+pt_m>front);
					con_z = ((z1D-pt_m)<d_max && (z2R+pt_m)>front)||((z2R-pt_m)<d_max && z1D+pt_m>front);
						
				if(con_x || (con_y || con_z)){
				histo_front_XY(SS,nodeX,nodeY,dis_nod,fabs(dx_nod),fabs(dy_nod),fabs(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w);
				}
							
						}
					}	
				}
			}
		}
	}
	#pragma omp critical
        {
            for( int a=0; a<bn; a++ ) {
                 *(XY+a)+=*(SS+a);
            }
        }
	}
}

//=================================================================== 
void NODE::histo_front_XX(unsigned int *PP, Node ***dat, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int _row, int _col, int _mom, int _u, int _v, int _w){
	int i, j;
	float dis_f,_dis,_d_x,_d_y,_d_z;
	float _x,_y,_z;
	//======================================================================
	// Si los puentos estás en las paredes laterales de X
	if( con_in_x ){
		// forma de calcular la distancia a las proyecciones usando la distancia entre puntos dentro de la caja
		dis_f = disn + ll - 2*dn_x*size_box;
		if (dis_f <= ddmax_nod){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
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
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
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
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
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
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
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
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
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
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
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
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				_x = dat[_row][_col][_mom].elements[i].x;
				_y = dat[_row][_col][_mom].elements[i].y;
				_z = dat[_row][_col][_mom].elements[i].z;
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
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
//=================================================================== 

void NODE::histo_front_XY(unsigned int *XY, Node ***dat, Node ***ran, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int _row, int _col, int _mom, int _u, int _v, int _w){
	/*
	Función para contar las distancias a las fronteras 
	por codiciones periodicas de frontera
	*/
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


NODE::~NODE(){
	
}
