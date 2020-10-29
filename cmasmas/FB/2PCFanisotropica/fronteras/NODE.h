
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
		void make_histoXX(unsigned int *, Node ***);
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
	int i, row, col, mom, partitions = (int)(ceil(size_box/size_node));
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

void NODE::make_histoXX(unsigned int *XX, Node ***nodeX){
	/*
	Función para crear los histogramas DD y RR.
	
	Argumentos
	DD: arreglo donde se creará el histograma DD.
	RR: arreglo donde se creará el histograma RR.
	
	*/
	int i, j, u, v, w, partitions = (int)(ceil(size_box/size_node));
	
	float dis, dis_nod;
	float pt_m = size_node/2;
	
	std::cout << "-> Estoy haciendo histograma XX..." << std::endl;
	
	#pragma omp parallel num_threads(4) private(i,j,u,v,w,dis,dis_nod) 
    	{
    	unsigned int *SS;
    	SS = new unsigned int[bn];
    	for (int k = 0; k < bn; k++){
		*(SS+k) = 0;
	}
    	#pragma omp for collapse(3)  schedule(dynamic) private(i,j,u,v,w,dis,dis_nod)
	for (int row = 0; row < partitions; row++){
		for (int col = 0; col < partitions; col++){
			for (int mom = 0; mom < partitions; mom++){
				
				//otras variables privadas:
				float x1D, y1D, z1D, x2D, y2D, z2D;
				float dx, dy, dz, dx_nod, dy_nod, dz_nod;
				bool con_x, con_y, con_z;
				
				// Distancias entre puntos del mismo nodo:
				//==================================================
				// Histograma DD
				for ( i= 0; i < nodeD[row][col][mom].len - 1; i++){
					for ( j = i+1; j < nodeD[row][col][mom].len; j++){
						dx = nodeD[row][col][mom].elements[i].x-nodeD[row][col][mom].elements[j].x;
						dy = nodeD[row][col][mom].elements[i].y-nodeD[row][col][mom].elements[j].y;
						r_pll = abs(nodeD[row][col][mom].elements[i].z-nodeD[row][col][mom].elements[j].z);
						if (r_pll <= d_max){
							r_ort = dx*dx + dy*dy 
							if(r_ort <= dd_max) 
								*(*(SS+(int)(r_pll*ds))+(int)(sqrt(r_ort)*ds))+=2;
						}
					}
				}
				// Distancias entre puntos del diferente nodo:
				//==================================================
				
				x1D = nodeD[row][col][mom].nodepos.x;
				y1D = nodeD[row][col][mom].nodepos.y;
				z1D = nodeD[row][col][mom].nodepos.z;
				
				//=======================================
				u = row;
				v = col;
				for ( w = mom + 1;  w < partitions ; w ++){	
					// Histograma DD
					x2D = nodeD[u][v][w].nodepos.x;
					y2D = nodeD[u][v][w].nodepos.y;
					z2D = nodeD[u][v][w].nodepos.z;
					dx_nod = x1D-x2D;
					dy_nod = y1D-y2D;
					dz_nod = z1D-z2D;
					dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
					if (sqrt(dis_nod)-corr <= d_max){
						for ( i = 0; i < nodeD[row][col][mom].len; i++){
							for ( j = 0; j < nodeD[u][v][w].len; j++){
							dx =  nodeD[row][col][mom].elements[i].x-nodeD[u][v][w].elements[j].x;
							dy =  nodeD[row][col][mom].elements[i].y-nodeD[u][v][w].elements[j].y;
							dz =  nodeD[row][col][mom].elements[i].z-nodeD[u][v][w].elements[j].z;
							dis = dx*dx + dy*dy + dz*dz;
							if (dis <= dd_max){
								*(SS + (int)(sqrt(dis)*ds)) += 2;
							}
							}
						}
					}
					// Distacia de los puntos frontera DD
					//=======================================
					
					//Condiciones de nodos en frontera:
					con_x = (x1D-pt_m<d_max && x2D+pt_m>front)||(x2D-pt_m<d_max && x1D+pt_m>front);
					con_y = (y1D-pt_m<d_max && y2D+pt_m>front)||(y2D-pt_m<d_max && y1D+pt_m>front);
					con_z = (z1D-pt_m<d_max && z2D+pt_m>front)||(z2D-pt_m<d_max && z1D+pt_m>front);
				
					if(con_x || (con_y || con_z) ){
					histo_front_XX(SS,nodeX,dis_nod,abs(dx_nod),abs(dy_nod),abs(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w);
					}
				}
				//=======================================
				for (v = col + 1; v < partitions ; v ++){
					for (w = 0; w < partitions ; w ++){		
						// Histograma DD
						x2D = nodeD[u][v][w].nodepos.x;
						y2D = nodeD[u][v][w].nodepos.y;
						z2D = nodeD[u][v][w].nodepos.z;
						dx_nod = x1D-x2D;
						dy_nod = y1D-y2D;
						dz_nod = z1D-z2D;
						dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
						if (sqrt(dis_nod)-corr <= d_max){
							for ( i = 0; i < nodeD[row][col][mom].len; i++){
								for ( j = 0; j < nodeD[u][v][w].len; j++){	
								dx =  nodeD[row][col][mom].elements[i].x-nodeD[u][v][w].elements[j].x;
								dy =  nodeD[row][col][mom].elements[i].y-nodeD[u][v][w].elements[j].y;
								dz =  nodeD[row][col][mom].elements[i].z-nodeD[u][v][w].elements[j].z;
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
						con_x = (x1D-pt_m<d_max && x2D+pt_m>front)||(x2D-pt_m<d_max && x1D+pt_m>front);
						con_y = (y1D-pt_m<d_max && y2D+pt_m>front)||(y2D-pt_m<d_max && y1D+pt_m>front);
						con_z = (z1D-pt_m<d_max && z2D+pt_m>front)||(z2D-pt_m<d_max && z1D+pt_m>front);
					
					if(con_x || (con_y || con_z)){ 
					histo_front_XX(SS,nodeX,dis_nod,abs(dx_nod),abs(dy_nod),abs(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w);
					}
					}
				}
				//=======================================
				for ( u = row + 1; u < partitions; u++){
					for ( v = 0; v < partitions; v++){
						for ( w = 0; w < partitions; w++){
							// Histograma DD
							x2D = nodeD[u][v][w].nodepos.x;
							y2D = nodeD[u][v][w].nodepos.y;
							z2D = nodeD[u][v][w].nodepos.z;
							dx_nod = x1D-x2D;
							dy_nod = y1D-y2D;
							dz_nod = z1D-z2D;
							dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
							if (sqrt(dis_nod)-corr <= d_max){
								for ( i = 0; i < nodeD[row][col][mom].len; i++){
									for ( j = 0; j < nodeD[u][v][w].len; j++){	
									dx =  nodeD[row][col][mom].elements[i].x-nodeD[u][v][w].elements[j].x;
									dy =  nodeD[row][col][mom].elements[i].y-nodeD[u][v][w].elements[j].y;
									dz =  nodeD[row][col][mom].elements[i].z-nodeD[u][v][w].elements[j].z;
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
							con_x = (x1D-pt_m<d_max && x2D+pt_m>front)||(x2D-pt_m<d_max && x1D+pt_m>front);
							con_y = (y1D-pt_m<d_max && y2D+pt_m>front)||(y2D-pt_m<d_max && y1D+pt_m>front);
							con_z = (z1D-pt_m<d_max && z2D+pt_m>front)||(z2D-pt_m<d_max && z1D+pt_m>front);
					
				if(con_x || (con_y || con_z)){
				histo_front_XX(SS,nodeX,dis_nod,abs(dx_nod),abs(dy_nod),abs(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w);
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
                 *(XX+a)+=*(SS+a);
            }
        }
	}
}
//=================================================================== 

void NODE::make_histoXY(unsigned int *XY, Node ***nodeX, Node ***nodeY){
	/*
	Función para crear el histograma DR. 
	
	Argumentos
	DR: arreglo donde se creará el histograma DR.
	
	*/
	int i, j, u, v, w, row, col, mom, partitions = (int)(ceil(size_box/size_node));
	float dx, dy, dz, dx_nod, dy_nod, dz_nod;
	float x1D, y1D, z1D, x2R, y2R, z2R;
	float dis, dis_nod;
	float pt_m = size_node/2;
	bool con_x, con_y, con_z;
	std::cout << "-> Estoy haciendo histograma XY..." << std::endl;
	
	#pragma omp parallel num_threads(4) private(i,j,u,v,w,dis,dis_nod)
    	{
    	unsigned int *SS;
    	SS = new unsigned int[bn];
    	for (int k = 0; k < bn; k++){
		*(SS+k) = 0;
	}
	#pragma omp for  collapse(3)  schedule(dynamic) private(i,j,u,v,w,dis,dis_nod)
	for (row = 0; row < partitions; row++){
		for (col = 0; col < partitions; col++){
			for (mom = 0; mom < partitions; mom++){
				// Distancias entre puntos de diferentes nodos de diferentes datos
				for ( u = 0; u < partitions; u++){
					for ( v = 0; v < partitions; v++){
						for ( w = 0; w < partitions; w++){
							// Histograma DR
							
							x1D = nodeX[row][col][mom].nodepos.x;
							y1D = nodeX[row][col][mom].nodepos.y;
							z1D = nodeX[row][col][mom].nodepos.z;
							x2R = nodeY[u][v][w].nodepos.x;
							y2R = nodeY[u][v][w].nodepos.y;
							z2R = nodeY[u][v][w].nodepos.z;
							
							dx_nod = x1D-x2R;
							dy_nod = y1D-y2R;
							dz_nod = z1D-z2R;
							
							dis_nod = dx_nod*dx_nod + dy_nod*dy_nod + dz_nod*dz_nod;
							if (sqrt(dis_nod)-corr <= d_max){
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
					con_x = (x1D-pt_m<d_max && x2R+pt_m>front)||(x2R-pt_m<d_max && x1D+pt_m>front);
					con_y = (y1D-pt_m<d_max && y2R+pt_m>front)||(y2R-pt_m<d_max && y1D+pt_m>front);
					con_z = (z1D-pt_m<d_max && z2R+pt_m>front)||(z2R-pt_m<d_max && z1D+pt_m>front);
						
				if(con_x || (con_y || con_z)){
				histo_front_XY(SS,nodeX,nodeY,dis_nod,abs(dx_nod),abs(dy_nod),abs(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w);
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
 
void NODE::histo_front_XX(unsigned int *XX, Node ***dat, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int _row, int _col, int _mom, int _u, int _v, int _w){
	int i, j;
	float dis_f, _dis, _d_x, _d_y, _d_z;
	
	if( con_in_x ){
		dis_f = sqrt(disn + ll - 2*dn_x*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-abs(dat[_row][_col][_mom].elements[i].x-dat[_u][_v][_w].elements[j].x);
					_d_y =  dat[_row][_col][_mom].elements[i].y-dat[_u][_v][_w].elements[j].y;
					_d_z =  dat[_row][_col][_mom].elements[i].z-dat[_u][_v][_w].elements[j].z;
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z; 
					if (_dis <= dd_max){
						*(XX+(int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
					
	if( con_in_y ){
		dis_f = sqrt(disn + ll - 2*dn_y*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  dat[_row][_col][_mom].elements[i].x-dat[_u][_v][_w].elements[j].x;
					_d_y =  size_box-abs(dat[_row][_col][_mom].elements[i].y-dat[_u][_v][_w].elements[j].y);
					_d_z =  dat[_row][_col][_mom].elements[i].z-dat[_u][_v][_w].elements[j].z;
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XX+(int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
			
	if( con_in_z ){
		dis_f = sqrt(disn + ll - 2*dn_z*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  dat[_row][_col][_mom].elements[i].x-dat[_u][_v][_w].elements[j].x;
					_d_y =  dat[_row][_col][_mom].elements[i].y-dat[_u][_v][_w].elements[j].y;
					_d_z =  size_box-abs(dat[_row][_col][_mom].elements[i].z-dat[_u][_v][_w].elements[j].z);
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XX+(int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
					
	if( con_in_x && con_in_y ){
		dis_f = sqrt(disn + 2*ll - 2*(dn_x+dn_y)*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-abs(dat[_row][_col][_mom].elements[i].x-dat[_u][_v][_w].elements[j].x);
					_d_y =  size_box-abs(dat[_row][_col][_mom].elements[i].y-dat[_u][_v][_w].elements[j].y);
					_d_z =  dat[_row][_col][_mom].elements[i].z-dat[_u][_v][_w].elements[j].z;
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XX+(int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
					
	if( con_in_x && con_in_z ){
		dis_f = sqrt(disn + 2*ll - 2*(dn_x+dn_z)*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-abs(dat[_row][_col][_mom].elements[i].x-dat[_u][_v][_w].elements[j].x);
					_d_y =  dat[_row][_col][_mom].elements[i].y-dat[_u][_v][_w].elements[j].y;
					_d_z =  size_box-abs(dat[_row][_col][_mom].elements[i].z-dat[_u][_v][_w].elements[j].z);
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XX+(int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
					
	if( con_in_y && con_in_z ){
		dis_f = sqrt(disn + 2*ll - 2*(dn_y+dn_z)*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  dat[_row][_col][_mom].elements[i].x-dat[_u][_v][_w].elements[j].x;
					_d_y =  size_box-abs(dat[_row][_col][_mom].elements[i].y-dat[_u][_v][_w].elements[j].y);
					_d_z =  size_box-abs(dat[_row][_col][_mom].elements[i].z-dat[_u][_v][_w].elements[j].z);
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XX+(int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
				
	if( con_in_x && con_in_y && con_in_z ){
		dis_f = sqrt(disn + 3*ll - 2*(dn_x+dn_y+dn_z)*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-abs(dat[_row][_col][_mom].elements[i].x-dat[_u][_v][_w].elements[j].x);
					_d_y =  size_box-abs(dat[_row][_col][_mom].elements[i].y-dat[_u][_v][_w].elements[j].y);
					_d_z =  size_box-abs(dat[_row][_col][_mom].elements[i].z-dat[_u][_v][_w].elements[j].z);
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XX+(int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
}

//=================================================================== 

void NODE::histo_front_XY(unsigned int *XY, Node ***dat, Node ***ran, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int _row, int _col, int _mom, int _u, int _v, int _w){
	int _pos, i, j;
	float dis_f, _dis, _d_x, _d_y, _d_z;
	
	if( con_in_x ){
		dis_f = sqrt(disn + ll - 2*dn_x*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-abs(dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x);
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
					
	if( con_in_y ){
		dis_f = sqrt(disn + ll - 2*dn_y*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x;
					_d_y =  size_box-abs(dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y);
					_d_z =  dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z;
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
			
	if( con_in_z ){
		dis_f = sqrt(disn + ll - 2*dn_z*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x;
					_d_y =  dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y;
					_d_z =  size_box-abs(dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z);
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
					
	if( con_in_x && con_in_y ){
		dis_f = sqrt(disn + 2*ll - 2*(dn_x+dn_y)*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-abs(dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x);
					_d_y =  size_box-abs(dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y);
					_d_z =  dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z;
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
					
	if( con_in_x && con_in_z ){
		dis_f = sqrt(disn + 2*ll - 2*(dn_x+dn_z)*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-abs(dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x);
					_d_y =  dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y;
					_d_z =  size_box-abs(dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z);
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
					
	if( con_in_y && con_in_z ){
		dis_f = sqrt(disn + 2*ll - 2*(dn_y+dn_z)*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x;
					_d_y =  size_box-abs(dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y);
					_d_z =  size_box-abs(dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z);
					_dis = _d_x*_d_x + _d_y*_d_y + _d_z*_d_z;
					if (_dis <= dd_max){
						*(XY+(int)(sqrt(_dis)*ds)) += 1;
					}
				}
			}
		}
	}
				
	if( con_in_x && con_in_y && con_in_z ){
		dis_f = sqrt(disn + 3*ll - 2*(dn_x+dn_y+dn_z)*size_box)-corr;
		if (dis_f <= d_max){
			for ( i = 0; i < dat[_row][_col][_mom].len; i++){
				for ( j = 0; j < dat[_u][_v][_w].len; j++){
					_d_x =  size_box-abs(dat[_row][_col][_mom].elements[i].x-ran[_u][_v][_w].elements[j].x);
					_d_y =  size_box-abs(dat[_row][_col][_mom].elements[i].y-ran[_u][_v][_w].elements[j].y);
					_d_z =  size_box-abs(dat[_row][_col][_mom].elements[i].z-ran[_u][_v][_w].elements[j].z);
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
