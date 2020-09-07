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
		float d_max;
		Point3D *dataD;
		// Derivados
		float dd_max;
		float ds;
	
	// Métodos de Clase:
	public:
		//Constructor de clase:
		NODE3P(int _bn, int _n_pts, float _size_box, float _d_max, Point3D *_dataD){
			bn = _bn;
			n_pts = _n_pts;
			size_box = _size_box;
			d_max = _d_max;
			dataD = _dataD;
			dd_max = d_max*d_max;
			ds = ((float)(bn))/d_max;
		}
		
		// Implementamos Método de mallas:
		void make_histoXXX(unsigned int ***);
		void symmetrize(unsigned int ***);
		//void histo_front_XX(unsigned int *, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int);
		//void histo_front_XY(unsigned int *, Node ***, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int);
		~NODE3P();
};

//=================================================================== 
//==================== Funciones ==================================== 
//=================================================================== 

void NODE3P::make_histoXXX(unsigned int ***XXX){
	/*
	Función para crear los histogramas DDD y RRR.
	
	Argumentos
	DDD: arreglo donde se creará el histograma DDD.
	RRR: arreglo donde se creará el histograma RRR.
	
	*/ 
	int i, j, k, row, col, mom, u, v, w, a ,b, c;
	float d12, d13, d23;
	float x1N, y1N, z1N, x2N, y2N, z2N, x3N, y3N, z3N;
	float x1, y1, z1, x2, y2, z2, x3, y3, z3;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod, dx_nod2, dy_nod2;
	bool con_x, con_y, con_z;
	//float d_max_pm = d_max + size_node/2, front_pm = front - size_node/2;
	
	std::cout << "-> Estoy haciendo histograma XXX..." << std::endl;

	for(i=0; i<n_pts-2; i++){
		x1 = dataD[i].x;
		y1 = dataD[i].y;
		z1 = dataD[i].z;
		for(j=i+1; j<n_pts-1; j++){
			x2 = dataD[j].x;
			y2 = dataD[j].y;
			z2 = dataD[j].z;
			dx = x2-x1;
			dy = y2-y1;
			dz = z2-z1;
			d12 = dx*dx + dy*dy + dz*dz;
			if(d12 <= dd_max){
			for(k=j+1; k<n_pts; k++){
				x3 = dataD[k].x;
				y3 = dataD[k].y;
				z3 = dataD[k].z;
				dx = x3-x1;
				dy = y3-y1;
				dz = z3-z1;
				d13 = dx*dx + dy*dy + dz*dz;
				if(d13 <= dd_max){
				dx = x3-x2;
				dy = y3-y2;
				dz = z3-z2;
				d23 = dx*dx + dy*dy + dz*dz;
				if(d23 <= dd_max){
				*(*(*(XXX+(int)(sqrt(d12)*ds))+(int)(sqrt(d13)*ds))+(int)(sqrt(d23)*ds))+=1;
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
/*
void NODE3P::histo_front_XXX(unsigned int *PP, Node ***dat, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int _row, int _col, int _mom, int _u, int _v, int _w){
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
						*(PP + (int)(sqrt(_dis)*ds)) += 2;
					}
				}
			}
		}
	}
}
*/
//=================================================================== 
NODE3P::~NODE3P(){
	
}
