#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <omp.h>

/*
 * Code for an octree that demonstrates insertion and search
 */
#include <iostream>
#include <vector>

using namespace std;

// Octantes de los nodos cubicos:
#define TLF 0    // top left front
#define TRF 1    // top right front
#define BRF 2    // bottom right front
#define BLF 3    // bottom left front
#define TLB 4    // top left back
#define TRB 5    // top right back
#define BRB 6    // bottom right back
#define BLB 7    // bottom left back

// Estructura de punto:
struct Point{
	float x;
	float y;
	float z;

	Point() : x(-1), y(-1), z(-1) {}
	
	Point(float a, float b, float c) : x(a), y(b), z(c) {}
};

// Clase:
class Octree{
	// if point == NULL, node is regional.
	// if point == (-1, -1, -1), node is empty.
	Point *point;
	
	// Representacion del espacio.
	Point *top_left_front
	Point *bottom_right_back;   
	
	std::vector<Octree *> children;

public:
	// Constructor
	Octree(){
		// Para declarar nodos vacios
		point = new Point();
    	}
	// Constructor con tres argumentos
	Octree(float x, float y, float z){
		// para declarar point node
		point = new Point(x, y, z);
	}
	// Constructor con seis argumentos
	Octree(float x1, float y1, float z1, float x2, float y2, float z2){
		/*
		 * Este usado para construir Octree  
		 * con límites definidos
		 */
		 
		if(x2 < x1 || y2 < y1 || z2 < z1){
		cout << "bounday poitns are not valid" << endl; 
		return;
		}
		
		// Asignamos nodos vacios para los hijos
		point = nullptr;
		top_left_front = new Point(x1, y1, z1);
		bottom_right_back = new Point(x2, y2, z2);
		children.assign(8, nullptr); // Ocho hijos (ocho cotantes divididos)
		
		for(int i = TLF; i <= BLB; ++i) children[i] = new Octree();
		
	}
	
	// Método de busqueda de puntos:
	void insert(float, float, float);
	bool find(float, float, float);
	~Octree();
};

void Octree::insert(float x, float y, float z){
	/*
	 *
	 * Función para insertar un punto en el octree.
	 *
	 */
	
	// Si el punto ya existe en el octree, este ya no se incluye:
	if(find(x, y, z)){
	cout << "Point already exist in the tree" << endl; 
	return;
        }
        
        // Definimos el punto medio del cubo raiz: 
	float midx = (top_left_front->x + bottom_right_back->x)/2;
	float midy = (top_left_front->y + bottom_right_back->y)/2;
	float midz = (top_left_front->z + bottom_right_back->z)/2;
	
	float pos = -1;
	
	// Busqueda para insertar el punto en uno de los octantes:
	if(x <= midx){
		if(y <= midy){
			if(z <= midz) pos = TLF;
			else pos = TLB;
		}
		else{
			if(z <= midz) pos = BLF;
			else pos = BLB;
		}
	}
	else{
		if(y <= midy){
			if(z <= midz) pos = TRF;
			else pos = TRB;
		}
		else{
			if(z <= midz) pos = BRF;
			else pos = BRB;
		}
	}
	
	// Si un nodo interno es encontrado:
	if(children[pos]->point == nullptr){
		// if region node
		children[pos]->insert(x, y, z);
		return;
	}
	// Si un nodo vacio es encontrado:
	else if(children[pos]->point->x == -1){
		// if empty node
		delete children[pos];
		children[pos] = new Octree(x, y, z);
		return;
	}
	else{
		float x_ = children[pos]->point->x;
		float y_ = children[pos]->point->y;
		float z_ = children[pos]->point->z;
		delete children[pos];
		children[pos] = nullptr;
		if(pos == TLF){
			children[pos] = new Octree(top_left_front->x, top_left_front->y, top_left_front->z, midx, midy, midz);
		}
		else if(pos == TRF){
			children[pos] = new Octree(midx + 1, top_left_front->y, top_left_front->z, bottom_right_back->x, midy, midz);
		}
		else if(pos == BRF){
			children[pos] = new Octree(midx + 1, midy + 1, top_left_front->z, bottom_right_back->x, bottom_right_back->y, midz);
		}
		else if(pos == BLF){
			children[pos] = new Octree(top_left_front->x, midy + 1, top_left_front->z, midx, bottom_right_back->y, midz);
		}
		else if(pos == TLB){
			children[pos] = new Octree(top_left_front->x, top_left_front->y, midz + 1, midx, midy, bottom_right_back->z);
		}
		else if(pos == TRB){
			children[pos] = new Octree(midx + 1, top_left_front->y, midz + 1, bottom_right_back->x, midy, bottom_right_back->z);
		}
		else if(pos == BRB){
			children[pos] = new Octree(midx + 1, midy + 1, midz + 1, bottom_right_back->x, bottom_right_back->y, bottom_right_back->z);
		}
		else if(pos == BLB){
			children[pos] = new Octree(top_left_front->x, midy + 1, midz + 1, midx, bottom_right_back->y, bottom_right_back->z);
		}
		children[pos]->insert(x_, y_, z_);
		children[pos]->insert(x, y, z);
	}
}
    
bool Octree::find(float x, float y, float z){
	/*
	 * Función que regresa true si el punto
	 * (x, y, z) existe en el octree.
	 */
	
	// Si el punto está fuera del límite:
	if(x < top_left_front->x || x > bottom_right_back->x || y < top_left_front->y || y > bottom_right_back->y || z < top_left_front->z || z > bottom_right_back->z) return 0;
	
	float midx = (top_left_front->x + bottom_right_back->x)/2;
	float midy = (top_left_front->y + bottom_right_back->y)/2;
	float midz = (top_left_front->z + bottom_right_back->z)/2;

	float pos = -1;

	if(x <= midx){
		if(y <= midy){
			if(z <= midz) pos = TLF;
			else pos = TLB;
		}
		else{
			if(z <= midz) pos = BLF;
			else pos = BLB;
		}
	}
	else{
		if(y <= midy){
			if(z <= midz) pos = TRF;
			else pos = TRB;
		}
		else{
			if(z <= midz) pos = BRF;
			else pos = BRB;
		}
	}
	// Si un nodo interno es encontrado:
	if(children[pos]->point == nullptr){
		// if region node
		return children[pos]->find(x, y, z);
	}
	// Si un nodo vacio es encontrado:
	else if(children[pos]->point->x == -1){
		// if empty node
		return 0;
	}
	else{
		if(x == children[pos]->point->x && y == children[pos]->point->y && z == children[pos]->point->z) return 1;
	}
	return 0;
}

//=================================================================== 
Octree::~Octree(){	
}


