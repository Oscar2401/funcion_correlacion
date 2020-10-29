#include <stdlib.h>
#include <cmath>

// Estructura de las componentes de un punto en 3-dimenciones
struct Point3D{
	float x;
	float y;
	float z;
};

// Definimos la calse iso2hist
class iso2hist{
	// Atributos de la Clase.
	private:
		int bin;		// Número de bins.
		int n_pts;	// Cantidad de puntos.
		float d_max;	// Distancia máxima de histogramas.
		float size_box;
		Point3D *data;	// Datos.
		Point3D *rand;	// Puntos aleatorios.
	// Métodos de la Clase.
	public:
		//Constructor de la clase.
		iso2hist(int _bin, int _n_pts, float _d_max, float _size_box, Point3D *_data, Point3D *_rand){
			bin = _bin;
			n_pts = _n_pts;
			d_max = _d_max;
			size_box = _size_box;
			data = _data;
			rand = _rand;
		}
		// Métodos para hacer histogramas.
		void make_histoXX(unsigned int *DD, unsigned int *RR){
			int pos; // Posición de apuntador.
			float dis, ds = (float)(bin)/d_max, dd_max = d_max*d_max, dx, dy, dz;
			float front = size_box - d_max, dis_f;
			float ll = size_box*size_box;
			bool con_x, con_y, con_z;
			float x1D, x2D, y1D, y2D, z1D, z2D, x1R, x2R, y1R, y2R, z1R, z2R; 
			std::cout << "Estoy haciendo histogramas DD" << std::endl;
			for(int i = 0; i < n_pts-1; i++){
				for(int j = i+1; j < n_pts; j++){
					//DATA
					//====================================
					x1D = data[i].x;
					x2D = data[j].x;
					y1D = data[i].y;
					y2D = data[j].y;
					z1D = data[i].z;
					z2D = data[j].z;
					
					dx = abs(x1D-x2D);
					dy = abs(y1D-y2D);
					dz = abs(z1D-z2D);
					
					dis = dx*dx + dy*dy + dz*dz;
					if(dis <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DD[pos] += 2;
					}
					//Condiciones para fronteras
					con_x = (x1D < d_max && x2D > front)||(x2D < d_max && x1D > front);
					con_y = (y1D < d_max && y2D > front)||(y2D < d_max && y1D > front);
					con_z = (z1D < d_max && z2D > front)||(z2D < d_max && z1D > front);
					// Distancias en frontera 
					if( con_x ){
						dis_f = dis + ll - 2*dx*size_box;
						if(dis_f <= dd_max){
							pos = (int)(sqrt(dis_f)*ds);
							DD[pos] += 2;
						}
					}
					
					if( con_y ){
						dis_f = dis + ll - 2*dy*size_box;
						if(dis_f <= dd_max){
							pos = (int)(sqrt(dis_f)*ds);
							DD[pos] += 2;
							//std::cout << DD[pos] << std::endl;
						}
					}
					
					if( con_z ){
						dis_f = dis + ll - 2*dz*size_box;
						if(dis_f <= dd_max){
							pos = (int)(sqrt(dis_f)*ds);
							DD[pos] += 2;
						}
					}
					
					if( con_x && con_y ){
						dis_f = dis + 2*ll - 2*(dx+dy)*size_box;
						if(dis_f <= dd_max){
							pos = (int)(sqrt(dis_f)*ds);
							DD[pos] += 2;
						}
					}
					
					if( con_x && con_z ){
						dis_f = dis + 2*ll - 2*(dx+dz)*size_box;
						if(dis_f <= dd_max){
							pos = (int)(sqrt(dis_f)*ds);
							DD[pos] += 2;
						}
					}
					
					if( con_y && con_z ){
						dis_f = dis + 2*ll - 2*(dy+dz)*size_box;
						if(dis_f <= dd_max){
							pos = (int)(sqrt(dis_f)*ds);
							DD[pos] += 2;
						}
					}
					
					if( con_x && con_y && con_z ){
						dis_f = dis + 3*ll - 2*(dx+dy+dz)*size_box;
						if(dis_f <= dd_max){
							pos = (int)(sqrt(dis_f)*ds);
							DD[pos] += 2;
						}
					}
				}
			}
		}
		
		void make_histoXY(unsigned int *DR){
			int pos;
			float dis, ds = (float)(bin)/d_max, dd_max = d_max*d_max, dx, dy, dz;
			float front = size_box - d_max, dis_f;
			float ll = size_box*size_box;
			bool con_x, con_y, con_z;
			float x1D, y1D, z1D, x2R, y2R, z2R;
			std::cout << "Estoy haciendo histograma DR..." << std::endl;
			for (int i = 0; i < n_pts; i++){
				for(int j = 0; j < n_pts; j++){
					
					x1D = data[i].x;
					y1D = data[i].y;
					z1D = data[i].z;
					x2R = rand[j].x;
					y2R = rand[j].y;
					z2R = rand[j].z;
					
					dx = abs(x1D-x2R);
					dy = abs(y1D-y2R);
					dz = abs(z1D-z2R);
					dis = dx*dx + dy*dy + dz*dz;
					if(dis <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DR[pos] += 1;
					}
					
					con_x = (x1D<d_max && x2R>front)||(x2R<d_max && x1D>front);
					con_y = (y1D<d_max && y2R>front)||(y2R<d_max && y1D>front);
					con_z = (z1D<d_max && z2R>front)||(z2R<d_max && z1D>front);
					
					// Distancias en frontera 
					if( con_x ){
					dis_f = dis + ll - 2*dx*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DR[pos] += 1;
					}
					}
					
					if( con_y ){
					dis_f = dis + ll - 2*dy*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis_f)*ds);
						DR[pos] += 1;
					}
					}
					
					if( con_z ){
					dis_f = dis + ll - 2*dz*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis_f)*ds);
						DR[pos] += 1;
					}
					}
					
					if( con_x && con_y ){
					dis_f = dis + 2*ll - 2*(dx+dy)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis_f)*ds);
						DR[pos] += 1;
					}
					}
					
					if( con_x && con_z ){
					dis_f = dis + 2*ll - 2*(dx+dz)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis_f)*ds);
						DR[pos] += 1;
					}
					}
					
					if( con_y && con_z ){
					dis_f = dis + 2*ll - 2*(dy+dz)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis_f)*ds);
						DR[pos] += 1;
					}
					}
					
					if( con_x && con_y && con_z ){
					dis_f = dis + 3*ll - 2*(dx+dy+dz)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis_f)*ds);
						DR[pos] += 1;
					}
					}
				}
			}
		}
		~iso2hist(){ // Destructor
		
		}
};
