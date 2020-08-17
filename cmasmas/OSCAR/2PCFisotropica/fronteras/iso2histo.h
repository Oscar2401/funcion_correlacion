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
			std::cout << "Estoy haciendo histogramas DD y RR..." << std::endl;
			int c = 0;
			for(int i = 0; i < n_pts-1; i++){
				for(int j = i+1; j < n_pts; j++){
					//DATA
					//====================================
					dx = data[i].x-data[j].x;
					dy = data[i].y-data[j].y;
					dz = data[i].z-data[j].z;
					dis = dx*dx + dy*dy + dz*dz;
					if(dis <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DD[pos] += 2;
					}
					
					con_x = (data[i].x < d_max || data[j].x < d_max)&&(data[i].x > front || data[j].x > front);
					con_y = (data[i].y < d_max || data[j].y < d_max)&&(data[i].y > front || data[j].y > front);
					con_z = (data[i].z < d_max || data[j].z < d_max)&&(data[i].z > front || data[j].z > front);
					
					// Distancias en frontera 
					if( con_x ){
					dis_f = dis + ll - 2*abs(dx)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DD[pos] += 4;
					}
					}
					
					if( con_y ){
					dis_f = dis + ll - 2*abs(dy)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DD[pos] += 4;
					}
					}
					
					if( con_z ){
					dis_f = dis + ll - 2*abs(dz)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DD[pos] += 4;
					}
					}
					
					if( con_x && con_y ){
					dis_f = dis + 2*ll - 2*(abs(dx)+abs(dy))*size_box;
					if(dis <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DD[pos] += 4;
					}
					}
					
					if( con_x && con_z ){
					dis_f = dis + 2*ll - 2*(abs(dx)+abs(dz))*size_box;
					if(dis <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DD[pos] += 4;
					}
					}
					
					if( con_y && con_z ){
					dis_f = dis + 2*ll - 2*(abs(dy)+abs(dz))*size_box;
					if(dis <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DD[pos] += 4;
					}
					}
					
					if( con_x && con_y && con_z ){
					dis_f = dis + 3*ll - 2*(abs(dx)+ abs(dy) + abs(dz))*size_box;
					if(dis <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DD[pos] += 4;
					}
					}
					// RANDOM
					//====================================
					dx = rand[i].x-rand[j].x;
					dy = rand[i].y-rand[j].y;
					dz = rand[i].z-rand[j].z;
					dis = dx*dx + dy*dy + dz*dz;
					if(dis <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						RR[pos] +=2;
					}
					
					
					con_x = (rand[i].x < d_max || rand[j].x < d_max)&&(rand[i].x > front || rand[j].x > front);
					con_y = (rand[i].y < d_max || rand[j].y < d_max)&&(rand[i].y > front || rand[j].y > front);
					con_z = (rand[i].z < d_max || rand[j].z < d_max)&&(rand[i].z > front || rand[j].z > front);
					
					// Distancias en frontera 
					if( con_x ){
					dis_f = dis + ll - 2*abs(dx)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						RR[pos] += 4;
					}
					}
					
					if( con_y ){
					dis_f = dis + ll - 2*abs(dy)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						RR[pos] += 4;
					}
					}
					
					if( con_z ){
					dis_f = dis + ll - 2*abs(dz)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						RR[pos] += 4;
					}
					}
					
					if( con_x && con_y ){
					dis_f = dis + 2*ll - 2*(abs(dx)+abs(dy))*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						RR[pos] += 4;
					}
					}
					
					if( con_x && con_z ){
					dis_f = dis + 2*ll - 2*(abs(dx)+abs(dz))*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						RR[pos] += 4;
					}
					}
					
					if( con_y && con_z ){
					dis_f = dis + 2*ll - 2*(abs(dy)+abs(dz))*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						RR[pos] += 4;
					}
					}
					
					if( con_x && con_y && con_z ){
					dis_f = dis + 3*ll - 2*(abs(dx)+abs(dy)+abs(dz))*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						RR[pos] += 4;
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
			std::cout << "Estoy haciendo histograma DR..." << std::endl;
			for (int i = 0; i < n_pts; i++){
				for(int j = 0; j < n_pts; j++){
					dx = data[i].x-rand[j].x;
					dy = data[i].y-rand[j].y;
					dz = data[i].z-rand[j].z;
					dis = dx*dx + dy*dy + dz*dz;
					if(dis <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DR[pos] += 1;
					}
					
					con_x = (rand[i].x < d_max || rand[j].x < d_max)&&(rand[i].x > front || rand[j].x > front);
					con_y = (rand[i].y < d_max || rand[j].y < d_max)&&(rand[i].y > front || rand[j].y > front);
					con_z = (rand[i].z < d_max || rand[j].z < d_max)&&(rand[i].z > front || rand[j].z > front);
					
					// Distancias en frontera 
					if( con_x ){
					dis_f = dis + ll - 2*abs(dx)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DR[pos] += 2;
					}
					}
					
					if( con_y ){
					dis_f = dis + ll - 2*abs(dy)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DR[pos] += 2;
					}
					}
					
					if( con_z ){
					dis_f = dis + ll - 2*abs(dz)*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DR[pos] += 2;
					}
					}
					
					if( con_x && con_y ){
					dis_f = dis + 2*ll - 2*(abs(dx)+abs(dy))*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DR[pos] += 2;
					}
					}
					
					if( con_x && con_z ){
					dis_f = dis + 2*ll - 2*(abs(dx)+abs(dz))*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DR[pos] += 2;
					}
					}
					
					if( con_y && con_z ){
					dis_f = dis + 2*ll - 2*(abs(dy)+abs(dz))*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DR[pos] += 2;
					}
					}
					
					if( con_x && con_y && con_z ){
					dis_f = dis + 3*ll - 2*(abs(dx)+abs(dy)+abs(dz))*size_box;
					if(dis_f <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DR[pos] += 2;
					}
					}
				}
			}
		}
		~iso2hist(){ // Destructor
		
		}
};
