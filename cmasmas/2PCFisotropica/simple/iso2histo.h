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
		Point3D *data;	// Datos.
		Point3D *rand;	// Puntos aleatorios.
	// Métodos de la Clase.
	public:
		//Constructor de la clase.
		iso2hist(int _bin, int _n_pts, float _d_max, Point3D *_data, Point3D *_rand){
			bin = _bin;
			n_pts = _n_pts;
			d_max = _d_max;
			data = _data;
			rand = _rand;
		}
		// Métodos para que usuario ingrece variables.
		void setBins(int _bin){
			bin = _bin;
		}
		void setNpts(int _n_pts){
			n_pts = _n_pts;
		}
		void setDmax(float _d_max){
			d_max = _d_max; 
		}
		void setData(Point3D *_data){
			data = _data;
		}
		void setRand(Point3D *_rand){
			rand = _rand;
		}
		// Método para obtener las variable ingresadas anteriormente.
		int getBins(){
			return bin;
		}
		int getNpts(){
			return n_pts;
		}
		float getDmax(){
			return d_max;
		}
		// Métodos para hacer histogramas.
		void make_histoXX(unsigned int *DD, unsigned int *RR){
			int pos; // Posición de apuntador.
			float dis, ds = (float)(bin)/d_max, dd_max = d_max*d_max, dx, dy, dz;
			std::cout << "Estoy haciendo histogramas DD y RR..." << std::endl;
			for(int i = 0; i < n_pts-1; i++){
				for(int j = i+1; j < n_pts; j++){
					dx = data[i].x-data[j].x;
					dy = data[i].y-data[j].y;
					dz = data[i].z-data[j].z;
					dis = dx*dx + dy*dy + dz*dz;
					if(dis <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						DD[pos] += 2;
					}
					dx = rand[i].x-rand[j].x;
					dy = rand[i].y-rand[j].y;
					dz = rand[i].z-rand[j].z;
					dis = dx*dx + dy*dy + dz*dz;
					if(dis <= dd_max){
						pos = (int)(sqrt(dis)*ds);
						RR[pos] +=2;
					}
				}
			}
		}
		void make_histoXY(unsigned int *DR){
			int pos;
			float dis, ds = (float)(bin)/d_max, dd_max = d_max*d_max, dx, dy, dz;
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
				}
			}
		}
		~iso2hist(){ // Destructor
		
		}
};
