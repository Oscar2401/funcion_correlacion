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
		void make_histoXX(float *DD, float *RR){
			int pos; // Posición de apuntador.
			float dis, aux, ds;
			ds = (float)(bin)/d_max;
			for(int i=0; i<n_pts; i++){
				for(int j=i+1; j<n_pts; j++){
					dis = sqrt((data[i].x-data[j].x)*(data[i].x-data[j].x) + (data[i].y-data[j].y)*(data[i].y-data[j].y) + (data[i].z-data[j].z)*(data[i].z-data[j].z));
					if(dis < d_max){
						pos = (int)(dis*ds);
						DD[pos] += 2;
					}
					dis = sqrt((rand[i].x-rand[j].x)*(rand[i].x-rand[j].x)+(rand[i].y - rand[j].y)*(rand[i].y - rand[j].y)+(rand[i].z - rand[j].z)*(rand[i].z - rand[j].z));
					if(dis < d_max){
						pos = (int)(dis*ds);
						RR[pos] +=2;
					}
				}
			}
		}
		void make_histoXY(float *DR){
			int pos;
			float dis, aux, ds;
			ds = (float)(bin)/d_max;
			for (int i=0; i<n_pts; i++){
				for(int j=0; j<n_pts; j++){
					dis = sqrt((data[i].x - rand[j].x)*(data[i].x - rand[j].x) + (data[i].y - rand[j].y)*(data[i].y - rand[j].y) + (data[i].z - rand[j].z)*(data[i].z - rand[j].z));
					if(dis<d_max){
						pos = (int)(dis*ds);
						DR[pos] += 1;
					}
				}
			}
		}
		~iso2hist(){ // Destructor
		
		}
};
