#include "distancias.h"
#include <stdlib.h> //funciona el new y delete y printf

//Structura que define un punto 2D
struct Punto{
    float x;
    float y;
};



//Definimos la clase:
class iso2PCF_mallas {
	// Atributos:
	private:
		int bin;	// bins de histograma
		int N_point;	// Cantidad de puntos datos
		float d_max;	// Distancia máxima de histograma
		float l_quad;	// Tamaño de cuadros en mallas
		Punto *data;	// Datos 
		Punto *rand;	// Puntos aleatorios
		Punto *mesh_d	// Malla de datos
		Punto *mesh_r	// Mlla de puntos aleatorios 
		Punto *PM_data	// Puntos medios en datos
		Punto *PM_rand // Puntos medios en puntos aleatorios
	// Métodos:
	public:
		iso2PCF_mallas(int _bin, int _N_point, float _d_max, float _l_quad, Punto *_data, Punto *_rand) {			
			bin = _bin;
			N_point = _N_point;
			d_max = _d_max;
			l_quad = _l_quad
			data = _data
			rand = _rand
		}
		//funciones para modificar los atributos 
		void setBin(int _bin){
			bin = _bin
		}
		void setPoints(int _N_point){
			N_point = _N_point;
		}
		void setDMAX(float _d_max){
			d_max = _d_max;
		}
		void setSidequad(int _l_quad){
			l_quad = _l_quad
		}
		void setData(Punto *_data){
			data = _data;
		}
		void setRand(Punto *_rand){
			rand = _rand;
		}
		//funciones para ver los atributos 
		int getBin(){
			return num_bins;
		}
		int getNumPoint(){
			return num_puntos;
		}
		float getDMAX(){
			return d_max;
		}
		void getDATA(){
			for (int i = 0; i < num_puntos; i++){
				printf("%f - %f \n",data[i].x, data[i].y);
			}
		}
		void getRAND(){
			for (int i = 0; i < num_puntos; i++){
				printf("%f - %f \n",rand[i].x, rand[i].y);
			}
		}
		void histogramasPuros(float *DD, float *RR){
		
		int i, j, pos;
		float dd, rr, s, aux;
		
		
		
		}

}


