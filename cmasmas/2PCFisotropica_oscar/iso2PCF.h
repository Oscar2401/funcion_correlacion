#include "distancias.h"
#include <stdlib.h> //funciona el new y delete y printf

//Structura que define un punto 3D
struct Punto{
    float x;
    float y;
    float z;
};

// clase
class iso2PCF{
    // Atributos
    private:
        int num_bins;
        int num_puntos;
        float d_max;
        Punto *data;
        Punto *rand;
    // Métodos
    public:
        //Constructor de la clase
        iso2PCF(int _num_bins, int _num_puntos, float _d_max, Punto *_data, Punto *_rand){
            num_bins = _num_bins;
            num_puntos = _num_puntos;
            d_max = _d_max;
            data = _data;
            rand = _rand;
        }
        void setBins(int _num_bins){
            num_bins =  _num_bins;
        }
        void setDMAX(float _d_max){
            d_max = _d_max;
        }
        void setData(Punto *_data){
            data = _data;
        }
        void setRand(Punto *_rand){
            rand = _rand;
        }
        void setPuntos(int _NPuntos){
            num_puntos = _NPuntos;
        }
        int getBins(){
            return num_bins;
        }
        int getNumPuntos(){
            return num_puntos;
        }
        float getDMAX(){
            return d_max;
        }
        void getDATA(){
            for (int i = 0; i < num_puntos; i++)
            {
                printf("%f - %f - %f \n",data[i].x, data[i].y, data[i].z);
            }
        }
        void getRAND(){
            for (int i = 0; i < num_puntos; i++)
            {
                printf("%f - %f - %f \n",rand[i].x, rand[i].y, rand[i].z);
            }
        }
        void histogramasPuros(float *DD, float *RR){
            int i,j, pos;
            float dd, rr, s, aux;
            s = (float)(num_bins)/d_max;
            for (i = 0; i < num_puntos-1; i++)
            {
                for (j = i+1; j < num_puntos; j++)
                {
                    // Distancia euclidea
                    dd = dist(data[i].x-data[j].x, data[i].y - data[j].y, data[i].z - data[j].z);
                    rr = dist(rand[i].x-rand[j].x, rand[i].y - rand[j].y, rand[i].z - rand[j].z);
                    if (dd < d_max)
                    {
                        pos = (int)(dd*s);
                        DD[pos] += 2;
                    }
                    if (rr < d_max)
                    {
                        pos = (int)(rr*s);
                        RR[pos] += 2;
                    }
                }
            }

        }
        void histogramasMixtos(float *DR){
            int i,j,pos;
            float dr, s, aux;
            s = (float)(num_bins)/d_max;
            for (i = 0; i < num_puntos; i++)
            {
                for (j = 0; j < num_puntos; j++)
                {
                    // Método aprox
                    //aux = distancia(data[i].x - rand[j].x, data[i].y - rand[j].y);
                    //dr = distancia(aux, data[i].z - rand[j].z);

                    // Distancia euclidea
                    dr = dist(data[i].x - rand[j].x, data[i].y - rand[j].y, data[i].z - rand[j].z);
                    if (dr < d_max)
                    {
                        pos = (int)(dr*s);
                        DR[pos] += 1;
                    }
                }
            }
        }
        ~iso2PCF(){ // destructor

        }
};
