
#include "distancias.h"

class ani2PCF{
    //Atributos
    private:
        Puntos *data;
        Puntos *rand;
        int num_bins;
        int num_datos;
        float d_max;        
    public:
        ani2PCF(Puntos*,Puntos*,int,int,float);
        void setData(Puntos*);
        void setRand(Puntos*);
        void setNum_Bins(int);
        void setNum_Puntos(int);
        void setD_Max(float);
        void calcular_histogramas_puros(float**,float**);
        void calcular_histogramas_mixtos(float**);
        ~ani2PCF();
};
//Constructor
ani2PCF::ani2PCF(Puntos *_data, Puntos *_rand, int _num_datos, int _num_bins, float _d_max)
{
    data = _data;
    rand = _rand;
    num_datos = _num_datos;
    num_bins = _num_bins;
    d_max = _d_max;
}

void ani2PCF::setData(Puntos *_data){
    data = _data;
}

void ani2PCF::setRand(Puntos *_rand){
    rand = _rand;
}

void ani2PCF::setNum_Puntos(int _num_puntos){
    num_datos = _num_puntos;
}

void ani2PCF::setNum_Bins(int _num_bins){
    num_bins = _num_bins;
}

void ani2PCF::setD_Max(float _d_max){
    d_max = _d_max;
}

void ani2PCF::calcular_histogramas_puros(float** DD, float** RR){
    int i,j,a,b;
    float ds = (float)(num_bins)/d_max , r_paralelo, r_ortogonal;
    for (i = 0; i < num_datos; i++)
    {
        for (j = i+1; j < num_datos-1; j++)
        {
            r_paralelo = abs(data[i].z - data[j].z);
            if (r_paralelo < d_max)
            {
                r_ortogonal = euclidean_dist2D(data[i].x - data[j].x, data[i].y - data[j].y);
                if (r_ortogonal < d_max)
                {
                    a = (int)(r_paralelo*ds);
                    b = (int)(r_ortogonal*ds);
                    *(*(DD+a)+b) += 2; //DD[a][b]
                }
            }
            r_paralelo = abs(rand[i].z - rand[j].z);
            if (r_paralelo < d_max)
            {
                r_ortogonal = euclidean_dist2D(rand[i].x - rand[j].x, rand[i].y - rand[j].y);
                if (r_ortogonal < d_max)
                {
                    a = (int)(r_paralelo*ds);
                    b = (int)(r_ortogonal*ds);
                    *(*(RR+a)+b) +=2; // DD[a][b]
                }
            }   
        }   
    }
}

void ani2PCF::calcular_histogramas_mixtos(float **DR){
    int i,j,a,b;
    float ds = (float)(num_bins)/d_max, r_paralelo, r_ortogonal;
    for (i = 0; i < num_datos; i++)
    {
        for (j = 0; j < num_datos; j++)
        {
            r_paralelo = abs(data[i].z - rand[j].z);
            if (r_paralelo < d_max)
            {
                r_ortogonal = euclidean_dist2D(data[i].x - rand[j].x, data[i].y - rand[j].y);
                if (r_ortogonal < d_max)
                {
                    a = (int)(r_paralelo*ds);
                    b = (int)(r_ortogonal*ds);
                    *(*(DR+a)+b) += 1; //DR[a][b]    
                }
                
            }
        }
    }
    
}

//Destructor
ani2PCF::~ani2PCF()
{
}
