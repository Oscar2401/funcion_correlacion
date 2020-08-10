#include "distancias.h"

class iso3PCF{
    //Atributos
    private:
        Puntos *DATA;
        Puntos *RAND;
        int num_puntos;
        int num_bins;
        float d_max;
    //Métodos
    private:
        void XXX(int, int, Puntos*, float ***);
        void XXY(int,int, Puntos*, Puntos*,float ***);
    public:
        void setData(Puntos*);
        void setRand(Puntos*);
        void setNum_Bins(int);
        void setNum_Puntos(int);
        void setD_Max(float);
        void calcular_Histogramas_Puros(float ***, float ***);
        void calcular_Histogramas_Mixtos(float ***, float ***);
        void simetrizar_Histograma(float ***);
        iso3PCF(Puntos*, Puntos*, int, int, float);
        ~iso3PCF();
};

iso3PCF::iso3PCF(Puntos *_data, Puntos *_rand, int _num_datos, int _num_bins, float _d_max){
    DATA = _data;
    RAND = _rand;
    num_puntos = _num_datos;
    num_bins = _num_bins;
    d_max = _d_max;
}

iso3PCF::~iso3PCF(){
}

void iso3PCF::XXX(int i,int j,Puntos *p, float ***HHH){
    int k,a,b,c;
    float l1,l2,l3, ds = (float)(num_bins)/d_max;
    l1 = euclidean_dist3D(p[i].x - p[j].x, p[i].y - p[j].y, p[i].z - p[j].z);
    if (l1 < d_max)
    {
        for (k = j+1; k < num_puntos; k++)
        {
            l2 = euclidean_dist3D(p[i].x - p[k].x, p[i].y - p[k].y, p[i].z - p[k].z);
            if (l2 < d_max)
            {
                l3 =  euclidean_dist3D(p[j].x - p[k].x, p[j].y - p[k].y, p[j].z - p[k].z);
                if (l3 < d_max)
                {
                    a = (int)(l1*ds);
                    b = (int)(l2*ds);
                    c = (int)(l3*ds);
                    *(*(*(HHH + a) + b) +c) += 1;                        
                }
            }
        }
    }
}

void iso3PCF::XXY(int i, int j, Puntos *p, Puntos *q, float ***HHH){
    int k,a,b,c;
    float ds = (float)(num_bins)/d_max,l1,l2,l3, inc = 1.0/3.0;
    l1 = euclidean_dist3D(p[i].x - p[j].x, p[i].y - p[j].y, p[i].z - p[j].z);
    if (l1 < d_max)
    {
        for (k = 0; k < num_puntos; k++)
        {
            l2 = euclidean_dist3D(p[i].x-q[k].x, p[i].y - q[k].y, p[i].z - q[k].z);
            if (l2 < d_max)
            {
                l3 = euclidean_dist3D(p[j].x - q[k].x, p[j].y - q[k].y, p[j].z - q[k].z);
                if (l3 < d_max)
                {
                    a = (int)(l1*ds);
                    b = (int)(l2*ds);
                    c = (int)(l3*ds);
                    *(*(*(HHH + a)+b)+c) += inc;
                }   
            }
        }
    }
}

void iso3PCF::setData(Puntos *_data){
    DATA = _data;
}

void iso3PCF::setRand(Puntos *_rand){
    RAND = _rand;
}

void iso3PCF::setNum_Puntos(int _num_puntos){
    num_puntos = _num_puntos;
}

void iso3PCF::setNum_Bins(int _num_bins){
    num_bins = _num_bins;
}

void iso3PCF::setD_Max(float _d_max){
    d_max = _d_max;
}

void iso3PCF::calcular_Histogramas_Puros(float ***DDD, float ***RRR){
    int i,j;
    for (i = 0; i < num_puntos-2; i++)
    {
        for (j = i + 1; j < num_puntos-1; j++)
        {
            // para DDD
            XXX(i,j,DATA,DDD);
            // para RRR
            XXX(i,j,RAND,RRR);
        }
    }
}

void iso3PCF::simetrizar_Histograma(float ***DDD){
    int i,j,k;
    float valor;
    for (i = 0; i < num_bins-2; i++)
    {
        for (j = i + 1; j < num_bins-1; j++)
        {
            for (k = j + 1; k < num_bins; k++)
            {
                valor = DDD[i][j][k] + DDD[i][k][j] + DDD[j][k][i] + DDD[j][i][k] + DDD[k][i][j] + DDD[k][j][i];
                DDD[i][k][j] = valor;
                DDD[j][k][i] = valor;
                DDD[j][i][k] = valor;
                DDD[k][i][j] = valor;
                DDD[k][j][i] = valor;
            }   
        }
    }
}

void iso3PCF::calcular_Histogramas_Mixtos(float ***DDR, float ***DRR){
    int i,j;
    for (i = 0; i < num_puntos -1 ; i++)
    {
        for (j = i+1; j < num_puntos; j++)
        {
            //Para DDR
            XXY(i,j,DATA, RAND,DDR);
            //para RRD
            XXY(i,j,RAND,DATA,DRR);
        }
    }
}
