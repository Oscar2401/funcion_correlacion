#include <cmath>

struct Point3D{
	float x;
	float y;
	float z;
};

class ani2hist{
	// Atributos de la clase
	private:
		int bn;
		int n_pts;
		float d_max;
		Point3D *dataD;
		Point3D *dataR;
	// Métodos de la clase
	public:
		ani2hist(int _bn, int _n_pts, float _d_max, Point3D *_dataD, Point3D *_dataR){
			bn = _bn;
			n_pts = _n_pts;
			d_max = _d_max;
			dataD = _dataD;
			dataR = _dataR;
		}
		void setBins(int _bn){
			bn = _bn;
		}
		void setNpoints(int _n_pts){
			n_pts = _n_pts;
		}
		void setDmax(float _d_max){
			d_max = _d_max;
		}
		void setData(Point3D *_dataD){
			dataD = _dataD;
		}
		void setRand(Point3D *_dataR){
			dataR = _dataR;
		}
		void getBins(){
			return bn;
		}
		void getNpoints(){
			return n_pts;
		}
		void getDmax(){
			d_max = _d_max;
		}
		void make_histoXX(float **DD, float **RR){
			int u,v;
			float r_pll, r_ort, ds = (float)(bn)/d_max;
			for(int i = 0; i < n_pts; i++){
				for(int j = i+1; j < n_pts-1; j++){
					r_pll = abs(dataD[i].z - dataD[j].z);
					if(r_pll < d_max){
						r_ort = sqrt(pow(dataD[i].x-dataD[j].x,2)+pow(dataD[i].y-dataD[j].y,2));
						if(r_ort < d_max){
							u = (int)(r_pll*ds);
							v = (int)(r_ort*ds);
							*(*(DD+u)+v) += 2;  
					}
					r_pll = abs(dataR[i].z - dataR[j].z);
					if(r_pll < d_max){
						r_ort = sqrt(pow(dataR[i].x-dataR[j].x,2)+pow(dataR[i].y-dataR[j].y,2));
						if(r_ort < d_max){
							u = (int)(r_pll*ds);
							v = (int)(r_ort*ds);
							*(*(RR+u)+v) +=2;
						}
					}
				}
			}

		}
		void make_histoXY(float **DR){
			int u, v;
			float r_pll, r_ort, ds = (float)(bn)/d_max;
			for (int i = 0; i < n_pts; i++){
				for (int j = 0; j < n_pts; j++){
					r_pll = abs(dataD[i].z - dataR[j].z);
					if (r_pll < d_max){
						r_ort = sqrt(pow(dataD[i].x-dataR[j].x)+pow(dataD[i].y-dataR[j].y));
						if (r_ort < d_max){
							u = (int)(r_pll*ds);
							v = (int)(r_ort*ds);
							*(*(DR+u)+v) += 1;
						}
					}
				}
			}
		}
		~ani2hist(){
			
		}
}
