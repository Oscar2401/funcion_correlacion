
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <omp.h>

struct Point3D{
	float x;
	float y; 
	float z;
};
struct PointW3D{
	float x;
	float y; 
	float z;
	float w;
};
struct Node{
	Point3D nodepos;	// Coordinates of the node (position of the node).
	int len;		// Number of elements in the node.
	PointW3D *elements;	// Node elements.
};

//=================================================================== 
//======================== Class ==================================== 
//=================================================================== 

class NODE2P{
	// Class attributes:
	private:
		// Assigned
		int bn;
		int n_pts;
		float size_box;
		float size_node;
		float d_max;
		Node ***nodeD;
		PointW3D *dataD;
		Node ***nodeR;
		PointW3D *dataR;
		// Derivatives
		float ll;
		float dd_max;
		float corr;
		float front;
		float ds;
		float ddmax_nod;
		
	private: 
		void make_nodos(Node ***, PointW3D *);
		void add(PointW3D *&, int&, float, float, float, float);
	
	// Class Methods:
	public:
		// Class constructor:
		NODE2P(int _bn, int _n_pts, float _size_box, float _size_node, float _d_max, PointW3D *_dataD, Node ***_nodeD, PointW3D *_dataR, Node ***_nodeR){
			
			// Assigned
			bn = _bn;
			n_pts = _n_pts;
			size_box = _size_box;
			size_node = _size_node;
			d_max = _d_max;
			dataD = _dataD;
			nodeD = _nodeD;
			dataR = _dataR;
			nodeR = _nodeR;
			
			// Derivatives
			ll = size_box*size_box;
			dd_max = d_max*d_max;
			front = size_box - d_max;
			corr = size_node*sqrt(3);
			ds = ((float)(bn))/d_max;
			ddmax_nod = (d_max+corr)*(d_max+corr);
			
			make_nodos(nodeD,dataD);
			make_nodos(nodeR,dataR); 
			std::cout << "I finished building nodes ..." << std::endl;
		}
		
		Node ***meshData(){
			return nodeD;
		};
		Node ***meshRand(){
			return nodeR;
		};
		
		// Implementamos Método de mallas:
		void make_histoXX(double*, Node ***);
		void histo_front_XX(double *, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int);
		void make_histoXY(double*, Node ***, Node ***);
		void histo_front_XY(double *, Node ***, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int);
		~NODE2P();
};

//=================================================================== 
//==================== Funciones ==================================== 
//===================================================================  

void NODE2P::make_nodos(Node ***nod, PointW3D *dat){
	/*
	This function classifies the data in the nodes
	
	Args
	nod: Node 3D array where the data will be classified
	dat: array of PointW3D data to be classified and stored in the nodes
	*/
	int i, row, col, mom, partitions = (int)((size_box/size_node)+1);
	float p_med = size_node/2;
	
	// First allocate memory as an empty node:
	for (row=0; row<partitions; row++){
	for (col=0; col<partitions; col++){
	for (mom=0; mom<partitions; mom++){
		nod[row][col][mom].nodepos.z = ((float)(mom)*(size_node))+p_med;
		nod[row][col][mom].nodepos.y = ((float)(col)*(size_node))+p_med;
		nod[row][col][mom].nodepos.x = ((float)(row)*(size_node))+p_med;
		nod[row][col][mom].len = 0;
		nod[row][col][mom].elements = new PointW3D[0];
	}}}
	
	// Classificate the ith elment of the data into a node and add that point to the node with the add function:
	for (i=0; i<n_pts; ++i){
		row = (int)(dat[i].x/size_node);
        	col = (int)(dat[i].y/size_node);
        	mom = (int)(dat[i].z/size_node);
		add( nod[row][col][mom].elements, nod[row][col][mom].len, dat[i].x, dat[i].y, dat[i].z, dat[i].w);
	}
}

//=================================================================== 
void NODE2P::add(PointW3D *&array, int &lon, float _x, float _y, float _z, float _w){
	lon++;
	PointW3D *array_aux = new PointW3D[lon];
	for (int i=0; i<lon-1; i++){
		array_aux[i].x = array[i].x;
		array_aux[i].y = array[i].y;
		array_aux[i].z = array[i].z;
		array_aux[i].w = array[i].w;
	}
	delete[] array;
	array = array_aux;
	array[lon-1].x = _x;
	array[lon-1].y = _y; 
	array[lon-1].z = _z;
	array[lon-1].w = _w; 
}

//=================================================================== 
void NODE2P::make_histoXX(double *XX, Node ***nodeX){
	/*
	Function to create the DD and RR histograms.

	Arguments
	XX: array where the DD histogram will be created.
	nodeX: array of nodes.
	
	*/
	
	int partitions = (int)((size_box/size_node)+1);
	float d_max_pm = d_max + size_node/2, front_pm = front - size_node/2;
	
	#pragma omp parallel num_threads(2) 
	{
    	
    	double *SS;
    	SS = new double[bn];
    	for (int k = 0; k < bn; k++) *(SS+k) = 0;
    	
    	// Private variables in threads:
	int i, j, row, col, mom, u, v, w;
	float dis, dis_nod;
	float x1D, y1D, z1D, x2D, y2D, z2D;
	float x, y, z, w1;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod;
	bool con_x, con_y, con_z;
	
	#pragma omp for collapse(3)  schedule(dynamic)
	for (row = 0; row < partitions; ++row){
	for (col = 0; col < partitions; ++col){
	for (mom = 0; mom < partitions; ++mom){
	x1D = nodeX[row][0][0].nodepos.x;
	y1D = nodeX[row][col][0].nodepos.y;
	z1D = nodeX[row][col][mom].nodepos.z;			
		//==================================================
		// Pairs of points in the same node:
		//==================================================
		for ( i= 0; i <nodeX[row][col][mom].len - 1; ++i){
		x = nodeX[row][col][mom].elements[i].x;
		y = nodeX[row][col][mom].elements[i].y;
		z = nodeX[row][col][mom].elements[i].z;
		w1 = nodeX[row][col][mom].elements[i].w;
			for ( j = i+1; j < nodeX[row][col][mom].len; ++j){
			dx = x-nodeX[row][col][mom].elements[j].x;
			dy = y-nodeX[row][col][mom].elements[j].y;
			dz = z-nodeX[row][col][mom].elements[j].z;
			dis = dx*dx+dy*dy+dz*dz;
			if (dis <= dd_max){
			*(SS + (int)(sqrt(dis)*ds)) += 2*w1*nodeX[row][col][mom].elements[j].w;
			}
			}
		}
		//==================================================
		// Pairs of points at different nodes
		//==================================================
		u = row;
		v = col;
		//=========================
		// N2 mobile in Z
		//=========================
		for (w=mom+1;  w<partitions ; ++w){	
		z2D = nodeX[u][v][w].nodepos.z;
		dz_nod = z1D-z2D;
		dis_nod = dz_nod*dz_nod;
		if (dis_nod <= ddmax_nod){
			for ( i = 0; i < nodeX[row][col][mom].len; ++i){
			x = nodeX[row][col][mom].elements[i].x;
			y = nodeX[row][col][mom].elements[i].y;
			z = nodeX[row][col][mom].elements[i].z;
			w1 = nodeX[row][col][mom].elements[i].w;
				for ( j = 0; j < nodeX[u][v][w].len; ++j){
				dx = x-nodeX[u][v][w].elements[j].x;
				dy = y-nodeX[u][v][w].elements[j].y;
				dz = z-nodeX[u][v][w].elements[j].z;
				dis = dx*dx+dy*dy+dz*dz;
				if (dis <= dd_max){
				*(SS + (int)(sqrt(dis)*ds)) += 2*w1*nodeX[u][v][w].elements[j].w;
				}
				}
				}
			}
			// ======================================= 
			// Distance of border points XX 
			// ======================================= 
			// Boundary node conditions:
			con_z = ((z1D<=d_max_pm)&&(z2D>=front_pm))||((z2D<=d_max_pm)&&(z1D>=front_pm));
			if(con_z){
			histo_front_XX(SS,nodeX,dis_nod,0.0,0.0,fabs(dz_nod),false,false,con_z,row,col,mom,u,v,w);
			}
		}
		//=========================
		// N2 mobile in ZY
		//=========================
		for (v = col + 1; v < partitions ; ++v){
		y2D = nodeX[u][v][0].nodepos.y;
		dy_nod = y1D-y2D;
		dy_nod *= dy_nod;
			for (w = 0; w < partitions ; ++w){		
			z2D = nodeX[u][v][w].nodepos.z;
			dz_nod = z1D-z2D;
			dz_nod *= dz_nod;
			dis_nod = dy_nod + dz_nod;
			if (dis_nod <= ddmax_nod){
				for ( i = 0; i < nodeX[row][col][mom].len; ++i){
				x = nodeX[row][col][mom].elements[i].x;
				y = nodeX[row][col][mom].elements[i].y;
				z = nodeX[row][col][mom].elements[i].z;
				w1 = nodeX[row][col][mom].elements[i].w;
					for ( j = 0; j < nodeX[u][v][w].len; ++j){	
					dx =  x-nodeX[u][v][w].elements[j].x;
					dy =  y-nodeX[u][v][w].elements[j].y;
					dz =  z-nodeX[u][v][w].elements[j].z;
					dis = dx*dx+dy*dy+dz*dz;
					if (dis <= dd_max){
					*(SS + (int)(sqrt(dis)*ds)) += 2*w1*nodeX[u][v][w].elements[j].w;
					}
					}
				}
			}
			// ======================================= 
			// Distance of border points XX 
			// ======================================= 
			// Boundary node conditions:
			con_y = ((y1D<=d_max_pm)&&(y2D>=front_pm))||((y2D<=d_max_pm)&&(y1D>=front_pm));
			con_z = ((z1D<=d_max_pm)&&(z2D>=front_pm))||((z2D<=d_max_pm)&&(z1D>=front_pm));
			if(con_y || con_z){ 
			histo_front_XX(SS,nodeX,dis_nod,0.0,sqrt(dy_nod),sqrt(dz_nod),false,con_y,con_z,row,col,mom,u,v,w);
			}
			}
		}
		//=========================
		// N2 mobile in  ZYX
		//=========================
		for ( u = row + 1; u < partitions; ++u){
		x2D = nodeX[u][0][0].nodepos.x;
		dx_nod = x1D-x2D;
		dx_nod *= dx_nod;	
			for ( v = 0; v < partitions; ++v){
			y2D = nodeX[u][v][0].nodepos.y;
			dy_nod = y1D-y2D;
			dy_nod *= dy_nod;	
				for ( w = 0; w < partitions; ++w){
				z2D = nodeX[u][v][w].nodepos.z;
				dz_nod = z1D-z2D;
				dz_nod *= dz_nod;
				dis_nod = dx_nod + dy_nod + dz_nod;
				if (dis_nod <= ddmax_nod){
					for ( i = 0; i < nodeX[row][col][mom].len; ++i){
					x = nodeX[row][col][mom].elements[i].x;
					y = nodeX[row][col][mom].elements[i].y;
					z = nodeX[row][col][mom].elements[i].z;
					w1 = nodeX[row][col][mom].elements[i].w;
						for ( j = 0; j < nodeX[u][v][w].len; ++j){	
						dx = x-nodeX[u][v][w].elements[j].x;
						dy = y-nodeX[u][v][w].elements[j].y;
						dz = z-nodeX[u][v][w].elements[j].z;
						dis = dx*dx + dy*dy + dz*dz;
						if (dis <= dd_max){
							*(SS + (int)(sqrt(dis)*ds)) += 2*w1*nodeX[u][v][w].elements[j].w;
						}
						}
					}
					}
				// ======================================= 
				// Distance of border points XX 
				// ======================================= 
				// Boundary node conditions:
				con_x = ((x1D<=d_max_pm)&&(x2D>=front_pm))||((x2D<=d_max_pm)&&(x1D>=front_pm));
				con_y = ((y1D<=d_max_pm)&&(y2D>=front_pm))||((y2D<=d_max_pm)&&(y1D>=front_pm));
				con_z = ((z1D<=d_max_pm)&&(z2D>=front_pm))||((z2D<=d_max_pm)&&(z1D>=front_pm));
				if(con_x || con_y || con_z){
				histo_front_XX(SS,nodeX,dis_nod,sqrt(dx_nod),sqrt(dy_nod),sqrt(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w);
				}	
				}	
			}
		}
	}
	}
	}
	#pragma omp critical
	for(int a=0; a<bn; a++) *(XX+a)+=*(SS+a);
	}
}

//=================================================================== 
void NODE2P::histo_front_XX(double *PP, Node ***dat, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int row, int col, int mom, int u, int v, int w){
	/*
	
	Function to add periodic boundary conditions
	
	*/
	
	int i, j;
	float dis_f,dis,d_x,d_y,d_z;
	float x,y,z,w1;
	//======================================================================
	if( con_in_x ){
	dis_f = disn + ll - 2*dn_x*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = y-dat[u][v][w].elements[j].y;
			d_z = z-dat[u][v][w].elements[j].z;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += 2*w1*dat[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================		
	if( con_in_y ){
	dis_f = disn + ll - 2*dn_y*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = x-dat[u][v][w].elements[j].x;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = z-dat[u][v][w].elements[j].z;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += 2*w1*dat[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================
	if( con_in_z ){
	dis_f = disn + ll - 2*dn_z*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = x-dat[u][v][w].elements[j].x;
			d_y = y-dat[u][v][w].elements[j].y;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += 2*w1*dat[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================	
	if( con_in_x && con_in_y ){
	dis_f = disn + 2*ll - 2*(dn_x+dn_y)*size_box;
	if (dis_f < ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = z-dat[u][v][w].elements[j].z;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += 2*w1*dat[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================			
	if( con_in_x && con_in_z ){
	dis_f = disn + 2*ll - 2*(dn_x+dn_z)*size_box;
	if (dis_f <= ddmax_nod){
		for ( i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = y-dat[u][v][w].elements[j].y;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += 2*w1*dat[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================		
	if( con_in_y && con_in_z ){
	dis_f = disn + 2*ll - 2*(dn_y+dn_z)*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = x-dat[u][v][w].elements[j].x;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += 2*w1*dat[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================		
	if( con_in_x && con_in_y && con_in_z ){
	dis_f = disn + 3*ll - 2*(dn_x+dn_y+dn_z)*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = d_x*d_x + d_y*d_y + d_z*d_z;
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += 2*w1*dat[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
}
//=================================================================== 
void NODE2P::make_histoXY(double *XY, Node ***nodeX, Node ***nodeY){
	/*
	Function to create the DD and RR histograms.

	Arguments
	XX: array where the DD histogram will be created.
	nodeX: array of nodes.
	
	*/
	
	int partitions = (int)((size_box/size_node)+1);
	float d_max_pm = d_max + size_node/2, front_pm = front - size_node/2;
	
	#pragma omp parallel num_threads(2) 
	{
    	
    	double *SS;
    	SS = new double[bn];
    	for (int k = 0; k < bn; k++) *(SS+k) = 0.0;
    	
    	// Private variables in threads:
	int i, j, row, col, mom, u, v, w;
	float dis, dis_nod;
	float x1D, y1D, z1D, x2R, y2R, z2R;
	float x, y, z, w1;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod;
	bool con_x, con_y, con_z;
	
	#pragma omp for collapse(3)  schedule(dynamic)
	for (row = 0; row < partitions; ++row){
	for (col = 0; col < partitions; ++col){
	for (mom = 0; mom < partitions; ++mom){
	x1D = nodeX[row][0][0].nodepos.x;
	y1D = nodeX[row][col][0].nodepos.y;
	z1D = nodeX[row][col][mom].nodepos.z;	
		//=========================
		// N2 mobile in ZYX
		//=========================
		for (u=0; u<partitions; ++u){
		x2R = nodeY[u][0][0].nodepos.x;
		dx_nod = x1D-x2R;
		dx_nod *= dx_nod;
			for (v=0; v<partitions; ++v){
			y2R = nodeY[u][v][0].nodepos.y;
			dy_nod = y1D-y2R;
			dy_nod *= dy_nod;
				for (w=0; w<partitions; ++w){
				z2R = nodeY[u][v][w].nodepos.z;
				dz_nod = z1D-z2R;
				dz_nod *= dz_nod; 
				dis_nod = dx_nod + dy_nod + dz_nod;
				if (dis_nod <= ddmax_nod){
					for (i=0; i<nodeX[row][col][mom].len; ++i){
					x = nodeX[row][col][mom].elements[i].x;
					y = nodeX[row][col][mom].elements[i].y;
					z = nodeX[row][col][mom].elements[i].z;
					w1 = nodeX[row][col][mom].elements[i].w;
						for (j=0; j<nodeY[u][v][w].len; ++j){	
						dx = x-nodeY[u][v][w].elements[j].x;
						dy = y-nodeY[u][v][w].elements[j].y;
						dz = z-nodeY[u][v][w].elements[j].z;
						dis = dx*dx + dy*dy + dz*dz;
						if (dis < dd_max){
						*(SS + (int)(sqrt(dis)*ds)) += w1*nodeY[u][v][w].elements[j].w;
						}
						}
					}
				}
				// ======================================= 
				// Distance of border points XX 
				// ======================================= 
				// Boundary node conditions:
				con_x = ((x1D<=d_max_pm)&&(x2R>=front_pm))||((x2R<=d_max_pm)&&(x1D>=front_pm));
				con_y = ((y1D<=d_max_pm)&&(y2R>=front_pm))||((y2R<=d_max_pm)&&(y1D>=front_pm));
				con_z = ((z1D<=d_max_pm)&&(z2R>=front_pm))||((z2R<=d_max_pm)&&(z1D>=front_pm));
				if(con_x || con_y || con_z){
				histo_front_XY(SS,nodeX,nodeY,dis_nod,sqrt(dx_nod),sqrt(dy_nod),sqrt(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w);
				}	
				}	
			}
		}
	}
	}
	}
	#pragma omp critical
	for(int a=0; a<bn; a++) *(XY+a)+=*(SS+a);
	}
}
//=================================================================== 
void NODE2P::histo_front_XY(double *PP, Node ***dat, Node ***ran, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int row, int col, int mom, int u, int v, int w){
	/*
	
	Function to add periodic boundary conditions
	
	*/
	
	int i, j;
	float dis_f,dis,d_x,d_y,d_z;
	float x,y,z,w1;
	//======================================================================
	if( con_in_x ){
	dis_f = disn + ll - 2*dn_x*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = fabs(x-ran[u][v][w].elements[j].x)-size_box;
			d_y = y-ran[u][v][w].elements[j].y;
			d_z = z-ran[u][v][w].elements[j].z;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================		
	if( con_in_y ){
	dis_f = disn + ll - 2*dn_y*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = x-ran[u][v][w].elements[j].x;
			d_y = fabs(y-ran[u][v][w].elements[j].y)-size_box;
			d_z = z-ran[u][v][w].elements[j].z;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================
	if( con_in_z ){
	dis_f = disn + ll - 2*dn_z*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = x-ran[u][v][w].elements[j].x;
			d_y = y-ran[u][v][w].elements[j].y;
			d_z = fabs(z-ran[u][v][w].elements[j].z)-size_box;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================	
	if( con_in_x && con_in_y ){
	dis_f = disn + 2*ll - 2*(dn_x+dn_y)*size_box;
	if (dis_f < ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = fabs(x-ran[u][v][w].elements[j].x)-size_box;
			d_y = fabs(y-ran[u][v][w].elements[j].y)-size_box;
			d_z = z-ran[u][v][w].elements[j].z;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================			
	if( con_in_x && con_in_z ){
	dis_f = disn + 2*ll - 2*(dn_x+dn_z)*size_box;
	if (dis_f <= ddmax_nod){
		for ( i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = fabs(x-ran[u][v][w].elements[j].x)-size_box;
			d_y = y-ran[u][v][w].elements[j].y;
			d_z = fabs(z-ran[u][v][w].elements[j].z)-size_box;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================		
	if( con_in_y && con_in_z ){
	dis_f = disn + 2*ll - 2*(dn_y+dn_z)*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = x-ran[u][v][w].elements[j].x;
			d_y = fabs(y-ran[u][v][w].elements[j].y)-size_box;
			d_z = fabs(z-ran[u][v][w].elements[j].z)-size_box;
			dis = d_x*d_x + d_y*d_y + d_z*d_z; 
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
	//======================================================================		
	if( con_in_x && con_in_y && con_in_z ){
	dis_f = disn + 3*ll - 2*(dn_x+dn_y+dn_z)*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		w1 = dat[row][col][mom].elements[i].w;
			for (j=0; j<ran[u][v][w].len; ++j){
			d_x = fabs(x-ran[u][v][w].elements[j].x)-size_box;
			d_y = fabs(y-ran[u][v][w].elements[j].y)-size_box;
			d_z = fabs(z-ran[u][v][w].elements[j].z)-size_box;
			dis = d_x*d_x + d_y*d_y + d_z*d_z;
			if (dis < dd_max){
			*(PP + (int)(sqrt(dis)*ds)) += w1*ran[u][v][w].elements[j].w;
			}
			}
		}
	}
	}
}
NODE2P::~NODE2P(){
	
}
