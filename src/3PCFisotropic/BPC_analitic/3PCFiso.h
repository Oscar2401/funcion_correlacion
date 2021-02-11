#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <omp.h>

struct Point3D{
	float x;
	float y; 
	float z;
	short int fx;
	short int fy;
	short int fz;
};

struct PointW3D{
	float x;
	float y; 
	float z;
	float w;
	short int fx;
	short int fy;
	short int fz;
};

struct Node{
	Point3D nodepos;	// Coordinates of the node (position of the node).
	int len;		// Number of points in the node
	PointW3D *elements;	// Points in the node
};

//=================================================================== 
//======================== Clase ==================================== 
//=================================================================== 


class NODE3P{
	// Class attributes:
	private:
		// Assigned
		int bn;
		int n_pts;
		int partitions;
		float size_box;
		float size_node;
		float d_max;
		Node ***nodeD;
		PointW3D *dataD;
		// Derivatives
		float ll;
		float dd_max;
		float corr;
		float front;
		double ds;
		float ddmax_nod;
		float d_max_pm;
		float front_pm;
		float size_box_2;
		
	private: 
		void make_nodos(Node ***, PointW3D *);
		void add(PointW3D *&, int&, float, float, float, float, short int, short int, short int);
	
	// Class methods:
	public:
		// Class constructor:
		NODE3P(int _bn, int _n_pts, float _size_box, float _size_node, float _d_max, PointW3D *_dataD, Node ***_nodeD){
			// Assigned
			bn = _bn;
			n_pts = _n_pts;
			size_box = _size_box;
			size_node = _size_node;
			d_max = _d_max;
			dataD = _dataD;
			nodeD = _nodeD;
			// Derivatives
			dd_max = d_max*d_max;
			front = size_box - d_max;
			corr = size_node*sqrt(3);
			ds = floor(((double)(bn)/d_max)*1000000)/1000000;
			ddmax_nod = (d_max+corr)*(d_max+corr);
			partitions = (int)(ceil(size_box/size_node));
			d_max_pm = d_max + size_node/2;
			front_pm = front - size_node/2;
			size_box_2 =  size_box/2;
			
			make_nodos(nodeD,dataD);
		}
		
		Node ***meshData(){
			return nodeD;
		};
		
		// Implementing grid method:
		void make_histoXXX(double ***, Node ***);
		void count_3_N111(int, int, int, double ***, Node ***);
		void count_3_N112(int, int, int, int, int, int, double ***, Node ***);
		void count_3_N123(int, int, int, int, int, int, int, int, int, double ***, Node ***);
		void symmetrize(double ***);
		void symmetrize_analitic(double ***);
		void BPC(double ***, PointW3D *);
		
		void front_node_112(int, int, int, int, int, int, bool, bool, bool, double ***, Node ***);
		void front_112(int,int,int,int,int,int,short int,short int,short int, double ***,Node ***);
		void front_node_123(int,int,int,int,int,int,int,int,int,bool,bool,bool, double ***,Node ***);
		void front_123(int,int,int,int,int,int,int,int,int,short int,short int,short int,short int,short int,short int,short int,short int,short int,double ***, Node ***);
		
		void make_histo_analitic(double ***, double ***, Node ***);
		void make_histoXX(double *, double *, Node ***, int);
		void histo_front_XX(double *, Node ***, float, float, float, float, bool, bool, bool, int, int, int, int, int, int, double);
		
		~NODE3P();
};

//=================================================================== 
//==================== Funciones ==================================== 
//=================================================================== 
void NODE3P::make_nodos(Node ***nod, PointW3D *dat){
	/*
	This function classifies the data in the nodes
	
	Args
	nod: Node 3D array where the data will be classified
	dat: array of PointW3D data to be classified and stored in the nodes
	*/
	
	float p_med = size_node/2;
	#pragma omp parallel num_threads(50) 
    	{
	int i, row, col, mom;
	float posx, posy, posz;
	#pragma omp for collapse(3)  schedule(dynamic)
	// First allocate memory as an empty node:
	for (row=0; row<partitions; ++row){
	for (col=0; col<partitions; ++col){
	for (mom=0; mom<partitions; ++mom){
		posx = ((float)(row)*(size_node))+p_med;
		posy = ((float)(col)*(size_node))+p_med;
		posz = ((float)(mom)*(size_node))+p_med;
		
		nod[row][col][mom].nodepos.x = posx;
		nod[row][col][mom].nodepos.y = posy;
		nod[row][col][mom].nodepos.z = posz;
		
		// We see if the node is on the border
		// border x:
		if(posx<=d_max_pm) nod[row][col][mom].nodepos.fx = -1;
		else if(posx>=front_pm) nod[row][col][mom].nodepos.fx = 1;
		else nod[row][col][mom].nodepos.fx = 0;
		// border y:
		if(posy<=d_max_pm) nod[row][col][mom].nodepos.fy = -1;
		else if(posy>=front_pm) nod[row][col][mom].nodepos.fy = 1;
		else nod[row][col][mom].nodepos.fy = 0;
		// border z:
		if(posz<=d_max_pm) nod[row][col][mom].nodepos.fz = -1;
		else if(posz>=front_pm) nod[row][col][mom].nodepos.fz = 1;
		else nod[row][col][mom].nodepos.fz = 0;
		
		nod[row][col][mom].len = 0;
		nod[row][col][mom].elements = new PointW3D[0];
	}
	}
	}
	}
	// Fill the nodes with the dat points:
	int i, row, col, mom;
	for (int i=0; i<n_pts; ++i){
		row = (int)(dat[i].x/size_node);
        	col = (int)(dat[i].y/size_node);
        	mom = (int)(dat[i].z/size_node);
		add(nod[row][col][mom].elements,nod[row][col][mom].len,dat[i].x,dat[i].y,dat[i].z,dat[i].w,dat[i].fx,dat[i].fy,dat[i].fz);
	}
}
//=================================================================== 
void NODE3P::add(PointW3D *&array, int &lon, float _x, float _y, float _z, float _w, short int _fx, short int _fy, short int _fz){
	lon++;
	PointW3D *array_aux = new PointW3D[lon];
	for (int i=0; i<lon-1; ++i){
		array_aux[i].x = array[i].x;
		array_aux[i].y = array[i].y;
		array_aux[i].z = array[i].z;
		array_aux[i].w = array[i].w;
		array_aux[i].fx = array[i].fx;
		array_aux[i].fy = array[i].fy;
		array_aux[i].fz = array[i].fz;
	}
	delete[] array;
	array = array_aux;
	array[lon-1].x = _x;
	array[lon-1].y = _y; 
	array[lon-1].z = _z;
	array[lon-1].w = _w; 
	array[lon-1].fx = _fx;
	array[lon-1].fy = _fy; 
	array[lon-1].fz = _fz;
}
//=================================================================== 
void NODE3P::make_histoXXX(double ***XXX, Node ***nodeX){
	/*
	Function to create the DDD and RRR histograms.
	
	Arg
	XXX: arrangement where the DDD y RRR histogram will be created.
	nodeX: array of nodes.

	*/
	std::cout << "Constructing DDD histogram..." << std::endl;
	#pragma omp parallel num_threads(50) 
	{
	
	int i, j, k, row, col, mom, u, v, w, a ,b, c;
	
	double ***SSS;
	
	SSS = new double**[bn];
	for (i=0; i<bn; i++){
	*(SSS+i) = new double*[bn];
		for (j = 0; j < bn; j++){
		*(*(SSS+i)+j) = new double[bn];
		}
	}
	
	for (i=0; i<bn; ++i){
	for (j=0; j<bn; ++j){
	for (k=0; k<bn; ++k) *(*(*(SSS+i)+j)+k) = 0.0;
	}
	}
	 
	double dis;
	float dis_nod, dis_nod2, dis_nod3;
	float x1N, y1N, z1N, x2N, y2N, z2N, x3N, y3N, z3N;
	float x, y, z;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod, dx_nod2, dy_nod2, dz_nod2, dx_nod3, dy_nod3, dz_nod3;
	float fx1N, fy1N, fz1N;
	bool con1_x, con1_y, con1_z;
	
	float fx2N, fy2N, fz2N;
	float fx3N, fy3N, fz3N;
	
	bool con2_x, con2_y, con2_z;
	bool con3_x, con3_y, con3_z;
	bool con12_x, con12_y, con12_z;
	bool con13_x, con13_y, con13_z;
	bool conx, cony, conz, conx_, cony_, conz_;
	
	#pragma omp for collapse(3)  schedule(dynamic)
	for (row=0; row<partitions; ++row){
	for (col=0; col<partitions; ++col){
	for (mom=0; mom<partitions; ++mom){
	x1N = nodeX[row][0][0].nodepos.x;
	fx1N = nodeX[row][0][0].nodepos.fx;
	con1_x = fx1N != 0;
	
	y1N = nodeX[row][col][0].nodepos.y;
	fy1N = nodeX[row][col][0].nodepos.fy;
	con1_y = fy1N != 0;
	
	z1N = nodeX[row][col][mom].nodepos.z;
	fz1N = nodeX[row][col][mom].nodepos.fz;
	con1_z = fz1N != 0;
	
	//==================================================
	// Triangles between points of the same node:
	//==================================================
	count_3_N111(row, col, mom, SSS, nodeX);		
	//==================================================
	// Triangles between points of the different node:
	//==================================================
	u = row;
	v = col;
	//=======================
	// Mobile node2 in Z:
	//=======================
	x2N = nodeX[u][0][0].nodepos.x;
	y2N = nodeX[u][v][0].nodepos.y;
	
	for (w=mom+1;  w<partitions; ++w){	
	z2N = nodeX[u][v][w].nodepos.z;
	dz_nod = z2N-z1N;
	dz_nod *= dz_nod;
	
	// inside box
	
	if (dz_nod < ddmax_nod){
	//==============================================
	// 2 points in node1 and 1 point in node2
	//==============================================
	count_3_N112(row, col, mom, u, v, w, SSS, nodeX);
	//==============================================
	// 1 point in node1, 1 point in node2 and 1 point in node3
	//==============================================
	a = u;
	b = v;
		x3N = nodeX[a][0][0].nodepos.x;
		y3N = nodeX[a][b][0].nodepos.y;
		//=======================
		// Mobile node3 in Z:
		//=======================
		for (c=w+1; c<partitions; ++c){
		z3N = nodeX[a][b][c].nodepos.z; 
		dz_nod2 = z3N-z1N;
		dz_nod2 *= dz_nod2;
		if (dz_nod2 < ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, SSS, nodeX);
		}
		//=======================
		// Mobile node3 in ZY:
		//=======================
		for (b=v+1; b<partitions; ++b){
		y3N = nodeX[a][b][0].nodepos.y;
		dy_nod2 = y3N-y1N;
		dy_nod2 *= dy_nod2;
			for (c=0; c<partitions; ++c){
			z3N = nodeX[a][b][c].nodepos.z;
			dz_nod2 = z3N-z1N;
			dz_nod2 *= dz_nod2;
			dis_nod2 = dy_nod2 + dz_nod2;
			if (dis_nod2 < ddmax_nod){
			dy_nod3 = y3N-y2N;
			dz_nod3 = z3N-z2N;
			dis_nod3 = dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
			if (dis_nod3 < ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, SSS, nodeX);
			}
			}
		}
		//=======================
		// Mobile node3 in ZYX:
		//=======================
		for (a=u+1; a<partitions; ++a){
		x3N = nodeX[a][0][0].nodepos.x;
		dx_nod2 = x3N-x1N;
		dx_nod2 *= dx_nod2;
			for (b=0; b<partitions; ++b){
			y3N = nodeX[a][b][0].nodepos.y;
			dy_nod2 = y3N-y1N;
			dy_nod2 *= dy_nod2;
				for (c=0; c<partitions; ++c){
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod2 = z3N-z1N;
				dz_nod2 *= dz_nod2;
				dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
				if (dis_nod2 < ddmax_nod){
				dx_nod3 = x3N-x2N;
				dy_nod3 = y3N-y2N;
				dz_nod3 = z3N-z2N;
				dis_nod3 = dx_nod3*dx_nod3 + dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
				if (dis_nod3 < ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, SSS, nodeX);
				}
				}
			}
		}
	}
	}
	//=======================
	// Mobile node2 in ZY:
	//=======================
	for (v=col+1; v<partitions ; ++v){
	y2N = nodeX[u][v][0].nodepos.y;
	dy_nod = y2N-y1N;
	dy_nod *= dy_nod;
		for (w=0; w<partitions ; ++w){		
		z2N = nodeX[u][v][w].nodepos.z;
		dz_nod = z2N-z1N;
		dz_nod *= dz_nod;
		dis_nod = dy_nod + dz_nod;
		if (dis_nod < ddmax_nod){
		//==============================================
		// 2 points in node1 and 1 point in node2
		//==============================================
		count_3_N112(row, col, mom, u, v, w, SSS, nodeX);
		//==============================================
		// 1 point in node1, 1 point in node2 and 1 point in node3
		//==============================================
		a = u;
		b = v;
		//=======================
		// Mobile node3 in Z:
		//=======================
		y3N = nodeX[a][b][0].nodepos.y;
		dy_nod2 = y3N-y1N;
		dy_nod2 *= dy_nod2;
			for (c=w+1;  c<partitions; ++c){
			z3N = nodeX[a][b][c].nodepos.z;
			dz_nod2 = z3N-z1N;
			dz_nod2 *= dz_nod2;
			dis_nod2 = dy_nod2 + dz_nod2;
			if (dis_nod2 < ddmax_nod){
			dz_nod3 = z3N-z2N;
			dis_nod3 = dz_nod3*dz_nod3;
			if (dis_nod3 < ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, SSS, nodeX);
			}
			}
			//=======================
			// Mobile node3 in ZY:
			//=======================	
			for (b=v+1; b<partitions; ++b){
			y3N = nodeX[a][b][0].nodepos.y;
			dy_nod2 = y3N-y1N;
			dy_nod2 *= dy_nod2;
				for (c=0;  c<partitions; ++c){
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod2 = z3N-z1N;
				dz_nod2 *= dz_nod2;
				dis_nod2 = dy_nod2 + dz_nod2;
				if (dis_nod2 < ddmax_nod){
				dy_nod3 = y3N-y2N;
				dz_nod3 = z3N-z2N;
				dis_nod3 = dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
				if (dis_nod3 < ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, SSS, nodeX);
				}
				}
			}
			//=======================
			// Mobile node3 in ZYX:
			//=======================
			for (a=u+1; a<partitions; ++a){
			x3N = nodeX[a][0][0].nodepos.x;
			dx_nod2 = x3N-x1N;
			dx_nod2 *= dx_nod2;
				for (b=0; b<partitions; ++b){
				y3N = nodeX[a][b][0].nodepos.y;
				dy_nod2 = y3N-y1N;
				dy_nod2 *= dy_nod2;
					for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					dz_nod2 *= dz_nod2;
					dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
					if (dis_nod2 < ddmax_nod){
					dx_nod3 = x3N-x2N;
					dy_nod3 = y3N-y2N;
					dz_nod3 = z3N-z2N;
					dis_nod3 = dx_nod3*dx_nod3 + dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
					if (dis_nod3 < ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, SSS, nodeX);
					}
					}
				}
			}
			}
		}	
	}			
	//=======================
	// Mobile node2 in ZYX:
	//=======================
	for (u=row+1; u<partitions; ++u){
	x2N = nodeX[u][0][0].nodepos.x;
	dx_nod = x2N-x1N;
	dx_nod *= dx_nod;
		for (v=0; v<partitions; ++v){
		y2N = nodeX[u][v][0].nodepos.y;
		dy_nod = y2N-y1N;
		dy_nod *= dy_nod;
			for (w=0; w<partitions; ++w){
			z2N = nodeX[u][v][w].nodepos.z;
			dz_nod = z2N-z1N;
			dz_nod *= dz_nod;
			dis_nod = dx_nod + dy_nod + dz_nod;
			if (dis_nod < ddmax_nod){
			//==============================================
			// 2 points in node1 and 1 point in node2
			//==============================================
			count_3_N112(row, col, mom, u, v, w, SSS, nodeX);
			//==============================================
			// 1 point in node1, 1 point in node2 and 1 point in node3
			//==============================================
			a = u;
			b = v;
			//=======================
			// Mobile node3 in Z:
			//=======================
			x3N = nodeX[a][0][0].nodepos.x;
			y3N = nodeX[a][b][0].nodepos.y;
			dx_nod2 = x3N-x1N;
			dx_nod2 *= dx_nod2; 
			dy_nod2 = y3N-y1N;
			dy_nod2 *= dy_nod2;
				for (c=w+1;  c<partitions; ++c){	
				z3N = nodeX[a][b][c].nodepos.z;
				dz_nod2 = z3N-z1N;
				dz_nod2 *= dz_nod2;
				dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
				if (dis_nod2 < ddmax_nod){
				dz_nod3 = z3N-z2N;
				dis_nod3 = dz_nod3*dz_nod3;
				if (dis_nod3 < ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, SSS, nodeX);
				}
				}
				//=======================
				// Mobile node3 in ZY:
				//=======================
				for (b=v+1; b<partitions; ++b){
				y3N = nodeX[a][b][0].nodepos.y;
				dy_nod2 = y3N-y1N;
				dy_nod2 *= dy_nod2;
					for (c=0;  c<partitions; ++c){
					z3N = nodeX[a][b][c].nodepos.z;
					dz_nod2 = z3N-z1N;
					dz_nod2 *= dz_nod2;
					dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
					if (dis_nod2 < ddmax_nod){
					dy_nod3 = y3N-y2N;
					dz_nod3 = z3N-z2N;
					dis_nod3 = dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
					if (dis_nod3 < ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, SSS, nodeX);
					}
					}
				}
				//=======================
				// Mobile node3 in ZYX:
				//=======================		
				for (a=u+1; a<partitions; ++a){
				x3N = nodeX[a][0][0].nodepos.x;
				dx_nod2 = x3N-x1N;
				dx_nod2 *= dx_nod2;
					for (b=0; b<partitions; ++b){
					y3N = nodeX[a][b][0].nodepos.y;
					dy_nod2 = y3N-y1N;
					dy_nod2 *= dy_nod2;
						for (c=0;  c<partitions; ++c){
						z3N = nodeX[a][b][c].nodepos.z;
						dz_nod2 = z3N-z1N;
						dz_nod2 *= dz_nod2;
						dis_nod2 = dx_nod2 + dy_nod2 + dz_nod2;
						if (dis_nod2 < ddmax_nod){
						dx_nod3 = x3N-x2N;
						dy_nod3 = y3N-y2N;
						dz_nod3 = z3N-z2N;
						dis_nod3 = dx_nod3*dx_nod3 + dy_nod3*dy_nod3 + dz_nod3*dz_nod3;
						if (dis_nod3 < ddmax_nod) count_3_N123(row, col, mom, u, v, w, a, b, c, SSS, nodeX);
						}
						}
					}
				}				
				}
			}	
		}
	}
	
	//  =========================  Border  =========================
	
	if (con1_x||con1_y||con1_z){ 
		
		u = row;
		v = col;
		
		fx2N = nodeX[u][0][0].nodepos.fx;
		con2_x = fx2N != 0;
		con12_x = fx1N == fx2N;
		conx_ = !con12_x && con1_x && con2_x;
		
		fy2N = nodeX[u][v][0].nodepos.fy;
		con2_y = fy2N != 0;
		con12_y = fy1N == fy2N;
		cony_ = !con12_y && con1_y && con2_y;
		
		//=======================
		// Mobile node2 in Z:
		//=======================
		for (w=mom+1;  w<partitions; ++w){	
		fz2N = nodeX[u][v][w].nodepos.fz;
		con2_z = fz2N != 0;
		con12_z = fz1N == fz2N;
		conz_ = !con12_z && con1_z && con2_z;
		
		if(con2_x || con2_y || con2_z){ // si todo es falso, entonces N2 no está en ninguna frontera
			
			// =======================================================================
			if (conx_||cony_||conz_) front_node_112(row,col,mom,u,v,w,conx_,cony_,conz_,SSS,nodeX);
			// =======================================================================
			
			a = u;
			b = v;
			
			fx3N = nodeX[a][0][0].nodepos.fx;
			con3_x = fx3N != 0;
			con13_x = fx1N == fx3N;
			conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
			
			fy3N = nodeX[a][b][0].nodepos.fy;
			con3_y = fy3N != 0;
			con13_y = fy1N == fy3N;
			cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
			
			for (c=w+1; c<partitions; ++c){
			fz3N = nodeX[a][b][c].nodepos.fz;
			con3_z = fz3N != 0;
			con13_z = fz1N == fz3N;
			conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
			if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,SSS,nodeX);
			
			}
			
			for (b=v+1; b<partitions; ++b){
			fy3N = nodeX[a][b][0].nodepos.fy;
			con3_y = fy3N != 0;
			con13_y = fy1N == fy3N;
			cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
				for (c=0; c<partitions; ++c){
				fz3N = nodeX[a][b][c].nodepos.fz;
				con3_z = fz3N != 0;
				con13_z = fz1N == fz3N;
				conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
				if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,SSS,nodeX);
				
				}
			}
			
			for (a=u+1; a<partitions; ++a){
			fx3N = nodeX[a][0][0].nodepos.fx;
			con3_x = fx3N != 0;
			con13_x = fx1N == fx3N;
			conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
				for (b=0; b<partitions; ++b){
				fy3N = nodeX[a][b][0].nodepos.fy;
				con3_y = fy3N != 0;
				con13_y = fy1N == fy3N;
				cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
					for (c=0;  c<partitions; ++c){
					fz3N = nodeX[a][b][c].nodepos.fz;
					con3_z = fz3N != 0;
					con13_z = fz1N == fz3N;
					conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
						
					if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,SSS,nodeX);
					
					}
				}
			}
		}
		}
		//=======================
		// Mobile node2 in ZY:
		//=======================
		for (v=col+1; v<partitions ; ++v){
		fy2N = nodeX[u][v][0].nodepos.fy;
		con2_y = fy2N != 0;
		con12_y = fy1N == fy2N;
		cony_ = !con12_y && con1_y && con2_y;
			for (w=0; w<partitions ; ++w){
			fz2N = nodeX[u][v][w].nodepos.fz;
			con2_z = fz2N != 0;
			con12_z = fz1N == fz2N;
			conz_ = !con12_z && con1_z && con2_z;
			
			if(con2_x || con2_y || con2_z){
			
				// =======================================================================
				if (conx_||cony_||conz_) front_node_112(row,col,mom,u,v,w,conx_,cony_,conz_,SSS,nodeX);
				// =======================================================================
				
				a = u;
				b = v;
				
				fx3N = nodeX[a][0][0].nodepos.fx;
				con3_x = fx3N != 0;
				con13_x = fx1N == fx3N;
				conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
				
				fy3N = nodeX[a][b][0].nodepos.fy;
				con3_y = fy3N != 0;
				con13_y = fy1N == fy3N;
				cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
				
				for (c=w+1; c<partitions; ++c){
				fz3N = nodeX[a][b][c].nodepos.fz;
				con3_z = fz3N != 0;
				con13_z = fz1N == fz3N;
				conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
				if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,SSS,nodeX);
				
				}
				
				for (b=v+1; b<partitions; ++b){
				fy3N = nodeX[a][b][0].nodepos.fy;
				con3_y = fy3N != 0;
				con13_y = fy1N == fy3N;
				cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
					for (c=0; c<partitions; ++c){
					fz3N = nodeX[a][b][c].nodepos.fz;
					con3_z = fz3N != 0;
					con13_z = fz1N == fz3N;
					conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
					if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,SSS,nodeX);
					
					}
				}
				
				for (a=u+1; a<partitions; ++a){
				fx3N = nodeX[a][0][0].nodepos.fx;
				con3_x = fx3N != 0;
				con13_x = fx1N == fx3N;
				conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
					for (b=0; b<partitions; ++b){
					fy3N = nodeX[a][b][0].nodepos.fy;
					con3_y = fy3N != 0;
					con13_y = fy1N == fy3N;
					cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
						for (c=0;  c<partitions; ++c){
						fz3N = nodeX[a][b][c].nodepos.fz;
						con3_z = fz3N != 0;
						con13_z = fz1N == fz3N;
						conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
						if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,SSS,nodeX);
						
						}
					}
				}
			}
			}
		}
		//=======================
		// Mobile node2 in ZYX:
		//=======================
		for (u=row+1; u<partitions ; ++u){
		fx2N = nodeX[u][0][0].nodepos.fx;
		con2_x = fx2N != 0;
		con12_x = fx1N == fx2N;
		conx_ = !con12_x && con1_x && con2_x;
			for (v=0; v<partitions ; ++v){
			fy2N = nodeX[u][v][0].nodepos.fy;
			con2_y = fy2N != 0;
			con12_y = fy1N == fy2N;
			cony_ = !con12_y && con1_y && con2_y;
				for (w=0; w<partitions ; ++w){
				fz2N = nodeX[u][v][w].nodepos.fz;
				con2_z = fz2N != 0;
				con12_z = fz1N == fz2N;
				conz_ = !con12_z && con1_z && con2_z;
				if(con2_x || con2_y || con2_z){
					
					// =======================================================================
					if (conx_||cony_||conz_) front_node_112(row,col,mom,u,v,w,conx_,cony_,conz_,SSS,nodeX);
					// =======================================================================
					
					a = u;
					b = v;
					
					fx3N = nodeX[a][0][0].nodepos.fx;
					con3_x = fx3N != 0;
					con13_x = fx1N == fx3N;
					conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
					
					fy3N = nodeX[a][b][0].nodepos.fy;
					con3_y = fy3N != 0;
					con13_y = fy1N == fy3N;
					cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
					
					for (c=w+1; c<partitions; ++c){
					fz3N = nodeX[a][b][c].nodepos.fz;
					con3_z = fz3N != 0;
					con13_z = fz1N == fz3N;
					conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
					
					if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,SSS,nodeX);
					
					}
					
					for (b=v+1; b<partitions; ++b){
					fy3N = nodeX[a][b][0].nodepos.fy;
					con3_y = fy3N != 0;
					con13_y = fy1N == fy3N;
					cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
						for (c=0; c<partitions; ++c){
						fz3N = nodeX[a][b][c].nodepos.fz;
						con3_z = fz3N != 0;
						con13_z = fz1N == fz3N;
						conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
						
						if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,SSS,nodeX);
						
						}
					}
					
					for (a=u+1; a<partitions; ++a){
					fx3N = nodeX[a][0][0].nodepos.fx;
					con3_x = fx3N != 0;
					con13_x = fx1N == fx3N;
					conx = !(con12_x && con13_x) && con1_x && con2_x && con3_x;
						for (b=0; b<partitions; ++b){
						fy3N = nodeX[a][b][0].nodepos.fy;
						con3_y = fy3N != 0;
						con13_y = fy1N == fy3N;
						cony = !(con12_y && con13_y) && con1_y && con2_y && con3_y;
							for (c=0;  c<partitions; ++c){
							fz3N = nodeX[a][b][c].nodepos.fz;
							con3_z = fz3N != 0;
							con13_z = fz1N == fz3N;
							conz = !(con12_z && con13_z) && con1_z && con2_z && con3_z;
							
							if (conx||cony||conz) front_node_123(row,col,mom,u,v,w,a,b,c,conx,cony,conz,SSS,nodeX);
							
							}
						}
					}
				}
				}
			}
		}	
	}
	}
	}
	}
	
	#pragma omp critical
	for (int i=0; i<bn; ++i){
	for (int j=0; j<bn; ++j){
	for (int k=0; k<bn; ++k) *(*(*(XXX+i)+j)+k)+=*(*(*(SSS+i)+j)+k);
	}
	}
	
	}
	//================================
	// Symmetrization:
	//================================
	symmetrize(XXX); 	
}
//=================================================================== 
void NODE3P::count_3_N111(int row, int col, int mom, double ***XXX, Node ***nodeS){
	/*
	Function to count the triangles in the same Node.

	Arg
	row, col, mom => Node1 position.
	XXX: arrangement where the DD yRR histogram will be created.
	nodeS: array of nodes.
	*/
	int i,j,k;
	float dx,dy,dz;
	double d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3;

	for (i=0; i<nodeS[row][col][mom].len-2; ++i){
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=i+1; j<nodeS[row][col][mom].len-1; ++j){
		x2 = nodeS[row][col][mom].elements[j].x;
		y2 = nodeS[row][col][mom].elements[j].y;
		z2 = nodeS[row][col][mom].elements[j].z;
		w2 = nodeS[row][col][mom].elements[j].w;
		dx = x2-x1;
		dy = y2-y1;
		dz = z2-z1;
		d12 = dx*dx+dy*dy+dz*dz;
		if (d12<=dd_max){
		d12 = sqrt(d12);
			for (k=j+1; k<nodeS[row][col][mom].len; ++k){ 
			x3 = nodeS[row][col][mom].elements[k].x;
			y3 = nodeS[row][col][mom].elements[k].y;
			z3 = nodeS[row][col][mom].elements[k].z;
			w3 = nodeS[row][col][mom].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz = z3-z1;
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<=dd_max){
			dx = x3-x2;
			dy = y3-y2;
			dz = z3-z2;
			d23 = dx*dx+dy*dy+dz*dz;
			if (d23<=dd_max){
				d13 = sqrt(d13);
				d23 = sqrt(d23);
				*(*(*(XXX+(int)(d12*ds))+(int)(d13*ds))+(int)(d23*ds))+=w1*w2*w3;
			}
			}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N112(int row, int col, int mom, int u, int v, int w, double ***XXX, Node ***nodeS){
	/*
	Function to count the triangles in two
	nodes with two points on Node1 and one point on Node2.
	
	
	Arg
	row, col, mom => Node1 position.
	u, v, w => Node2 position.
	XXX: arrangement where the DD yRR histogram will 
		be created.
	nodeS: array of nodes.
	
	*/
	int i,j,k;
	float dx,dy,dz;
	double d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3;

	for (i=0; i<nodeS[u][v][w].len; ++i){
	x1 = nodeS[u][v][w].elements[i].x;
	y1 = nodeS[u][v][w].elements[i].y;
	z1 = nodeS[u][v][w].elements[i].z;
	w1 = nodeS[u][v][w].elements[i].w;
		for (j=0; j<nodeS[row][col][mom].len; ++j){
		x2 = nodeS[row][col][mom].elements[j].x;
		y2 = nodeS[row][col][mom].elements[j].y;
		z2 = nodeS[row][col][mom].elements[j].z;
		w2 = nodeS[row][col][mom].elements[j].w;
		dx = x2-x1;
		dy = y2-y1;
		dz = z2-z1;
		d12 = dx*dx+dy*dy+dz*dz;
		if (d12<dd_max){
		d12 = sqrt(d12);
			for (k=j+1; k<nodeS[row][col][mom].len; ++k){
			x3 = nodeS[row][col][mom].elements[k].x;
			y3 = nodeS[row][col][mom].elements[k].y;
			z3 = nodeS[row][col][mom].elements[k].z;
			w3 = nodeS[row][col][mom].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz = z3-z1;
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<dd_max){
			dx = x3-x2;
			dy = y3-y2;
			dz = z3-z2;
			d23 = dx*dx+dy*dy+dz*dz;
			if (d23<dd_max){
				d13 = sqrt(d13);
				d23 = sqrt(d23);
				*(*(*(XXX+(int)(d12*ds))+(int)(d13*ds))+(int)(d23*ds))+=w1*w2*w3;
			}
			}
			}
			for (k=i+1; k<nodeS[u][v][w].len; ++k){
			x3 = nodeS[u][v][w].elements[k].x;
			y3 = nodeS[u][v][w].elements[k].y;
			z3 = nodeS[u][v][w].elements[k].z;
			w3 = nodeS[u][v][w].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz = z3-z1;
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<dd_max){
			dx = x3-x2;
			dy = y3-y2;
			dz = z3-z2;
			d23 = dx*dx+dy*dy+dz*dz;
			if (d23<dd_max){
				d13 = sqrt(d13);
				d23 = sqrt(d23);
				*(*(*(XXX+(int)(d12*ds))+(int)(d13*ds))+(int)(d23*ds))+=w1*w2*w3;
			}
			}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::count_3_N123(int row, int col, int mom, int u, int v, int w, int a, int b, int c, double ***XXX, Node ***nodeS){
	/*
	Function to count the triangles in three
	nodes with a points in Node1, a point in Node2
	and a point in Node3.
	
	Arg
	row, col, mom => Node1 position.
	u, v, w => Node2 position.
	a, b, c => Node3 position.
	XXX: arrangement where the DD yRR histogram will 
		be created.
	nodeS: array of nodes.
	
	*/
	int i,j,k;
	float dx,dy,dz;
	double d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3;
	for (i=0; i<nodeS[row][col][mom].len; ++i){
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeS[a][b][c].len; ++j){
		x3 = nodeS[a][b][c].elements[j].x;
		y3 = nodeS[a][b][c].elements[j].y;
		z3 = nodeS[a][b][c].elements[j].z;
		w3 = nodeS[a][b][c].elements[j].w;
		dx = x3-x1;
		dy = y3-y1;
		dz = z3-z1;
		d13 = dx*dx+dy*dy+dz*dz;
		if (d13<dd_max){
		d13 = sqrt(d13);
			for (k=0; k<nodeS[u][v][w].len; ++k){
			x2 = nodeS[u][v][w].elements[k].x;
			y2 = nodeS[u][v][w].elements[k].y;
			z2 = nodeS[u][v][w].elements[k].z;
			w2 = nodeS[u][v][w].elements[k].w;
			dx = x3-x2;
			dy = y3-y2;
			dz = z3-z2;
			d23 = dx*dx+dy*dy+dz*dz;
			if (d23<dd_max){
			dx = x2-x1;
			dy = y2-y1;
			dz = z2-z1;
			d12 = dx*dx+dy*dy+dz*dz;
			if (d12<dd_max){
				d12 = sqrt(d12);
				d23 = sqrt(d23);
				*(*(*(XXX+(int)(d12*ds))+(int)(d13*ds))+(int)(d23*ds))+=w1*w2*w3;
			}
			}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::symmetrize(double ***XXX){
	/*
	Function to symmetrize histogram

	Arg
	XXX: array to symmetrize
	*/ 
	int i,j,k;
	float elem;
	for (i=0; i<bn; i++){
	for (j=i; j<bn; j++){
	for (k=j; k<bn; k++){
		elem = XXX[i][j][k] + XXX[k][i][j] + XXX[j][k][i] + XXX[j][i][k] + XXX[k][j][i] + XXX[i][k][j];
		XXX[i][j][k] = elem;
		XXX[k][i][j] = elem;
		XXX[j][k][i] = elem;
		XXX[j][i][k] = elem;
		XXX[k][j][i] = elem;
		XXX[i][k][j] = elem;
	}   
	}
	}
}

//=========================================================
//================ Functions to BPC =======================
//=========================================================
void NODE3P::front_node_112(int row, int col, int mom, int u, int v, int w, bool conx, bool cony, bool conz, double ***XXX, Node ***nodeS){
	
	float x1N = nodeS[row][0][0].nodepos.x;
	float y1N = nodeS[row][col][0].nodepos.y;
	float z1N = nodeS[row][col][mom].nodepos.z;
	float fx2N = nodeS[u][0][0].nodepos.fx, x2N = nodeS[u][0][0].nodepos.x;
	float fy2N = nodeS[u][v][0].nodepos.fy, y2N = nodeS[u][v][0].nodepos.y;
	float fz2N = nodeS[u][v][w].nodepos.fz, z2N = nodeS[u][v][w].nodepos.z;
	
	float dx, dy, dz;
	
	if (conx){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,fx2N,0,0,XXX,nodeS);
	}
	if (cony){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,0,fy2N,0,XXX,nodeS);
	}
	if (conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz =(z2N-(fz2N*size_box)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,0,0,fz2N,XXX,nodeS);
	}
	if (conx && cony){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,fx2N,fy2N,0,XXX,nodeS);
	}
	if (conx && conz){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = (z2N-fz2N*size_box) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,fx2N,0,fz2N,XXX,nodeS);
	}
	if (cony && conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = (z2N-fz2N*size_box) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,0,fy2N,fz2N,XXX,nodeS);
	}
	if (conx && cony && conz){
		dx = (x2N-fx2N*size_box) - x1N;
		dx *= dx;
		dy = (y2N-fy2N*size_box) - y1N;
		dy *= dy; 
		dz = (z2N-fz2N*size_box) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod) front_112(row,col,mom,u,v,w,fx2N,fy2N,fz2N,XXX,nodeS);
	}
}
//=================================================================== 
void NODE3P::front_112(int row,int col,int mom,int u,int v,int w,short int fx2,short int fy2,short int fz2,double ***XXX,Node ***nodeS){

	int i,j,k;
	int n,m,l;
	float dx,dy,dz;
	double d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3,w1,w2,w3;

	for (i=0; i<nodeS[row][col][mom].len; ++i){
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeS[u][v][w].len; ++j){
		x2 = nodeS[u][v][w].elements[j].x - fx2*size_box;
		y2 = nodeS[u][v][w].elements[j].y - fy2*size_box;
		z2 = nodeS[u][v][w].elements[j].z - fz2*size_box;
		w2 = nodeS[u][v][w].elements[j].w;
		dx = x2-x1;
		dy = y2-y1;
		dz = z2-z1;
		d12 = dx*dx+dy*dy+dz*dz;
		if (d12<dd_max){
		d12 = sqrt(d12);
			for (k=i+1; k<nodeS[row][col][mom].len; ++k){
			x3 = nodeS[row][col][mom].elements[k].x;
			y3 = nodeS[row][col][mom].elements[k].y;
			z3 = nodeS[row][col][mom].elements[k].z;
			w3 = nodeS[row][col][mom].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz = z3-z1;
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<dd_max){
			dx = x3-x2;
			dy = y3-y2;
			dz = z3-z2;
			d23 = dx*dx+dy*dy+dz*dz;
			if (d23<dd_max){
				d13 = sqrt(d13);
				d23 = sqrt(d23);
				n = (int)(d12*ds);
				m = (int)(d13*ds);
				l = (int)(d23*ds);
				*(*(*(XXX+n)+m)+l)+=w1*w2*w3;
			}
			}
			}
			for (k=j+1; k<nodeS[u][v][w].len; ++k){
			x3 = nodeS[u][v][w].elements[k].x - fx2*size_box;
			y3 = nodeS[u][v][w].elements[k].y - fy2*size_box;
			z3 = nodeS[u][v][w].elements[k].z - fz2*size_box;
			w3 = nodeS[u][v][w].elements[k].w;
			dx = x3-x1;
			dy = y3-y1;
			dz = z3-z1;
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<dd_max){
			dx = x3-x2;
			dy = y3-y2;
			dz = z3-z2;
			d23 = dx*dx+dy*dy+dz*dz;
			if (d23<dd_max){
				d13 = sqrt(d13);
				d23 = sqrt(d23);
				n = (int)(d12*ds);
				m = (int)(d13*ds);
				l = (int)(d23*ds);
				*(*(*(XXX+n)+m)+l)+=w1*w2*w3;
			}
			}
			}
			}
		}
	}
}
//=================================================================== 
void NODE3P::front_node_123(int row, int col, int mom, int u, int v, int w, int a, int b, int c, bool conx, bool cony, bool conz, double ***XXX, Node ***nodeS){
	
	float fx1N = nodeS[row][0][0].nodepos.fx, x1N = nodeS[row][0][0].nodepos.x;
	float fy1N = nodeS[row][col][0].nodepos.fy, y1N = nodeS[row][col][0].nodepos.y;
	float fz1N = nodeS[row][col][mom].nodepos.fz, z1N = nodeS[row][col][mom].nodepos.z;
	float fx2N = nodeS[u][0][0].nodepos.fx, x2N = nodeS[u][0][0].nodepos.x;
	float fy2N = nodeS[u][v][0].nodepos.fy, y2N = nodeS[u][v][0].nodepos.y;
	float fz2N = nodeS[u][v][w].nodepos.fz, z2N = nodeS[u][v][w].nodepos.z;
	float fx3N = nodeS[a][0][0].nodepos.fx, x3N = nodeS[a][0][0].nodepos.x;
	float fy3N = nodeS[a][b][0].nodepos.fy, y3N = nodeS[a][b][0].nodepos.y;
	float fz3N = nodeS[a][b][c].nodepos.fz, z3N = nodeS[a][b][c].nodepos.z;
	
	float dx,dy,dz;
	
	if (conx){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = y3N-y1N;
			dy *= dy; 
			dz = z3N-z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = y3N-y2N;
				dy *= dy; 
				dz = z3N-z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,fx1N,0,0,fx2N,0,0,fx3N,0,0,XXX,nodeS);
				}
			}
		}
	
	}
	if (cony){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = x3N-x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = z3N-z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = x3N-x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = z3N-z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,0,fy1N,0,0,fy2N,0,0,fy3N,0,XXX,nodeS);
				}
			}
		}
	}
	if (conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz =(z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = x3N-x1N;
			dx *= dx;
			dy = y3N-y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = x3N-x2N;
				dx *= dx;
				dy = y3N-y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,0,0,fz1N,0,0,fz2N,0,0,fz3N,XXX,nodeS);
				}
			}
		}
	}
	if (conx && cony){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = z2N-z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = z3N-z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = z3N-z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,fx1N,fy1N,0,fx2N,fy2N,0,fx3N,fy3N,0,XXX,nodeS);
				}
			}
		}
	}
	if (conx && conz){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = y2N-y1N;
		dy *= dy; 
		dz = (z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = y3N-y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = y3N-y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,fx1N,0,fz1N,fx2N,0,fz2N,fx3N,0,fz3N,XXX,nodeS);
				}
			}
		}
	}
	if (cony && conz){
		dx = x2N-x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = (z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = x3N-x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = x3N-x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,0,fy1N,fz1N,0,fy2N,fz2N,0,fy3N,fz3N,XXX,nodeS);
				}
			}
		}
	}
	if (conx && cony && conz){
		dx = (x2N-((fx2N-fx1N)*size_box_2)) - x1N;
		dx *= dx;
		dy = (y2N-((fy2N-fy1N)*size_box_2)) - y1N;
		dy *= dy; 
		dz = (z2N-((fz2N-fz1N)*size_box_2)) - z1N;
		dz *= dz;
		if (dx+dy+dz < ddmax_nod){
			dx = (x3N-((fx3N-fx1N)*size_box_2)) - x1N;
			dx *= dx;
			dy = (y3N-((fy3N-fy1N)*size_box_2)) - y1N;
			dy *= dy; 
			dz = (z3N-((fz3N-fz1N)*size_box_2)) - z1N;
			dz *= dz;
			if (dx+dy+dz < ddmax_nod){
				dx = (x3N-((fx3N-fx2N)*size_box_2)) - x2N;
				dx *= dx;
				dy = (y3N-((fy3N-fy2N)*size_box_2)) - y2N;
				dy *= dy; 
				dz = (z3N-((fz3N-fz2N)*size_box_2)) - z2N;
				dz *= dz;
				if (dx+dy+dz < ddmax_nod){
					front_123(row,col,mom,u,v,w,a,b,c,fx1N,fy1N,fz1N,fx2N,fy2N,fz2N,fx3N,fy3N,fz3N,XXX,nodeS);
				}
			}
		}
	}	
}
//=================================================================== 
void NODE3P::front_123(int row,int col,int mom,int u,int v,int w,int a,int b,int c,short int fx1N,short int fy1N,short int fz1N,short int fx2N,short int fy2N,short int fz2N,short int fx3N,short int fy3N,short int fz3N,double ***XXX, Node ***nodeS){
	
	int i,j,k;
	int n,m,l;
	float dx,dy,dz;
	double d12,d13,d23;
	float x1,y1,z1,x2,y2,z2,x3,y3,z3;
	float w1,w2,w3;
	for (i=0; i<nodeS[row][col][mom].len; ++i){
	x1 = nodeS[row][col][mom].elements[i].x;
	y1 = nodeS[row][col][mom].elements[i].y;
	z1 = nodeS[row][col][mom].elements[i].z;
	w1 = nodeS[row][col][mom].elements[i].w;
		for (j=0; j<nodeS[u][v][w].len; ++j){
		x2 = nodeS[u][v][w].elements[j].x;
		y2 = nodeS[u][v][w].elements[j].y;
		z2 = nodeS[u][v][w].elements[j].z;
		w2 = nodeS[u][v][w].elements[j].w;
		dx = (x2-((fx2N-fx1N)*size_box_2))-x1;
		dy = (y2-((fy2N-fy1N)*size_box_2))-y1;
		dz = (z2-((fz2N-fz1N)*size_box_2))-z1;
		d12 = dx*dx+dy*dy+dz*dz;
		if (d12<dd_max){
		d12 = sqrt(d12);
			for (k=0; k<nodeS[a][b][c].len; ++k){
			x3 = nodeS[a][b][c].elements[k].x;
			y3 = nodeS[a][b][c].elements[k].y;
			z3 = nodeS[a][b][c].elements[k].z;
			w3 = nodeS[a][b][c].elements[k].w;
			dx = (x3-((fx3N-fx1N)*size_box_2))-x1;
			dy = (y3-((fy3N-fy1N)*size_box_2))-y1;
			dz = (z3-((fz3N-fz1N)*size_box_2))-z1;
			d13 = dx*dx+dy*dy+dz*dz;
			if (d13<dd_max){
			dx = (x3-((fx3N-fx2N)*size_box_2))-x2;
			dy = (y3-((fy3N-fy2N)*size_box_2))-y2;
			dz = (z3-((fz3N-fz2N)*size_box_2))-z2;
			d23 = dx*dx+dy*dy+dz*dz;
			if (d23<dd_max){
				d13 = sqrt(d13);
				d23 = sqrt(d23);
				n = (int)(d12*ds);
				m = (int)(d13*ds);
				l = (int)(d23*ds);
				*(*(*(XXX+n)+m)+l)+=w1*w2*w3;
			}
			}
			}
			}
		}
	}
}

//========================================================================== 
//================ Functions for analytical formulas =======================
//==========================================================================

void NODE3P::make_histo_analitic(double ***XXY, double ***XXX, Node ***nodeX){
	/*
	Function to construct histograms with analytical functions

	Arguments
	XXX_: random array
	XXY: mixed arrangement
	nodeX: data mesh
	
	*/ 
	//=======================================
	// RRR, DDR and DRR Histogram (ANALYTICS)
	//=======================================
	
	std::cout << "Constructing RRR and DDR histograms..." << std::endl;
	
	int i, j, k, u, v, w, a, b, c;
	
	//====================================================
	// Constants for RRR
	double dr = d_max/(double)(bn);
	double ri, rj, rk;
	double V = size_box*size_box;
	double beta = n_pts/V;
	double gama  = 8*(3.14159265359*3.14159265359)*(n_pts*beta*beta);
	gama /= V;
	double alph = gama*dr*dr*dr;
	
	// Refinement for RRR
	int ptt_r = 200;
	int bn_ref_r = ptt_r*bn; 
	double dr_ref_r = dr/(double)(ptt_r);
	
	//====================================================
	// Constants for DDR
	int ptt = 1000;
	int bins = ptt*bn;
	double dr_ptt = d_max/(double)(bins);
	double *DD;
	double *RR;
	DD = new double[bins];
	RR = new double[bins];
	for (i=0; i<bins; ++i){
		*(DD+i) = 0.0; 
		*(RR+i) = 0.0;
	}
	
	// Make the DD and RR histograms
	//====================================================
	make_histoXX(DD, RR, nodeX, bins);
	
	//Imprimir DD y RR
	std::cout << " "  << std::endl;
	std::cout << "   DD   " << std::endl;
	int stp = bins/15;
	
	for (i=0; i<stp; i++){
		for (j=0; j<15; j++) std::cout << *(DD+(i*15)+j) << " ";
		std::cout << " "  << std::endl;
	}
	std::cout << " "  << std::endl;
	std::cout << "   RR   "  << std::endl;
	for (i=0; i<stp; i++){
		for (j=0; j<15; j++) std::cout << *(RR+(i*15)+j) << " ";
		std::cout << " "  << std::endl;
	}
	
	//====================================================
	// Initiate an arrangement for the function f
	double *f;
	f = new double[bins];
	for (i=0; i<bins; ++i) *(f+i) = 0.0;
	
	// Calculate the f
	#pragma omp parallel num_threads(50) shared(DD,RR)
	{
	double *SS;
    	SS = new double[bins];
    	for (int k=0; k<bins; k++) *(SS+k) = 0.0;
    	
    	double ri;
    	
	#pragma omp for schedule(dynamic)
	for(int i=0; i<bins; ++i){
		ri = (i+0.5)*dr_ptt;
		*(SS+i) += ri*((*(DD+i)/(*(RR+i))) - 1);
	}
	#pragma omp critical
	for(int a=0; a<bins; ++a) *(f+a) += *(SS+a);
	delete[] SS;
	}
	
	delete[] DD;
	delete[] RR;
	
	//====================================================
	// Calculate the f_averrage
	double *f_v;
	f_v = new double[bn];
	for (i=0; i<bn; ++i) *(f_v+i) = 0.0;
	
	int i_;
	double ff;
	
	for(i=0; i<bn; ++i){
	i_ = i*ptt;
	ff = 0.0;
		#pragma omp parallel num_threads(50) shared(ff)
		#pragma omp for reduction(+:ff)
		for(int j=0; j<ptt; ++j)  ff += *(f+i_+j);
	*(f_v+i) += ff/(double)(ptt);
	}
	
	//====================================================
	int j_;
	int n;
	double alph_ref = gama*dr_ref_r*dr_ref_r*dr_ref_r;
	double dr2 = dr/2;
	double dr_ptt2 = dr_ptt/2;
	double dr_ref2 = dr_ref_r/2;
	
	std::cout << "Constructing RRR and DDR histograms..." << std::endl;
	
    	int short v_in;
	double r1, r2, r3;
	double ru, rv, rw;
	double S;
	double ff_ave;
	bool con;
	double c_RRR;
	
	for(i=0; i<bn; ++i) {
	ri = i*dr;
	i_ = i*ptt;
	for(j=i; j<bn; ++j) {
	rj = j*dr;
	j_ = j*ptt;
	for(k=j; k<bn; ++k) {
	rk = k*dr;
		
		v_in = 0;
		
		// Check vertices of the 
		// cube to make refinement
		for (a = 0; a < 2; ++a){
		r1 = ri + (a*dr);
		for (b = 0; b < 2; ++b){
		r2 = rj + (b*dr);
		for (c = 0; c < 2; ++c){
		r3 = rk + (c*dr);	
			if (r1 + r2 > r3 && r1 + r3 > r2 && r2 + r3 > r1) ++v_in; 
		}
		}
		}
		
		// bins inequality of the triangle true
		if (v_in==8){
			S = alph*(ri+dr2)*(rj+dr2)*(rk+dr2);
			ff = 1;
			ff += *(f_v+i)/(3*(ri+dr2));
			ff += *(f_v+j)/(3*(rj+dr2));
			ff += *(f_v+k)/(3*(rk+dr2));
			ff *= S;
			*(*(*(XXY+i)+j)+k) += ff;
			*(*(*(XXX+i)+j)+k) += S;
		}
		
		// bins inequality of the triangle false
		else if (v_in < 8 && v_in > 0){
			
			//====================================================================
			
			con = false;
			c_RRR = 0.0;
			
			for(int u=0; u<ptt_r; ++u) {
			ru = ri + (u*dr_ref_r);
			for(int v=0; v<ptt_r; ++v) {
			rv = rj + (v*dr_ref_r);
			for(int w=0; w<ptt_r; ++w) {
			rw = rk + (w*dr_ref_r);
			
				v_in = 0;
				
				// Check vertices of the 
				// cube to make refinement
				for (int a=0; a<2; ++a){
				r1 = ru + (a*dr_ref_r);
				for (int b=0; b<2; ++b){
				r2 = rv + (b*dr_ref_r);
				for (int c=0; c<2; ++c){
				r3 = rw + (c*dr_ref_r);	
					if (r1 + r2 > r3 && r1 + r3 > r2 && r2 + r3 > r1) ++v_in;
				}
				}
				}
				
				// bins inequality of the triangle true
				if (v_in==8){
					c_RRR += (ru+dr_ref2)*(rv+dr_ref2)*(rw+dr_ref2);
					con = true;
				}
			}
			}
			}
			if (con) *(*(*(XXX+i)+j)+k) += alph_ref*c_RRR;
			
			//====================================================================
			
			con = false;
			ff_ave = 0.0;
			n = 0;
			
			for(int u=0; u<ptt; ++u) {
			ru = ri + (u*dr_ptt);
			for(int v=0; v<ptt; ++v) {
			rv = rj + (v*dr_ptt);
			for(int w=0; w<ptt; ++w) {
			rw = rk + (w*dr_ptt);
			
				v_in = 0;
				
				// Check vertices of the 
				// cube to make refinement
				for (int a=0; a<2; ++a){
				r1 = ru + (a*dr_ptt);
				for (int b=0; b<2; ++b){
				r2 = rv + (b*dr_ptt);
				for (int c=0; c<2; ++c){
				r3 = rw + (c*dr_ptt);	
					if (r1 + r2 > r3 && r1 + r3 > r2 && r2 + r3 > r1) ++v_in;
				}
				}
				}
				
				// bins inequality of the triangle true
				if (v_in==8){
					ff = 1;
					ff += *(f+i_+u)/(3*(ru+dr_ptt2));
					ff += *(f+j_+v)/(3*(rv+dr_ptt2));
					ff += *(f+(k*ptt)+w)/(3*(rw+dr_ptt2));
					ff_ave += ff;
					con = true;
					n++; 
				}
			}
			}
			}
			if (con) *(*(*(XXY+i)+j)+k) += *(*(*(XXX+i)+j)+k)*(ff_ave/n);
		}
	}
	}
	}
	
	symmetrize_analitic(XXX);
	symmetrize_analitic(XXY); 
	
}
//=================================================================== 
void NODE3P::make_histoXX(double *XX, double *YY, Node ***nodeX, int bins){
	/*
	Function to create the DD and RR histograms.

	Arg:
	DD: array where the DD histogram will be created.
	RR: array where the RR histogram will be created.
	
	*/
	
	std::cout << "Constructing DD..." << std::endl;
	
	double ds_new = ((float)(bins))/d_max;
	
	#pragma omp parallel num_threads(50)
	{
	
	int i, j, k, row, col, mom, u, v, w;
	
	double *SS;
    	SS = new double[bins];
    	for (k = 0; k < bins; ++k) *(SS+k) = 0.0;
    	
	// Private variables in threads:
	double dis;
	float dis_nod;
	float x1D, y1D, z1D, x2D, y2D, z2D;
	float x, y, z, a;
	float dx, dy, dz, dx_nod, dy_nod, dz_nod;
	bool con_x, con_y, con_z;
	
	#pragma omp for collapse(3) schedule(dynamic)
	for (row = 0; row < partitions; ++row){
	for (col = 0; col < partitions; ++col){
	for (mom = 0; mom < partitions; ++mom){
	x1D = nodeX[row][0][0].nodepos.x;
	y1D = nodeX[row][col][0].nodepos.y;
	z1D = nodeX[row][col][mom].nodepos.z;			
		//==================================================
		// Pairs of points in the same node:
		//==================================================
		for (i= 0; i<nodeX[row][col][mom].len-1; ++i){
		x = nodeX[row][col][mom].elements[i].x;
		y = nodeX[row][col][mom].elements[i].y;
		z = nodeX[row][col][mom].elements[i].z;
		a = nodeX[row][col][mom].elements[i].w;
			for (j=i+1; j<nodeX[row][col][mom].len; ++j){
			dx = x-nodeX[row][col][mom].elements[j].x;
			dy = y-nodeX[row][col][mom].elements[j].y;
			dz = z-nodeX[row][col][mom].elements[j].z;
			dis = dx*dx+dy*dy+dz*dz;
			if (dis < dd_max) *(SS + (int)(sqrt(dis)*ds_new)) += 2*a*nodeX[row][col][mom].elements[j].w;
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
		dz_nod *= dz_nod;
		if (dz_nod <= ddmax_nod){
			for (i=0; i<nodeX[row][col][mom].len; ++i){
			x = nodeX[row][col][mom].elements[i].x;
			y = nodeX[row][col][mom].elements[i].y;
			z = nodeX[row][col][mom].elements[i].z;
			a = nodeX[row][col][mom].elements[i].w;
				for (j=0; j<nodeX[u][v][w].len; ++j){
				dx = x-nodeX[u][v][w].elements[j].x;
				dy = y-nodeX[u][v][w].elements[j].y;
				dz = z-nodeX[u][v][w].elements[j].z;
				dis = dx*dx+dy*dy+dz*dz;
				if (dis < dd_max) *(SS + (int)(sqrt(dis)*ds_new)) += 2*a*nodeX[u][v][w].elements[j].w;
				}
			}
			}
			// ======================================= 
			// Distance of border points XX 
			// ======================================= 
			// Boundary node conditions:
			con_z = ((z1D<=d_max_pm)&&(z2D>=front_pm))||((z2D<=d_max_pm)&&(z1D>=front_pm));
			if(con_z){
			histo_front_XX(SS,nodeX,dis_nod,0.0,0.0,fabs(dz_nod),false,false,con_z,row,col,mom,u,v,w,ds_new);
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
				a = nodeX[row][col][mom].elements[i].w;
					for ( j = 0; j < nodeX[u][v][w].len; ++j){	
					dx =  x-nodeX[u][v][w].elements[j].x;
					dy =  y-nodeX[u][v][w].elements[j].y;
					dz =  z-nodeX[u][v][w].elements[j].z;
					dis = dx*dx+dy*dy+dz*dz;
					if (dis < dd_max) *(SS + (int)(sqrt(dis)*ds_new)) += 2*a*nodeX[u][v][w].elements[j].w;
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
			histo_front_XX(SS,nodeX,dis_nod,0.0,fabs(dy_nod),fabs(dz_nod),false,con_y,con_z,row,col,mom,u,v,w,ds_new);
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
					a = nodeX[row][col][mom].elements[i].w;
						for ( j = 0; j < nodeX[u][v][w].len; ++j){	
						dx = x-nodeX[u][v][w].elements[j].x;
						dy = y-nodeX[u][v][w].elements[j].y;
						dz = z-nodeX[u][v][w].elements[j].z;
						dis = dx*dx + dy*dy + dz*dz;
						if (dis < dd_max) *(SS + (int)(sqrt(dis)*ds_new)) += 2*a*nodeX[u][v][w].elements[j].w;
						}
					}
				}
				//======================================= 
				// Distance of border points XX 
				//======================================= 
				// Boundary node conditions:
				con_x = ((x1D<=d_max_pm)&&(x2D>=front_pm))||((x2D<=d_max_pm)&&(x1D>=front_pm));
				con_y = ((y1D<=d_max_pm)&&(y2D>=front_pm))||((y2D<=d_max_pm)&&(y1D>=front_pm));
				con_z = ((z1D<=d_max_pm)&&(z2D>=front_pm))||((z2D<=d_max_pm)&&(z1D>=front_pm));
			if(con_x || con_y || con_z){
			histo_front_XX(SS,nodeX,dis_nod,fabs(dx_nod),fabs(dy_nod),fabs(dz_nod),con_x,con_y,con_z,row,col,mom,u,v,w,ds_new);
			}	
				}	
			}
		}
	}
	}
	}
	#pragma omp critical
	for(int a=0; a<bins; ++a) *(XX+a)+=*(SS+a);
	}
	
	std::cout << "Constructing RR..." << std::endl;
	
	//======================================
	// Histogram RR (ANALYTICAL)
	//======================================
	
	double dr = d_max/(double)(bins);
	double V = size_box;
	double beta1 = n_pts/V;
	beta1 *= n_pts/size_box;
	beta1 /= size_box;
	double alph = 4*(2*acos(0.0))*(beta1)/3;
	alph *= dr*dr*dr;
	
	#pragma omp parallel num_threads(50)
	{
	
	double *SS;
    	SS = new double[bins];
    	for (int k = 0; k < bins; k++) *(SS+k) = 0.0;
	
	double r1, r2;
	#pragma omp for schedule(dynamic)	
	for(int a=0; a<bins; ++a) {
		r2 = (double)(a);
		r1 = r2+1;
        	*(SS+a) += alph*((r1*r1*r1)-(r2*r2*r2));
	}
	#pragma omp critical
	for(int a=0; a<bins; ++a) *(YY+a)+=*(SS+a);
	
	delete[] SS;
	}
}

//=================================================================== 
void NODE3P::histo_front_XX(double *PP, Node ***dat, float disn, float dn_x, float dn_y, float dn_z, bool con_in_x, bool con_in_y, bool con_in_z, int row, int col, int mom, int u, int v, int w, double ds_new){
	int i, j;
	double dis_f, dis;
	float d_x, d_y, d_z;
	float x, y, z, a;
	//======================================================================
	if( con_in_x ){
	dis_f = disn + ll - 2*dn_x*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = y-dat[u][v][w].elements[j].y;
			d_z = z-dat[u][v][w].elements[j].z;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis < dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += 2*a*dat[u][v][w].elements[j].w;
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
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = x-dat[u][v][w].elements[j].x;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = z-dat[u][v][w].elements[j].z;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis < dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += 2*a*dat[u][v][w].elements[j].w;
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
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = x-dat[u][v][w].elements[j].x;
			d_y = y-dat[u][v][w].elements[j].y;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis < dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += 2*a*dat[u][v][w].elements[j].w;
			}
		}
	}
	}
	//======================================================================			
	if( con_in_x && con_in_y ){
	dis_f = disn + 2*ll - 2*(dn_x+dn_y)*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = z-dat[u][v][w].elements[j].z;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis < dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += 2*a*dat[u][v][w].elements[j].w;
			}
		}
	}
	}
	//======================================================================			
	if( con_in_x && con_in_z ){
	dis_f = disn + 2*ll - 2*(dn_x+dn_z)*size_box;
	if (dis_f <= ddmax_nod){
		for (i=0; i<dat[row][col][mom].len; ++i){
		x = dat[row][col][mom].elements[i].x;
		y = dat[row][col][mom].elements[i].y;
		z = dat[row][col][mom].elements[i].z;
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = y-dat[u][v][w].elements[j].y;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis < dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += 2*a*dat[u][v][w].elements[j].w;
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
		a = dat[row][col][mom].elements[i].w;
			for ( j = 0; j < dat[u][v][w].len; ++j){
			d_x = x-dat[u][v][w].elements[j].x;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = (d_x*d_x) + (d_y*d_y) + (d_z*d_z); 
			if (dis < dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += 2*a*dat[u][v][w].elements[j].w;
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
		a = dat[row][col][mom].elements[i].w;
			for (j=0; j<dat[u][v][w].len; ++j){
			d_x = fabs(x-dat[u][v][w].elements[j].x)-size_box;
			d_y = fabs(y-dat[u][v][w].elements[j].y)-size_box;
			d_z = fabs(z-dat[u][v][w].elements[j].z)-size_box;
			dis = d_x*d_x + d_y*d_y + d_z*d_z;
			if (dis < dd_max) *(PP + (int)(sqrt(dis)*ds_new)) += 2*a*dat[u][v][w].elements[j].w;
			}
		}
	}
	}
}
//=================================================================== 
void NODE3P::symmetrize_analitic(double ***XXX){
	/*
	Function to symmetrize histogram

	Arg
	XXX: array to symmetrize
	*/ 

	int i,j,k;
	float elem;
	for (i=0; i<bn; i++){
	for (j=0; j<bn; j++){
	for (k=0; k<bn; k++){
		elem = XXX[i][j][k];
		XXX[k][i][j] = elem;
		XXX[j][k][i] = elem;
		XXX[j][i][k] = elem;
		XXX[k][j][i] = elem;
		XXX[i][k][j] = elem;
	}   
	}
	}
}
//=================================================================== 
NODE3P::~NODE3P(){	
}


