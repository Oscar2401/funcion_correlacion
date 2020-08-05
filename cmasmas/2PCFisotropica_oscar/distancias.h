#include <cmath>


// Distancia euclidea
template <typename TDG>
TDG dist(TDG x, TDG y, TDG z){
	return sqrt(x*x + y*y);
}

// Máxima distancia en el x
template <typename TDG1>
TDG1 max_x(TDG1 data, TDG1 num_puntos){
	float max;
	for (int i = 1; i < num_puntos; i++){
		if (data[i].x > data[i-1].x){
			max = data[i].x;
		}
		else{
			max = data[i-1].x;
		}
	}
	return max;
}

// Máxima distancia en el y
template <typename TDG2>
TDG2 max_y(TDG2 data, TDG2 num_puntos){
	float max;
	for (int i = 1; i < num_puntos; i++){
		if (data[i].y > data[i-1].y){
			max = data[i].y;
		}
		else{
			max = data[i-1].y;
		}
	}
	return max;
}

// Mínima distancia en el x
template <typename TDG3>
TDG3 min_x(TDG3 data, TDG3 num_puntos){
	float min;
	for (int i = 1; i < num_puntos; i++){
		if (data[i].x < data[i-1].x){
			min = data[i].x;
		}
		else{
			min = data[i-1].x;
		}
	}
	return min;
}

// Mínima distancia en el y
template <typename TDG4>
TDG4 min_y(TDG4 data, TDG4 num_puntos){
	float min;
	for (int i = 1; i < num_puntos; i++){
		if (data[i].y < data[i-1].y){
			min = data[i].y;
		}
		else{
			min = data[i-1].y;
		}
	}
	return min;
}

template <typename TDG5>
TDG5 make_mesh(TDG5 data, TDG5 mesh, TDG5 PM, TDG5 num_puntos, TDG5 tam){
	
	float d_max_x = ceil(max_x( data, num_puntos));
	float d_max_y = ceil(max_y( data, num_puntos));
	float d_min_x = floor(min_x( data, num_puntos));
	float d_min_y = floor(min_y( data, num_puntos));
	
	std::cout << d_max_x << std::endl;
	
	int n_x = ceil((d_max_x-d_min_x)/tam);
	int n_y = ceil((d_max_y-d_min_y)/tam);
	
	for (int i = 0; i < num_puntos-1; i++){
		for (int u = 0;  u < n_x-1; u++){
			for (int v = 0;  v < n_y-1; v++){
				if (data[i] > u*d_min_x && data[i] < (u+1)*d_max_x && data[i] > v*d_min_y && data[i] < (v+1)*d_max_y){
					mesh[u+(v*n_x)] = data[i];
				}
			} 
		}
	}
}

