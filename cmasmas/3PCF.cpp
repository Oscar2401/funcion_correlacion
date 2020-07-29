// jul 27 Oscar De la Cruz Echeveste
// Acomodar nuestros en un arreglo 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;


int main()
{
	//Variables y vectores
	float inputData;
	vector<float> data_v;
	vector<float> rand_v;
	vector<vector<float>> data_v2;
	vector<vector<float>> rand_v2;

	//Declaramos una variable ifstream
	ifstream data, rand;
	
	//Leemos los datos
	data.open("/home/echeveste/mis_trabajos/correlation_f/data/data_500.dat");
	rand.open("/home/echeveste/mis_trabajos/correlation_f/data/ran0_500.dat");

	if (!data && !rand) {
		cout << "No fue posible abrir el archivo \n";
		exit(1); // Terminar con error
	}

	while(data>>inputData){ //Lee un dato a la vez
		data_v.push_back(inputData); // Adiere cada elemento al vector
	}

	while(rand>>inputData){ //Lee un dato a la vez
                rand_v.push_back(inputData); // Adiere cada elemento al vector   
	}

	data.close();
	rand.close();
	
	//reacomodamos los datos en un arreglo de 4x500
	data_v2.resize(4,data_v);
	rand_v2.resize(4,rand_v);
	cout << data_v2[3][3] << endl;
	
	vector<vector<vector<float>>> DDD;
	DDD.assign(5,50);
	
	
	
	return 0;
}
