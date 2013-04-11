#include <iostream>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "FastGraphCluster.h"
using namespace std;

int main(int argc,char *argv[])
{
	string file, clusterFile;
	size_t min_cluster_size = 2;
	double min_density = 0.5;
	double min_increase = 0.5;
	char mode = '0';	// mode of graph representation
	string type,value;

	if ( argc < 5 )
	{
		if(argc == 2 && (value = argv[1]) == "-h")
		{
			cout << "\nExtremely fast graph clustering algorithm. (Peng Jiang 2009)\n" << endl;
			cout << "Usage: spici [OPTIONS]... [FILES]...\n"<<endl;
			cout << "\t-h Print help and exit\n"<<endl;
			cout << "Main:\n"<<endl;
			cout << "\t-i\tInput graph file"<<endl;
			cout << "\t-o\tOutput cluster results\n"<<endl;
			cout << "Parameters:\n"<<endl;
			cout << "\t-d\tminimum density threshold. Default: 0.5"<<endl;
			cout << "\t-s\tminimum cluster size. Default: 2"<<endl;
			cout << "\t-g\tminimum increment ratio. Default: 0.5"<<endl;
			cout << "\t-m\tGraph mode. Default:0.\n"<<endl;
			cout << "\t\t0: sparse graph"<<endl;
			cout << "\t\t1: dense graph"<<endl;
			cout << "\t\t2: large sparse graph"<<endl;
			cout << "\nReport bugs to peng.jiang.software@gmail.com\n"<<endl;
			exit(0);
		}else{
			cerr << "Wrong arguments input. \n-h for help"<<endl;
			exit(1);
		}
	}

	// read in all parameters
	size_t i,parseCnt = (argc - 1)/2;

	for(i=0;i<parseCnt;i++){
		type = argv[2*i+1];
		value = argv[2*i+2];

		if(type == "-i"){
			file = value;
		}else if (type == "-o"){
			clusterFile = value;
		}else if (type == "-d"){
			min_density = atof(value.c_str());
		}else if (type == "-s"){
			min_cluster_size = atoi(value.c_str());
		}else if (type == "-g"){
			min_increase = atof(value.c_str());
		}else if (type == "-m"){
			mode = value[0];
		}else{
			cerr << "Cannot recognize:\t"<< type <<endl;
			cerr << "Wrong arguments input. \n-h for help"<<endl;
			exit(1);
		}
	}

	if ( mode != '0' && mode != '1' && mode != '2' ){
		cerr << "Invalid graph mode! can only be 0:sparse\t1:nearly complete graph\t2:somewhere between 0 and 1"<<endl;
		exit(1);
	}

	FastGraphCluster cluster(file,min_density,min_cluster_size,min_increase,mode);
	int ncluster = cluster.fastCluster(clusterFile);

	cout <<"There are "<<ncluster<<" clusters."<<endl;
	return 0;
}

