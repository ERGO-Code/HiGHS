//#include <cmath>
//#include <cstring>
//#include <cstdio>
//#include <fstream>
//#include <iostream>
//using namespace std;

int readToy_LP_cpp(const char *filename, int* m_p, int* n_p, int* maxmin, double* offset,
		 double ** A,
		 double ** b, double ** c, double ** lb, double ** ub);

int readToy_MIP_cpp(const char *filename, int* numRow_p, int* numCol_p, int* objSense_p, double* objOffset_p,
		  double ** A_cw_p,
		  double ** rhs_p, double ** cost_p, double ** lb_p, double ** ub_p,
		  int** integerColumn);
