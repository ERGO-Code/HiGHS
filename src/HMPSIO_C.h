
int readMPS_dense_fc(const char *filename, int* numRow_p, int* numCol_p, int* objSense_p, double* objOffset_p,
		    double ** A_cw_p,
		    double ** rhs_p, double ** cost_p, double ** lb_p, double ** ub_p,
		    int** integerColumn_p);

