// Header file for input output functions
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define kB 0.0000015
#define mag(x) sqrt(x[0]*x[0] + x[1]*x[1])


typedef struct _Neutron Neutron;
struct _Neutron {
	float X[2];
	float V[2];
};

void printN(Neutron n){
	printf("Neutron:\n");
	printf("  X: %e, %e\n", n.X[0], n.X[1]);
	printf("  V: %e, %e\n", n.V[0], n.V[1]);
}

void printNE(Neutron n[], int LEN_n){
	printf("Ensamble:\n");
	for (int i = 0; i < LEN_n; i++){
		printf("  ");
		printN(n[i]);
	}
}

void motion(Neutron* n, float dt){
	n->X[0] += n->V[0]*dt;
	n->X[1] += n->V[1]*dt;
}

void ensamble_motion(Neutron* n, int LEN_n, float dt){
	for (int i=0; i<LEN_n; i++){
		motion(&n[i], dt);
	}
}

void scatter(Neutron* n, float A, float T){ //elastic scattering
	float vn[2] = {n->V[0], n->V[1]};
	float vt_mag = sqrt(2*(- kB*T*log(1 - ((float)rand()/(float)RAND_MAX)))/A);
	float vt_dir = 2*M_PI*(float)rand()/(float)RAND_MAX;
	float vt[2] = {vt_mag*cos(vt_dir), vt_mag*sin(vt_dir)};
	float vcm[2] = {(vn[0] + A*vt[0]) / (A+1), (vn[1] + A*vt[1]) / (A+1)};

	float Vn[2] = {-vn[0] + 2*vcm[0], -vn[1] + 2*vcm[1]};
	n->V[0] = Vn[0];
	n->V[1] = Vn[1];
}

void ensamble_scatter(Neutron* n, int LEN_n, float A, float T){ //ensamble elastic scattering
	for (int i=0; i<LEN_n; i++){
		scatter(&n[i], A, T);
	}
}


// interaction probability
// microscopic cross-section s and number density N
// macroscopic cross-section S = sN
// probability of interaction at path of length L is P(L) = 1 - exp(-SL)
// all that's required are a list of s and N, or just a table of S

// event:
// 	absorption:
// 		Fission
// 		Capture
// 	scatter:
// 		Elastic
// 		Inelastic


int main() {
	
	Neutron n_test = {.X={0, 0}, .V={20, 0}};
	
	printN(n_test);
	
	scatter(&n_test,2,500000000);

	printN(n_test);

	printf("==================\n");

	int L_n = 5;

	Neutron ns[L_n];
	for (int i = 0; i < L_n; i++){
		ns[i].X[0] = 0;
		ns[i].X[1] = 0;
		ns[i].V[0] = 1;
		ns[i].V[1] = 1;
	}

	printNE(ns, L_n);

	ensamble_scatter(ns, L_n, 5, 500);

	printNE(ns, L_n);

	ensamble_scatter(ns, L_n, 5, 500);

	printNE(ns, L_n);
	
	printf("==================\n");
	Neutron n = {.X={0, 0}, .V={0, 20e6}};
	printN(n);
	FILE *fp = fopen("out.txt", "w");
	if (!fp) printf("NO FILE");
	fprintf(fp, "H,X,Y\n");
	fprintf(fp, "%f,%e,%e\n", mag(n.V), n.X[0], n.X[1]);
	for (int i = 0; i < 100; i++){
		scatter(&n,238,300);
		n.X[0] += n.V[0]*0.0001;
		n.X[1] += n.V[1]*0.0001;
		fprintf(fp, "%f,%e,%e\n", mag(n.V), n.X[0], n.X[1]);
	}
	fclose(fp);

	printf("==================\n");
	



	return 0;
}
