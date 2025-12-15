// Header file for input output functions
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define kB 0.0000015
#define mag(x) sqrt(x[0]*x[0] + x[1]*x[1])

typedef struct _Neutron Neutron;
struct _Neutron {
	float X[2]; 	// position 
	float V[2]; 	// velocity
	float L;	// to interaction 
};

void printN(Neutron n){
	printf("Neutron:\n");
	printf("  X: %e, %e\n", n.X[0], n.X[1]);
	printf("  V: %e, %e\n", n.V[0], n.V[1]);
}

void printNE(Neutron n[], int LEN_n){
	printf("Ensemble:\n");
	for (int i = 0; i < LEN_n; i++){
		printf("  ");
		printN(n[i]);
	}
}

void move(Neutron* n){
	float Vmag = mag(n->V);
       	n->X[0] += n->L*n->V[0]/Vmag;
       	n->X[1] += n->L*n->V[1]/Vmag;
}

void sample_inter_dist(Neutron* n, float s[], int LEN_s){
	for (int i = 0; i < LEN_s; i++){
        	n->L = -log(1-((float) rand()) / ((float) RAND_MAX))/s[i];
        }
	printf("%e\n",n->L);
}

void scatter(Neutron* n, float A, float T){ //elastic scattering
	float vn[2] = {n->V[0], n->V[1]};
	float vt_mag = sqrt(2*(- kB*T*log(1 - ((float)rand()/(float)RAND_MAX)))/A);
	float vt_dir = 2*M_PI*(float)rand()/(float)RAND_MAX;
	float vt[2] = {vt_mag*cos(vt_dir), vt_mag*sin(vt_dir)};
	float vcm[2] = {(vn[0] + A*vt[0]) / (A+1), (vn[1] + A*vt[1]) / (A+1)};
	printf("x %f, y %f\n", n->V[0], n->V[1]);
	float Vn[2] = {-vn[0] + 2*vcm[0], -vn[1] + 2*vcm[1]};
	n->V[0] = Vn[0];
	n->V[1] = Vn[1];
}

void ensemble_scatter(Neutron* n, int LEN_n, float A, float T){ //ensemble elastic scattering
	for (int i=0; i<LEN_n; i++){
		scatter(&n[i], A, T);
	}
}

// event:
// 	absorption:
// 		Fission
// 		Capture
// 	scatter:
// 		Elastic
// 		Inelastic

// TODO
// 	radius of nucleus for fast neutrons R ~ 1e-15*A^0.3333
// 	energy cross-section fit as s ~ pi*(R + h/(mv))^2
// 	add target temp dependence as s(T) = s0 * sqrt(T0 / T)
//	figure out a lookup table of cross-sections
//	implement fission and absorption

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

	ensemble_scatter(ns, L_n, 5, 500);

	printNE(ns, L_n);

	ensemble_scatter(ns, L_n, 5, 500);

	printNE(ns, L_n);
	
	printf("==================\n");
	
	Neutron n = {.X={0, 0}, .V={0, 2e6}, .L=0};
	printN(n);
	FILE *fp = fopen("out.txt", "w");
	if (!fp) printf("NO FILE");
	fprintf(fp, "H,X,Y,L\n");
	float tmass = 238;
	float T = 300;
	float s_scat[1] = {5};
	fprintf(fp, "%f,%e,%e,%e\n", mag(n.V), n.X[0], n.X[1], n.L);
	for (int i = 0; i < 1000000; i++){
		sample_inter_dist(&n, s_scat, 1);
		move(&n);
		scatter(&n, tmass, T);
		fprintf(fp, "%f,%e,%e,%e\n", mag(n.V), n.X[0], n.X[1], n.L);
	}
	
	fclose(fp);

	printf("==================\n");
	



	return 0;
}
