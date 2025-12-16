// Header file for input output functions
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define NUM_INTER 3
#define kB 0.0000015
#define mag(x) sqrt(x[0]*x[0] + x[1]*x[1])
#define mag2(x) x[0]*x[0] + x[1]*x[1]

int genID() {
	static int num = 0; 
	return num++;      
}

enum inter {
  SCATTER,
  ABSORB,
};

typedef struct _Neutron Neutron;
struct _Neutron {
	float X[2]; 	// position 
	float V[2]; 	// velocity
	float L;	// to interaction 
	enum inter s;	// next interaction type
	int id;		// identifier for the neutron (to enable reuse of destroyed neutrons by renaming)
	bool alive;	// if false will stop tracking neutron and consider it absorbed
};

void printN(Neutron n){
	printf("Neutron:\n");
	printf("  X: %e, %e\n", n.X[0], n.X[1]);
	printf("  V: %e, %e\n", n.V[0], n.V[1]);
}

void move(Neutron* n){
	float Vmag = mag(n->V);
       	n->X[0] += n->L*n->V[0]/Vmag;
       	n->X[1] += n->L*n->V[1]/Vmag;
}

void sample_inter_dist(Neutron* n, float s[], int LEN_s){
	for (int i = 0; i < LEN_s; i++){
        	n->L = -log(1-((float) rand()) / ((float) RAND_MAX))/s[i];
		n->s = i;
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

float s_a(Neutron n, float A){
	static const float a = 38;
	static const float alfa = 0.12;
	static const float b = 0.85;
	static const float k1 = 0.07;
	static const float k2 = 0.44;
	
	float A03 = pow(A, 0.333333333);
	float K = k2*A03;
	float e = sqrt(a + b*mag2(n.V)/2)-sqrt(mag2(n.V)/2);
	float R = 1e-15*A03;
	float L = 1/mag(n.V);

	return 2*M_PI*(R+L)*(R+L)*(1 - alfa*cos(K*(e + k1*e*e)));
}

void s_temp(float* s, int LEN_S, float T){
	static const float T0 = 300;
	float rootT0_T = sqrt(T0/T);
	for (int i=0; i<LEN_S; i++) s[i] *= rootT0_T;
}


// event:
// 	absorption:
// 		Fission		cross-section s = 2 pi (R + h/mv)^2 (1 - a cosB) then sampled from sf/sa
// 		Capture		cross-section s = 2 pi (R + h/mv)^2 (1 - a cosB)
// 	scatter:
// 		Elastic		constant cross-section
// 		Inelastic	constant cross-section but maybe sampled randomly afterwards (TBD)

// TODO
// 	radius of nucleus for fast neutrons R ~ 1e-15*A^0.3333
// 	energy cross-section fit as s ~ pi*(R + h/(mv))^2 ~ pi*(1e-15*A^0.3333 + h/(mv))^2
// 	add target temp dependence as s(T) = s0 * sqrt(T0 / T)
//	figure out a lookup table of cross-sections
//	implement fission and absorption
//	Figure out nucleus temperature because it might be useful
//	data needed: A, sf/sa, nu

//Nucleus,	A,	Thermal ν, Fast ν, 	  Thermal Ratio (σf/σa), 	Fast Ratio (σf/σa)
//Thorium-232,	232.04,	—,	   2.46,	  ~0.00,		        0.07
//Uranium-233,	233.04,	2.50,	   2.65,	  0.920,		        0.94
//Uranium-235,	235.04,	2.44,	   2.61,	  0.855,		        0.82
//Uranium-238,	238.05,	—,	   2.82,	  ~0.00,		        0.15
//Plutonium-239,239.05,	2.88,	   3.12,	  0.735,		        0.90
//Plutonium-241,241.06,	2.95,	   3.14,	  0.737,		        0.91

int main() {
	
	Neutron n = {.X={0, 0}, .V={0, 2e6}, .L=0, .alive=true, .id=genID()};
	printN(n);
	FILE *fp = fopen("out.txt", "w");
	if (!fp) printf("NO FILE");
	fprintf(fp, "H,X,Y,L\n");
	float A = 238;
	float T = 300;
	int LEN_S = 2;
	float s[LEN_S];
	s[0] = 5;
	s[1] = s_a(n, A);
	fprintf(fp, "%f,%e,%e,%e\n", mag(n.V), n.X[0], n.X[1], n.L);
	for (int i = 0; i < 1000000; i++){
		sample_inter_dist(&n, s, LEN_S);
		move(&n);
		scatter(&n, A, T);
		fprintf(fp, "%f,%e,%e,%e\n", mag(n.V), n.X[0], n.X[1], n.L);
	}
	
	fclose(fp);

	printf("==================\n");
	
	printf("Neutron struct size: %ld\n", sizeof(Neutron));


	return 0;
}
