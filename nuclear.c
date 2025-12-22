// Header file for input output functions
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define NUM_INTER 3
#define kB 0.0000015
#define mag(x) sqrt(x[0]*x[0] + x[1]*x[1])
#define mag2(x) x[0]*x[0] + x[1]*x[1]
#define lgst(x) 1/(1+exp(-x))

int genID() {
	static int num = 0; 
	return num++;      
}

enum inter {
  SCATTER,
  ABSORB,
};

typedef struct {
	float X[2]; 	// position 
	float V[2]; 	// velocity
	float L;	// to interaction 
	enum inter s;	// next interaction type
	int id;		// identifier for the neutron (to enable reuse of destroyed neutrons by renaming)
	bool alive;	// if false will stop tracking neutron and consider it absorbed
} Neutron;

typedef struct {
	float A; 	// atomic mass
	float N;	// atomic density
	float s_scat;	// scattering cross-section
	float s_abs_th;	// thermal absorption cross-section
	float s_abs_f;	// fast absorption cross-section
	float s_fis_th;	// thermal fission cross-section
	float s_fis_f;	// fast fission cross-section
	float nu_th;	// thermal neutron production
	float nu_f;	// fast neutron production
} Material;

void printN(Neutron n){
	printf("Neutron:\n");
	printf("  X: %e, %e\n", n.X[0], n.X[1]);
	printf("  V: %e, %e\n", n.V[0], n.V[1]);
}

void move(Neutron* n){
	if (n->alive){
		float Vmag = mag(n->V);
       		n->X[0] += n->L*n->V[0]/Vmag;
       		n->X[1] += n->L*n->V[1]/Vmag;
	}
	if (n->s == ABSORB){
		n->alive = false;
	}
}

void sample_inter_dist(Neutron* n, float s[], int LEN_s){
	float L = 1000000000;
	float sample;
	int reaction_type;
	for (int i = 0; i < LEN_s; i++){
        	sample = -log(1-((float) rand()) / ((float) RAND_MAX))/s[i];
		if (sample < L) {
			L = sample;
			reaction_type = i;
		}
        }
	n->L = L;
	n->s = reaction_type;
		//printf("%d %f\n", i, n->L);
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

float s_a(Neutron n, float A, float sref_t,  float sref_f){
	float s_t = sref_t*0.158114/mag(n.V);
	//float a = lgst(0.1*(mag(n.V) - sref_t*0.158114/sref_f));
	//return s_t*(1-a) + sref_f*(a);
	return fmax(s_t, sref_f);
}

void s_temp(float* s, int LEN_S, float T){
	static const float T0 = 300;
	float rootT0_T = sqrt(T0/T);
	for (int i=0; i<LEN_S; i++) s[i] *= rootT0_T;
}

void step(Neutron* n, Material m, float T){
	// calculate s_abs
	float s[2];
	s[SCATTER] = m.s_scat;
	s[ABSORB] = s_a(*n, m.A, m.s_abs_th, m.s_abs_f);
	// sample s
	sample_inter_dist(n, s, 2);
	// move
	move(n);
	// if not absorbed scatter
	if (n->alive) {
		scatter(n, m.A, T);
	}
	// else sample for fission event

}

// event:
// 	absorption:
// 		Fission		cross-section s = 2 pi (R + h/mv)^2 (1 - a cosB) then sampled from sf/sa
// 		Capture		cross-section s = 2 pi (R + h/mv)^2 (1 - a cosB)
// 	scatter:
// 		Elastic		constant cross-section
// 		Inelastic	constant cross-section but maybe sampled randomly afterwards (TBD)

// TODO
// 	add target temp dependence as s(T) = s0 * sqrt(T0 / T)
//	figure out a lookup table of cross-sections
//	implement fission and absorption
//	Figure out nucleus temperature because it might be useful
//	data needed: A, sf/sa, nu

//Nucleus,	A,	Thermal ν, Fast ν, Thermal Ratio (σf/σa), Fast Ratio (σf/σa)
//Thorium-232,	232.04,	—,	   2.46,   0.000,		  0.07
//Uranium-233,	233.04,	2.50,	   2.65,   0.920,		  0.94
//Uranium-235,	235.04,	2.44,	   2.61,   0.855,		  0.82
//Uranium-238,	238.05,	—,	   2.82,   0.000,		  0.15
//Plutonium-239,239.05,	2.88,	   3.12,   0.735,		  0.90
//Plutonium-241,241.06,	2.95,	   3.14,   0.737,		  0.91

int main() {
	
	Neutron n = {.X={0, 0}, .V={0, 2e6}, .L=0, .alive=true, .id=genID()};
	Material m = {.A=238, .s_scat=9, .s_abs_th=2, .s_abs_f=0.07};
	printN(n);
	FILE *fp = fopen("out.txt", "w");
	if (!fp) printf("NO FILE");
	fprintf(fp, "H,X,Y,L,s_s,s_a\n");
	float A = 238;
	float T = 300;
	int LEN_S = 2;
	float s[LEN_S];
	s[SCATTER] = 9;
	s[ABSORB] = s_a(n, A, 2, 0.07);
	printf("%f, %f\n", s[0], s[1]);
	fprintf(fp, "%f,%e,%e,%e,%e,%e\n", mag(n.V), n.X[0], n.X[1], n.L, s[0], s[1]);
	for (int i = 0; i < 1000000; i++){
		s[SCATTER] = 9;
		s[ABSORB] = s_a(n, A, 2, 0.07);
		sample_inter_dist(&n, s, LEN_S);
		move(&n);
		if (n.alive) {
			fprintf(fp, "%f,%e,%e,%e,%e,%e\n", mag(n.V), n.X[0], n.X[1], n.L, s[0], s[1]);
			scatter(&n, A, T);
		}
		else{
			break;
			printf("dead at: %d\n", i);
		}
	}
	
	fclose(fp);

	printf("==================\n");

	// to plot the absorption cross-section spectrum
	FILE *fp2 = fopen("abs.txt", "w");
	fprintf(fp2, "v,s\n");
	for (int i = 0; i < 10000; i++){
		n.V[0] = 0;
		n.V[1] = (float)i*0.01;
		fprintf(fp2, "%f,%e\n", (float)i*0.01, s_a(n, 235, 99, 0.09));
	}
	for (int i = 1; i < 10000000; i++){
		n.V[0] = 0;
		n.V[1] = (float)i*100;
		fprintf(fp2, "%d,%e\n", i*100, s_a(n, 235, 99, 0.09));
	}
	fclose(fp2);

	printf("Neutron struct size: %ld\n", sizeof(Neutron));


	return 0;
}
