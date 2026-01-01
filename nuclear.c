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
	bool caused_fis;// if true it means it caused fission and then will revert to false 
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

// basic dynammic array to manage future multiple neutrons
// inspired from https://stackoverflow.com/questions/3536153/c-dynamically-growing-array

typedef struct {
	Neutron *array;
	size_t used;
	size_t size;
} NeutronArray;

void initArray(NeutronArray* a, size_t initialSize) {
	a->array = malloc(initialSize * sizeof(Neutron));
	a->used = 0;
	a->size = initialSize;
}

void insertArray(NeutronArray* a, Neutron element) {
	// a->used is the number of used entries, because a->array[a->used++] updates a->used only *after* the array has been accessed.
	// Therefore a->used can go up to a->size 
	if (a->used == a->size) {
		a->size += 1;
		a->array = realloc(a->array, a->size * sizeof(Neutron));
	}
	a->array[a->used++] = element;
}

bool allDead(NeutronArray* a){
	bool check = true;
	for (int i = 0; i < a->used; i++){
		check = check || a->array[i].alive;
	}
	return !check;
}

int countAlive(NeutronArray* a){
	int alive = 0;
	for (int i = 0; i < a->used; i++){
		alive += (int)a->array[i].alive;
	}
	return alive;
}


void fissileNeutrons(NeutronArray* a, Material m){
	float v2, nu, ang;
	for (int i = 0; i < a->used; i++){
		v2 = mag2(a->array[i].V);
		nu = m.nu_th*( v2 < 1e6 ) + m.nu_f*( v2 >= 1e6 );
		if (a->array[i].caused_fis){
			// reset velocity to MeV (placeholder value)
			float vnew = 2e6/mag(a->array[i].V);
			a->array[i].V[0] *= vnew;
			a->array[i].V[1] *= vnew;
			a->array[i].caused_fis = false;
			nu -= 1;
			//creates new neutron (a single one for now) and adds some velocity spread
			for (float f = nu; f >=1; f-=1){
				ang = 2*M_PI*(float)rand()/(float)RAND_MAX;
				insertArray(a, a->array[i]); 
				a->array[a->used-1].V[0] = 2e6*cos(ang);
				a->array[a->used-1].V[1] = 2e6*sin(ang);
			}
			if (nu > (float)rand()/(float)RAND_MAX){
				ang = 2*M_PI*(float)rand()/(float)RAND_MAX;
				insertArray(a, a->array[i]); 
				a->array[a->used-1].V[0] = 2e6*cos(ang);
				a->array[a->used-1].V[1] = 2e6*sin(ang);
			}
		}
	}
	
}

void freeArray(NeutronArray *a) {
	free(a->array);
	a->array = NULL;
	a->used = a->size = 0;
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
	return fmax(s_t, sref_f);
}

void fission(Neutron* n, float sf_sa_th, float sf_sa_f){
	float p = sf_sa_th*(mag2(n->V) < 1e6) + sf_sa_f*(mag2(n->V) >= 1e6);// will be kept simple for now, constant in thermal and constant in fast regimes
	if (p > (float)rand()/(float)RAND_MAX){
		n->alive = true;
		n->caused_fis = true;
	}
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
	else{
		if (m.s_fis_th > 0 || m.s_fis_f > 0) fission(n, m.s_fis_th/m.s_abs_th, m.s_fis_f/m.s_abs_f);
	}

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
	for (int i = 0; i < 10000; i++){
		step(&n, m, T);
		if (n.alive) {
			fprintf(fp, "%f,%e,%e,%e,%e,%e\n", mag(n.V), n.X[0], n.X[1], n.L, s[0], s[1]);
		}
		else{
			break;
			printf("dead at: %d\n", i);
		}
	}
	
	fclose(fp);

	printf("==================\n");
	
	NeutronArray N;
	Material M = {.A=235, .s_scat=9, .s_abs_th=683, .s_abs_f=1.09, .s_fis_th=583, .s_fis_f=1, .nu_th=2.44, .nu_f=2.61};
	int LEN_N = 1;
	int MAX_N = 10;
	initArray(&N, LEN_N); //factor of 2 extra because might as well since the first fission event will trigger an array doubling event anyway
	for (short i = 0; i < LEN_N; i++){
		float ang = 2*M_PI*(float)rand()/(float)RAND_MAX;
		Neutron _n = {.X={0, 0}, .V={2e6*cos(ang), 2e6*sin(ang)}, .L=0, .alive=true, .caused_fis=false, .id=genID()};
		insertArray(&N, _n);
	}
	FILE *fp3 = fopen("out2.txt", "w");
	if (!fp3) printf("NO FILE");
	for (short i = 0; i < 2*MAX_N+3; i++){
		if (i > 0) fprintf(fp3, ",");
		fprintf(fp3, "X%d,Y%d", i, i);
	}
	fprintf(fp3, "\n");
	for (short i = 0; i < N.used; i++){
		if (i > 0) fprintf(fp3, ",");
		fprintf(fp3, "%e,%e", N.array[i].X[0],  N.array[i].X[1]);
	}
	fprintf(fp3, "\n");

	for (int j = 0; j < 100000; j++){
		printf("number of neutrons: %d, (%ld)\n", countAlive(&N), N.used);
		for (short i = 0; i < N.used; i++){
			if (N.array[i].alive) {
				if (N.array[i].caused_fis){
					printf("Fission!\n");
					fissileNeutrons(&N, M);
				}
				if (i > 0) fprintf(fp3, ",");
				fprintf(fp3, "%e,%e", N.array[i].X[0],  N.array[i].X[1]);
				step(&N.array[i], M, T);
			}
		}
		if (allDead(&N)) {
			break;
		}
		if (N.used >= MAX_N) break;
		fprintf(fp3, "\n");
	}

	
	fclose(fp3);
	printf("==================\n");


	return 0;
}
