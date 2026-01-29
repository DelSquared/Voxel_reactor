#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "nuclear.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif




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
	float s[2];
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
		printf("number of neutrons: %d, (%ld)\n", countAlive(&N), (long)N.used);
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
