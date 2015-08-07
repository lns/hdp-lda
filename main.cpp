#include <cstdio>
#include <cstdlib>
#include "rng.hpp"
#include "hdp.hpp"
#include "pct.hpp"

int main(int argc, char* argv[])
{
	uint32_t Ndoc=0, Nword=0;
	double alpha=1, beta=0.5, gamma=1;
	char * dat = NULL;
	char * outdir = (char*)"./";
	uint32_t max_iter = 100, out_iter = 0;
	uint32_t seed = 0, verbosity = 1;

	if(1==argc) {
		fprintf(stderr," Usage: %s [OPTIONS]\n",argv[0]);
		fprintf(stderr,"	-data		Datafile in lda-c format\n");
		fprintf(stderr,"	-ndoc		Number of doc in data\n");
		fprintf(stderr,"	-nword		Number of vocabulary in data\n");
		fprintf(stderr,"	-outdir		Directory for output (./)\n");
		fprintf(stderr,"	-alpha		2nd-level effective sample size in HDP (1.0)\n");
		fprintf(stderr,"	-beta		Prior for Dirichlet dist on words (0.5)\n");
		fprintf(stderr,"	-gamma		1st-level effective sample size in HDP (1.0)\n");
		fprintf(stderr,"	-max_iter	Max iteration for CRF procedure (100)\n");
		fprintf(stderr,"	-out_iter	Output iteration for CRF procedure (max_iter)\n");
		fprintf(stderr,"	-seed		Random seed (0)\n");
		fprintf(stderr,"	-verbosity	Verbosity ranges from 0 to 2. (1)\n");
		return 0;
	} else if(0==(argc%2)) {
		fprintf(stderr,"Incorrect number of arguments %d\n.",argc);
		return -1;
	}
	for(int i=1;i<argc;i++)
	{
		if(0==strcmp(argv[i],"-data"))
			dat = argv[++i];
		else if(0==strcmp(argv[i],"-outdir"))
			outdir = argv[++i];
		else if(0==strcmp(argv[i],"-ndoc"))
			Ndoc = strtol(argv[++i],NULL,10);
		else if(0==strcmp(argv[i],"-nword"))
			Nword = strtol(argv[++i],NULL,10);
		else if(0==strcmp(argv[i],"-alpha"))
			alpha = strtod(argv[++i],NULL);
		else if(0==strcmp(argv[i],"-beta"))
			beta = strtod(argv[++i],NULL);
		else if(0==strcmp(argv[i],"-gamma"))
			gamma = strtod(argv[++i],NULL);
		else if(0==strcmp(argv[i],"-max_iter"))
			max_iter = strtol(argv[++i],NULL,10);
		else if(0==strcmp(argv[i],"-out_iter"))
			out_iter = strtol(argv[++i],NULL,10);
		else if(0==strcmp(argv[i],"-seed"))
			seed = strtol(argv[++i],NULL,10);
		else if(0==strcmp(argv[i],"-verbosity"))
			verbosity = strtol(argv[++i],NULL,10);
		else
			fprintf(stderr,"Unknown option %s\n",argv[i++]);
	}
	if(out_iter==0)
		out_iter = max_iter;

	HDP hdp(Ndoc,Nword); //#{doc}, #{vocab}
	hdp.read_data(dat);
	hdp.config(alpha,beta,gamma);
	if(verbosity>0)
	{
		fprintf(stderr,"alpha:	%8lf\n",alpha);
		fprintf(stderr,"beta:	%8lf\n",beta);
		fprintf(stderr,"gamma:	%8lf\n",gamma);
	}
	hdp.init();
	hdp.summary(verbosity);
	lcg64(seed);
	for(uint32_t i=1;i<=max_iter;i++) {
		hdp.gibbs_table();
		hdp.remove_empty();
		hdp.gibbs_menu();
		hdp.remove_empty();
		if(verbosity>0)
			printf("iter: %3u\t",i);
		if(i%out_iter==0) {
			char topic_fname[4096]; //overflow?
			char assign_fname[4096]; //overflow?
			sprintf(topic_fname,"%s/%04d_topics.txt",outdir,i);
			sprintf(assign_fname,"%s/%04d_assignments.txt",outdir,i);
			FILE* topic_f = fopen(topic_fname,"w");
			FILE* assign_f = fopen(assign_fname,"w");
			hdp.output_topics(topic_f);
			hdp.output_assignments(assign_f);
			fclose(topic_f);
			fclose(assign_f);
		}
		hdp.summary(verbosity);
	}
	return 0;
}

