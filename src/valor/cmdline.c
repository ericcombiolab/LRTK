#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "structural_variation.h"
#include "vector.h"
#include "valor.h"
#include "cmdline.h"


sv_type parse_svs(char * optt){

	vector_t *sv_strings = dang_string_tokenize(optt, ",");
	sv_type sv_to_find = 0;

	int i;
	for(i=0;i<sv_strings->size;i++){
		sv_to_find|=atosv(vector_get(sv_strings,i));
	}
	vector_free(sv_strings);
	return sv_to_find;
}

int parse_command_line( int argc, char** argv, parameters* params)
{
	static struct option long_options[] = 
	{
		{"input"  , required_argument,   0, 'i'},
		{"sonic"  , required_argument,   0, 's'},
		{"threads", required_argument,   0, 't'},
		{"help"   , no_argument,         0, 'h'},
		{"version", no_argument,         0, 'v'},
		{"out"    , required_argument,	 0, 'o'},
		{"svs_to_find", required_argument,0,'f'},
		{"ploidy", required_argument,0,'p'},
		{"single-copy-chr", required_argument,0,'y'},
		{"low_mem", no_argument,	 0, 'm'},
		{"log_file", required_argument,	 0, 'l'},
		{"contig_count", required_argument,	 0, 'c'},
		{0        , 0,                   0,  0 }
	};
  
	if( argc == 1)
	{
		print_help();
		return -1;
	}
 	int o;
	int index;
	vector_t *haplotype_chrs = vector_init(sizeof(char)*128,128);
	char *buffer;
	while( ( o = getopt_long( argc, argv, "i:s:t:hvo:f:p:y:ml:c:", long_options, &index)) != -1)
	{
		switch( o)
		{
			case 'i':
				set_str( &( params->bam_file), optarg);
			break;
		
			case 'l':
				set_str( &( params->logfile), optarg);
			break;
		
			
			case 'p':
				params->ploidy= atoi(optarg);
			break;
			case 'y':
				buffer = malloc(128);
				strncpy(buffer,optarg,127);
				vector_soft_put(haplotype_chrs,buffer);
			break;
			case 'f':
				params->svs_to_find = parse_svs(optarg);
			break;

			case 's':
				set_str( &( params->sonic_file), optarg);
			break;
			case 'm':
				params->low_mem = 1;
			break;
			
			case 'o':
				set_str( &( params->outprefix), optarg);
			break;
			
			case 't':
				params->threads = atoi( optarg);
			break;
			case 'h':
				print_help();
				return -1;
			break;
			case 'c':
				params->chromosome_count = atoi(optarg);
			break;
			case 'v':
				fprintf( stderr, "\nVALOR: VAriation using LOng Range information.\n");
				fprintf( stderr, "Version %s\n\tLast update: %s, build date: %s\n\n", VALOR_VERSION, VALOR_UPDATE, BUILD_DATE);
				return -1;
			break; 
		}
	}
 	int ret = RETURN_SUCCESS;
	/* TODO: check parameter validity */
  
	/* check algorithms to run; run_all is the default */
	/* check if outprefix is given */
	if( params->outprefix == NULL)
	{
		fprintf( stderr, "[VALOR CMDLINE ERROR] Please enter the output directory name prefix using the --out option.\n");
		ret |=RETURN_ERROR;
	}
	
	if( params->svs_to_find == 0)
	{
		fprintf( stderr, "[VALOR CMDLINE ERROR] Please enter type of SV's to discover -f or --svs_to_find [DUP,IDUP,INV].\n");
		ret |=RETURN_ERROR;
	}
	if( params->bam_file ==NULL){
		fprintf( stderr, "[VALOR CMDLINE ERROR] Please enter the input BAM file path using the -i or --input option.\n");
		ret |=RETURN_ERROR;
	}

	if( params->sonic_file == NULL){
		fprintf( stderr, "[VALOR CMDLINE ERROR] Please enter the sonic file path using the -s or --sonic option.\n");
		ret|= RETURN_ERROR;
	
	}
	if( params->logfile == NULL){
	        char *tmp_logfilename = (char *) malloc(sizeof(char *) * (strlen(VALOR_DEFAULT_LOG_FILE)+strlen(params->outprefix)+2));
		sprintf( tmp_logfilename, "%s-%s", params->outprefix, VALOR_DEFAULT_LOG_FILE);
		set_str( &( params->logfile), tmp_logfilename);
		free( tmp_logfilename);
	}
	int i;
	sonic *snc = sonic_load(params->sonic_file);
	params->chr_copy_count = malloc(sizeof(int) * params->chromosome_count);
	for(i=0;i<params->chromosome_count;i++){
		params->chr_copy_count[i] = params->ploidy;
	}
	int j;
	for(i=0;i<haplotype_chrs->size;i++){
		for(j=0;j<snc->number_of_chromosomes;j++){
			if(strncmp(vector_get(haplotype_chrs,i),snc->chromosome_names[j],128)  ==0){
				params->chr_copy_count[j] =1;
				break;
			}
		}
		if(j==snc->number_of_chromosomes){
			fprintf(stderr,"contig %s does not exist in the reference",(char *)vector_get(haplotype_chrs,i));
			ret |= RETURN_ERROR;
		}
	}
	return ret;
}

void print_help( void)
{  
	fprintf( stdout, "\nVALOR: VAriation using LOng Range information.\n");

	fprintf( stdout, "Version %s\n\tLast update: %s, build date: %s\n\n", VALOR_VERSION, VALOR_UPDATE, BUILD_DATE);	
	fprintf( stdout, "Required Parameters:\n");

	fprintf( stdout, "\t-i, --input [BAM files]        : Input files in sorted BAM format.\n");
	fprintf( stdout, "\t-o, --out   [output folder]    : Folder to put stuff in\n");

	fprintf( stdout, "\t-s, --sonic  [sonic file]      : Sonic file. Check: https://github.com/calkan/sonic.\n");

	fprintf( stdout, "\t-p, --ploidy  [Number of chromosome sets]      : Default is 2.\n");
	fprintf( stdout, "\t-y, --single-copy-chr  [Chromosome name]      : sets copy number of a chromosome to 1 (example -y X -y Y for male humans).\n");

	fprintf( stdout, "\t-f, --svs_to_find   [sv type]: Comma separated list of SV types (i.e. DUP,IDUP,INV).\n");
	fprintf( stdout, "Optional Parameters:\n");

	fprintf( stdout, "\t-c, --contig_count   [Number of contigs to run VALOR (First N contigs in sorted BAM) ]: default is 24 (Male Human)\n");
	fprintf( stdout, "\t-t, --threads   [Number of threads to run VALOR]: default is 1\n");
	fprintf( stdout, "\t-m, --low_mem	[Use disk to reduce memory usage]: Not Implemented\n");
	fprintf( stdout, "\t-l, --log_file [logfile name]: default is valor.log\n");
	fprintf( stdout, "Help Parameters:\n");
	fprintf( stdout, "\t-v, --version                  : Print version and exit.\n");
	fprintf( stdout, "\t-h, --help                     : Print this help screen and exit.\n\n");

}

