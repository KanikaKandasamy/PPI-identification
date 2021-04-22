#include "global_parameters.h"
#include "PtoHSP.h"
#include "scoreing_matrix.h"

int main(int argc, char * argv[]) {

	cout<<"-------------------------------------------------------------------\n";
	string error_msg = "To predict interactions\n";
	cout<<error_msg;
	cout<<"-------------------------------------------------------------------\n";
	for(int a = 0; a < argc; a ++){
		if(!strcmp(argv[a], "-p")){
			PROTEIN_FN = argv[a+1];
		}
		if(!strcmp(argv[a], "-h")){
			HSP_FN = argv[a+1];
		}
		if(!strcmp(argv[a], "-tr")){
			TRAIN_FN = argv[a+1];
		}
		if(!strcmp(argv[a], "-pos")){
			TEST_FN_POS = argv[a+1];
		}
		if(!strcmp(argv[a], "-neg")){
			TEST_FN_NEG = argv[a+1];
		}
		if(!strcmp(argv[a], "-o")){
			OUTPUT_FN = argv[a+1];
			OUTPUT_pos_FN = OUTPUT_FN + ".pos";
			OUTPUT_neg_FN = OUTPUT_FN + ".neg";
		}
		if(!strcmp(argv[a], "-e")){
			PROTEOME = 1;
		}
		if(!strcmp(argv[a], "-Thc")){
			T_hsp_max = atoi(argv[a+1]);
		}
	}
	cout<<"PROTEIN_FN: "<<PROTEIN_FN<<endl;
	cout<<"HSP_FN: "<<HSP_FN<<endl;
	cout<<"TRAIN_FN: "<<TRAIN_FN<<endl;
	cout<<"TEST_FN_POS: "<<TEST_FN_POS<<endl;
	cout<<"TEST_FN_NEG: "<<TEST_FN_NEG<<endl;
	cout<<"OUTPUT_FN: "<<OUTPUT_FN<<endl;
	cout<<"OUTPUT_pos_FN: "<<OUTPUT_pos_FN<<endl;
	cout<<"OUTPUT_neg_FN: "<<OUTPUT_neg_FN<<endl;
	cout<<"PROTEOME: "<<PROTEOME<<endl;

	load_protein(PROTEIN_FN);
	load_BLOSUM_convert(BLOSUM_convert);
	init_flag();
	PtoHSP hsp;
	if(DEBUG){
		cout<<"loading hsp finished, printing hsp flag"<<endl;
		print_flag();
		cout<<"printing hsp flag finished"<<endl;
	}
	SCORING_MATRIX score_matrix(hsp);
	cout<<"Loading training set"<<endl;
	score_matrix.load_traing(TRAIN_FN, hsp);
	if(!PROTEOME){
		score_matrix.load_test(TEST_FN_POS, 'p');
		score_matrix.load_test(TEST_FN_NEG, 'n');
	}	
	else{
		cout<<"printing the entire proteome score"<<endl;
		score_matrix.print_entire_final_score_matrix();
	}	
    cout<<"----end------"<<endl;
	return 0;
}
