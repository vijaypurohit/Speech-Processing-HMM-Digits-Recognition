
/**************************************************************************************
	To Read Lambda(A,B,Pi) and Observation Sequence from File
	Input: WordNames Array Index
**************************************************************************************/
void readLambdaABPi(int d, unsigned int model_type){
	FILE *fp_ind;
	char completePathInd[200]; 
	char completePathInput[200];
	char model_type_name[20];
	char lambda_folder_name[50];

	if(model_type==1) // use old model to read from input folder
	{
		sprintf(completePathInput, "%s%s %s/", input_folder, WordFolderName, WordNames[d]);
		sprintf(lambda_folder_name, "%s %s/", WordFolderName, WordNames[d]);
		strcpy(model_type_name, "Input Folder Model");
	}
	else if(model_type==2) //use new converged model to read from output folder
	{
		sprintf(completePathInput, "%s%s/%s/", output_folder, output_folder_Model_name, WordNames[d]); 
		sprintf(lambda_folder_name, "%s/%s/", output_folder_Model_name, WordNames[d]);
		strcpy(model_type_name, "Output Folder Model");
	}
	
	//printf("\n ----------------------- ----------------------- > Reading Lambda from File: %s %s (%s) < ----------------------- \n", WordFolderName, WordNames[d], model_type_name); 
	//fprintf(fp_console, "\n ----------------------- ----------------------- > Reading Lambda from File: %s %s (%s) < ----------------------- \n", WordFolderName, WordNames[d], model_type_name);				
		
	//fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Reading Lambda(A): %s/%s%s.txt\n\n", lambda_folder_name,LambdaFileNames[0], WordNames[d]);
		
		sprintf(completePathInd, "%s%s%s.txt", completePathInput,LambdaFileNames[0], WordNames[d]);  // {completePathInput}+{A_}+{1}+".txt 
		fp_ind = fopen(completePathInd, "r");				//to read input A from lamda
			if(fp_ind == NULL){ 
				perror("\nError: ");
				printf("\n File Name is: %s\n", completePathInd);
				getch();
				return;
			}	
		for (int si = 0; si < N; ++si){
			//fprintf(fp_console, "S[%d]\t", si+1);
			for (int sj = 0; sj < N; ++sj){
					fscanf(fp_ind,"%Lf",&A[si][sj]);
					//fprintf(fp_console, "%.16g(%d)\t",A[si][sj], sj+1);
				}
			//fprintf(fp_console, "\n");
		}
		//for (int i = 0; i < N; ++i)
		//{
		//	for (int j = 0; j < N; ++j)
		//		printf("%0.20g\t", A[i][j]);
		//	printf("\n");
		//}
		fflush(fp_ind);  fclose(fp_ind);


	//fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Reading Lambda(B): %s/%s%s.txt\n", lambda_folder_name,LambdaFileNames[1], WordNames[d]);
		sprintf(completePathInd, "%s%s%s.txt", completePathInput,LambdaFileNames[1], WordNames[d]);  // {completePathInput}+{B_}+{1}+".txt 
		fp_ind = fopen(completePathInd, "r");				//to read input B from lamda
			if(fp_ind == NULL){ 
				perror("\nError: ");
				printf("\n File Name is: %s\n", completePathInd);
				getch();
				return;
			}	
		
		for (int si = 0; si < N; ++si){
			//fprintf(fp_console, "S[%d]\t", si+1);
			for (int m = 0; m < M; ++m){
					fscanf(fp_ind,"%Lf",&B[si][m]);
					//fprintf(fp_console, "%.16g(%d)\t",B[si][m], m+1);
				}
			//fprintf(fp_console, "\n");
		}
		//for (int i = 0; i < N; ++i)
		//{
		//	for (int j = 0; j < M; ++j)
		//		printf("%.20g\t", B[i][j]);
		//	printf("\n");
		//}	
		fflush(fp_ind);  fclose(fp_ind); 
				
	//fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Reading Lambda(PI): %s/%s%s.txt\n",lambda_folder_name, LambdaFileNames[2], WordNames[d]);
		sprintf(completePathInd, "%s%s%s.txt", completePathInput,LambdaFileNames[2], WordNames[d]);  // {completePathInput}+{Pi_}+{1}+".txt 
		fp_ind = fopen(completePathInd, "r");				//to read input Pi from lamda	
			if(fp_ind == NULL){ 
				perror("\nError: ");
				printf("\n File Name is: %s\n", completePathInd);
				getch();
				return;
			}	
			
		for (int si = 0; si < N; ++si){
			fscanf(fp_ind,"%Lf",&PI[si]);
			//fprintf(fp_console, "%0.16g(S%d)\t",PI[si], si+1);
		}
		//fprintf(fp_console, "\n");

		fflush(fp_ind); fclose(fp_ind);
			
		//printf("\n -----------------------> Reading Completed: %s %s <----------------------- \n", WordFolderName, WordNames[d]); 
}//readLambdaABPi

void readObsSeq(int d, unsigned short int seq_type){
	
	unsigned short NumOfSequences;
	char lambda_obs_seq_file_name[50];
	char seq_type_name[20];

	FILE *fp_obs ;			//to read observation seq
	char completePathObs[200] ;

	if(seq_type == 1)
	{	
		NumOfSequences = On;
		strcpy(lambda_obs_seq_file_name, LambdaFileNames[3]);
		strcpy(seq_type_name, "TRAINING");
	}
	else if(seq_type == 2)
	{	
		NumOfSequences = Ot;
		strcpy(lambda_obs_seq_file_name, LambdaFileNames[4]);
		strcpy(seq_type_name, "TESTING");
	}

	
	sprintf(completePathObs, "%s%s %s/%s%s.txt", input_folder, WordFolderName, WordNames[d], lambda_obs_seq_file_name, WordNames[d]);  // {input_lamda/}+{WordFolderName}+" "+{1}+"/"+{obs_seq__}+{1}+".txt 

	fp_obs = fopen(completePathObs, "r"); //to save compelete observation sequence in one file
	//printf("\n File Names is :  %s \n", completePathObs);
	if(fp_obs == NULL ){ 
		perror("\n Error: ");
		printf("\n No Observation Sequence: File Names is :  %s \n", completePathObs);
		printf("\n\t Therfore Generating Sequence Now: \n");
		sequence_generation(seq_type);
		getch();
		//return;
	}
		
	printf("\n ---#---#---#---#---#---#---#---> Reading %s Observation Sequence From File: %s %s/%s%s.txt <---#---#---#---#---#---#---#--- \n", seq_type_name, WordFolderName, WordNames[d], lambda_obs_seq_file_name, WordNames[d]); 
			fprintf(fp_console, "\n ---#---#---#---#---#---#---#---> Reading %s Observation Sequence From File: %s %s/%s%s.txt <---#---#---#---#---#---#---#--- \n", seq_type_name, WordFolderName, WordNames[d], lambda_obs_seq_file_name, WordNames[d]);

		fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> %s Obs Seq: %s %s/%s%s.txt \n", seq_type_name, WordFolderName, WordNames[d],lambda_obs_seq_file_name, WordNames[d]);
				char skipDashes[1024];
				 int dd=0;
				  int num_frames_in_seq=0;
				for (int i = 0; i < NumOfSequences; ++i){

					fscanf(fp_obs,"%s %d:%d %s",&skipDashes,&dd, &num_frames_in_seq, &skipDashes);
					OFmax[i] = num_frames_in_seq;	// store max frames of each sequence

					fprintf(fp_console, "O[%d]:%i:\t", i+1,OFmax[i]);

					for (int t = 0; t < num_frames_in_seq; ++t){
							fscanf(fp_obs,"%d",&O[i][t]);
							fprintf(fp_console, "%d (%d)\t",O[i][t], t+1);
						}
					fprintf(fp_console, "\n");

					
				}// for each sequences
				
		printf("\n -----------------------> Reading Completed %s Observation Sequence: %s %s/%s%s.txt <----------------------- \n", seq_type_name, WordFolderName, WordNames[d], lambda_obs_seq_file_name, WordNames[d]); 
		printf("\n-----------------------------------------------------\n");
}//readTrainingTestingObsSeq

void readObsSeqOfLiveRecording(){

	unsigned short NumOfSequences=1;

	FILE *fp_obs ;			//to read observation seq
	char completePathObs[200] ;

	sprintf(completePathObs, "%s%s_obs_seq_.txt", input_live_voice_folder, liveRecordingFileName); 

	fp_obs = fopen(completePathObs, "r"); //to save compelete observation sequence in one file
	//printf("\n File Names is :  %s \n", completePathObs);
	if(fp_obs == NULL ){ 
		perror("\n Error: ");
		printf("\n No Observation Sequence: File Names is :  %s \n", completePathObs);
		printf("\n\t Therfore Generating Sequence Now: \n");
		live_sequence_generation();
		getch();
		//return;
	}
		
	printf("\n ---~---> Reading %s Observation Sequence From File: %s <---~--- \n", "LIVE VOICE RECORDING", completePathObs); 
		fprintf(fp_console, "\n ---~---> Reading %s Observation Sequence From File: %s <---~--- \n", "LIVE VOICE RECORDING", completePathObs);
		char skipDashes[1024];
			int dd=0;
			int num_frames_in_seq=0;
		for (int i = 0; i < NumOfSequences; ++i){

			fscanf(fp_obs,"%s %d:%d %s",&skipDashes,&dd, &num_frames_in_seq, &skipDashes);
			OFmax[i] = num_frames_in_seq;	// store max frames of each sequence

			fprintf(fp_console, "O[%d]:%i:\t", i+1,OFmax[i]);

			for (int t = 0; t < num_frames_in_seq; ++t){
					fscanf(fp_obs,"%d",&O[i][t]);
					fprintf(fp_console, "%d (%d)\t",O[i][t], t+1);
				}
			fprintf(fp_console, "\n");
		}// for each sequences
				
		//printf("\n ---~---~---~---~---~---> Reading Completed %s Observation Sequence From File: %s <---~---~---~---~---~--- \n", "LIVE VOICE RECORDING", completePathObs); 
		printf("\n-----------------------------------------------------\n");
}//readTrainingTestingObsSeq
/**************************************************************************************
	P1 Evaluation Problem (Scoring Problem) | SOLUTION: FORWARD PROCEDURE
	Input: Observation Sequnce, Observation Count
	Output: Probability(Observation sequence given lambda)
**************************************************************************************/
long double P1_Forward_Procedure(int *Oi, int o){
	long double alpha_t_sum =0;
	long double prob_alpha=0;

	int T = OFmax[o];
	//S1: Intialization //t=0 (or t1)
	for(int s=0; s<N; ++s)
		alpha[0][s] = PI[s]*B[s][Oi[0]-1];

	//S2: Induction
	for(int t=0; t<T-1; ++t){
		for(int sj=0; sj<N; ++sj){
				alpha_t_sum=0;
			for(int si=0; si<N; ++si){
					alpha_t_sum+=alpha[t][si]*A[si][sj];
			}// for si
			 alpha[t+1][sj] = alpha_t_sum*B[sj][Oi[t+1]-1];
		}//for sj
	}//for t

	//S3: Termination
	for(int s=0; s<N; ++s)
		 prob_alpha += alpha[T-1][s];

	//printing in file
	//fprintf(fp_console, "\n -----------------------> ALPHA MATRIX for ~~~~~~~~~~~~~ O[%d]\n", o+1);
	//for (int s = 0; s < N; ++s)
	//{
		//fprintf(fp_console, "S[%d]:: \t", s+1);
		//for (int t = 0; t < T; ++t)
		//{
			//fprintf(fp_console, "%.5g(%d)\t",alpha[t][s],t+1);
		//}
		//fprintf(fp_console, "\n");
	//}
	//fprintf(fp_console, "\n :::: Alpha Probabiltiy :  %g\n", prob_alpha);

	probability_alpha = prob_alpha;
	return prob_alpha;
}//P1_Forward_Procedure

void PRINT_P1_Forward_Procedure(int o){
	fprintf(fp_console, "\n -----------------------> ALPHA MATRIX for ~~~~~~~~~~~~~ O[%d]\n", o+1);
	int T = OFmax[o];
	for (int s = 0; s < N; ++s)
	{
		fprintf(fp_console, "S[%d]:: \t", s+1);
		for (int t = 0; t < T; ++t)
		{
			fprintf(fp_console, "%.5g(%d)\t",alpha[t][s],t+1);
		}
		fprintf(fp_console, "\n");
	}
	fprintf(fp_console, "\n :::: Alpha Probabiltiy :  %g\n", probability_alpha);

}//PRINT_P1_Forward_Procedure

/**************************************************************************************
	P1 Evaluation Problem (Scoring Problem) | SOLUTION: BACKWARD PROCEDURE
	Input: Observation Sequnce, Observation Count
	Output: Probability(Observation sequence given lambda)
**************************************************************************************/
long double P1_Backward_Procedure(int *Oi, int o){
	long double beta_t_sum =0;
	long double prob_beta=0;
	int T = OFmax[o];
	//S1: Intialization //t=0 (or t1)
	for(int s=0; s<N; ++s)
		beta[T-1][s] = 1;
	
	//S2: Induction
	for(int t=T-2; t>=0; --t){
		for(int si=0; si<N; ++si){
				beta_t_sum=0;
			for(int sj=0; sj<N; ++sj){
					beta_t_sum += A[si][sj]*B[sj][Oi[t+1]-1]*beta[t+1][sj];
			}// for sj
			 beta[t][si] = beta_t_sum;
		}//for si
	}//for t

	//S3: Termination
	for(int s=0; s<N; ++s)
		 prob_beta += beta[0][s];

	//printing in file
	//fprintf(fp_console, "\n -----------------------> BETA MATRIX for ~~~~~~~~~~~~~ O[%d]\n", o+1);
	//for (int s = 0; s < N; ++s)
	//{
	//	fprintf(fp_console, "S[%d]:: \t", s+1);
	//	for (int t = 0; t < T; ++t)
	//	{
	//		fprintf(fp_console, "%.5g(%d)\t",beta[t][s],t+1);
	//	}
	//	fprintf(fp_console, "\n");
	//}

	//fprintf(fp_console, "\n :::: Beta Probabiltiy :  %g\n", prob_beta);
	probability_beta=prob_beta;
	return prob_beta;
}//P1_Backward_Procedure

void PRINT_P1_Backward_Procedure(int o){
	fprintf(fp_console, "\n -----------------------> BETA MATRIX for ~~~~~~~~~~~~~ O[%d]\n", o+1);
	int T = OFmax[o];
	for (int s = 0; s < N; ++s)
	{
		fprintf(fp_console, "S[%d]:: \t", s+1);
		for (int t = 0; t < T; ++t)
		{
			fprintf(fp_console, "%.5g(%d)\t",beta[t][s],t+1);
		}
		fprintf(fp_console, "\n");
	}

	fprintf(fp_console, "\n :::: Beta Probabiltiy :  %g\n", probability_beta);

}//PRINT_P1_Backward_Procedure

/**************************************************************************************
	P2 Uncovering the Problem | SOLUTION: VITERBI ALGO
	Input: Observation Sequnce, Observation Count
**************************************************************************************/
void P2_Viterbi_Algo(int *Oi, int o){
	
	//long double temp_delta=0, temp_delta_max=0;
	int argmax=0; //argument which gives the max probability of delta
	int T = OFmax[o];
	//S1: Intialization //t=0 (or t1)
	for(int s=0; s<N; ++s){
		delta[0][s] = PI[s]*B[s][Oi[0]-1];
		psi[0][s]=-1;
	}

	//S2: Recursion
	for(int t=1; t<T; ++t){
		
		for(int sj=0; sj<N; ++sj){
				argmax=0;				// argument i for which delta[t][j] is maximum, let first state is max
			for(int si=1; si<N; ++si){
				//temp_delta = delta[t-1][si]*A[si][sj];
				//temp_delta_max = delta[t-1][argmax]*A[argmax][sj];
					if((delta[t-1][si]*A[si][sj]) > (delta[t-1][argmax]*A[argmax][sj]))
						argmax=si;
			}// for si
			 delta[t][sj] = delta[t-1][argmax]*A[argmax][sj]*B[sj][Oi[t]-1];
			 psi[t][sj]=argmax;		//argument index which gave the max probability.
		}//for sj
	
	}//for t
	
	//S3: Termination
	argmax=0;
	for(int sj=1; sj<N; ++sj){
			if(delta[T-1][sj] > delta[T-1][argmax])
				argmax=sj;
	}//for sj
	Pstar = delta[T-1][argmax];
	Qstar[T-1]=argmax;		//argument index which gave the max probability.

	//S4: Backtracking, State Sequence Path.
	for(int t=T-2; t>=0; --t)
	{
		Qstar[t]=psi[t+1][Qstar[t+1]];  
	}//for t

	//fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Value of Pstar: %g\n", Pstar);
	//fprintf(fp_console, ":::::O[%d]:\t",o+1);
	//for (int t = 0; t < T; ++t){
	//		fprintf(fp_console, "%d (%d)\t",Oi[t], t+1);
	//}	
	//fprintf(fp_console, "\n:::::Q[%d]:\t",o+1);
	//for (int t = 0; t < T; ++t){
	//		fprintf(fp_console, "%d (%d)\t",Qstar[t]+1, t+1);
	//}
	//fprintf(fp_console, "\n");	

}//P2_Viterbi_Algo

void PRINT_P2_Viterbi_Algo(int o){

	int T = OFmax[o];

	//printing in file
	fprintf(fp_console, "\n -----------------------> DELTA & PSI MATRIX ~~~~~~~~~~~~~ O[%d]\n", o+1);
	for (int s = 0; s < N; ++s)
	{
		fprintf(fp_console, "S[%d]:: \t", s+1);
		for (int t = 0; t < T; ++t)
		{
			fprintf(fp_console, "%.5g(%d)\t",delta[t][s],t+1);
		}
		fprintf(fp_console, "\n");
	}
	fprintf(fp_console, "\n");
	for (int s = 0; s < N; ++s)
	{
		fprintf(fp_console, "Psi[%d]\t", s+1);
		for (int t = 0; t < T; ++t)
		{
			fprintf(fp_console, "%d(%d)\t",psi[t][s]+1,t+1);
		}
		fprintf(fp_console, "\n");
	}

}//PRINT_P2_Viterbi_Algo

void PRINT_P2_PStar_State_Sequence(int *Oi, int o){
	int T = OFmax[o];
	fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Value of Pstar: %g\n", Pstar);
	fprintf(fp_console, ":::::O[%d]:\t",o+1);
	for (int t = 0; t < T; ++t){
			fprintf(fp_console, "%d (%d)\t",Oi[t], t+1);
	}	
	fprintf(fp_console, "\n:::::Q[%d]:\t",o+1);
	for (int t = 0; t < T; ++t){
			fprintf(fp_console, "%d (%d)\t",Qstar[t]+1, t+1);
	}
	fprintf(fp_console, "\n");
}//PRINT_P2_PStar_State_Sequence

/**************************************************************************************************
	P2 Uncovering the Problem | SOLUTION: GAMMA PRCOEDURE 
	Probability of being in state Si at time t given the observation sequnece and model lambda
**************************************************************************************************/
long double P2_Gamma_Procedure(int o){

	int T = OFmax[o];

	int Q[Tmax];	 // Gamma Procedure, Best State Sequence Path which is most likely state for each time point t
	int argmax=0;
	long double denominator_sum=0;
	long double prob_gamma=1;

	for(int t=0; t<T; t++){

		denominator_sum=0;
		for(int si=0; si<N; si++)
			denominator_sum += alpha[t][si]*beta[t][si]; 
		
		argmax=0;
		for(int si=0; si<N; si++){
			gamma[t][si] = alpha[t][si]*beta[t][si]/denominator_sum;
			if(gamma[t][si]>gamma[t][argmax])
				argmax=si;
		}//for si
		Q[t]=argmax;
	}// for t

	for(int t=0; t<T; t++)
		prob_gamma*=gamma[t][Q[t]];

	return prob_gamma;
}//P2_Gamma_Procedure

void PRINT_P2_Gamma_Procedure(int o){
	int T = OFmax[o];
	//printing in file
	fprintf(fp_console, "\n -----------------------> GAMMA MATRIX ~~~~~~~~~~~~~ O[%d]\n", o+1);
	for (int s = 0; s < N; ++s)
	{
		fprintf(fp_console, "S[%d]:: \t", s+1);
		for (int t = 0; t < T; ++t)
		{
			fprintf(fp_console, "%.5g(%d)\t",gamma[t][s],t+1);
		}
		fprintf(fp_console, "\n");
	}
	fprintf(fp_console, "\n");

}//PRINT_P2_Gamma_Procedure

/**************************************************************************************************
	P3 Re-estimation Problem | SOLUTION: Baum Welch Method
	XI: Prob of being in state Si at time t and state Sj at time t+1 given the model lambda
		and observation sequence.
	Input: Observation Sequnce, Observation Count
**************************************************************************************************/
void P3_Baum_Welch_Procedure(int *Oi, int o){
	int T = OFmax[o];
	long double denominator_sum=0; 

	for(int t=0; t<T-1; t++) // For all T-1 
	{
		denominator_sum=0;
		for(int si=0; si<N; si++){
			for(int sj=0; sj<N; sj++)
				denominator_sum += alpha[t][si]*A[si][sj]*B[sj][Oi[t+1]-1]*beta[t+1][sj];
		}
		for(int si=0; si<N; si++){
			for(int sj=0; sj<N; sj++)
				XI[t][si][sj]=(alpha[t][si]*A[si][sj]*B[sj][Oi[t+1]-1]*beta[t+1][sj])/denominator_sum;
		}
	}
}//P3_Baum_Welch_Procedure

/**************************************************************************************************
	P3 Re-estimation Procedure | Calculating Abar Bbar Pibar 
	XI: Prob of being in state Si at time t and state Sj at time t+1 given the model lambda
		and observation sequence.
	Input: Observation Sequnce, Observation Count
**************************************************************************************************/
void P3_Reestimation_Procedure(int *Oi, int o){

	int T = OFmax[o];

	// Re-estimation of Pi as Pi_bar
	for(int si=0; si<N; si++) {
		PI_BAR[si]=gamma[0][si];
	}//Pi as Pi_bar

	// Re-estimation of A as A_bar
	long double nume_exp_num_of_transitionf_from_si_to_sj=0;
	long double deno_exp_num_of_transitionf_from_si=0;

	for(int si=0; si<N; si++) 
	{
		for(int sj=0; sj<N; sj++)
		{
			nume_exp_num_of_transitionf_from_si_to_sj=0;
			deno_exp_num_of_transitionf_from_si=0;
			for(int t=0; t<T-2; t++)  // from (1 to T-1)
			{
				nume_exp_num_of_transitionf_from_si_to_sj += XI[t][si][sj];
				deno_exp_num_of_transitionf_from_si += gamma[t][si];
			}

			A_BAR[si][sj]=nume_exp_num_of_transitionf_from_si_to_sj/deno_exp_num_of_transitionf_from_si;
		}
	}//A as A_bar

	// Re-estimation of B as B_bar
	long double nume_exp_num_of_times_in_sj_observing_vk =0;
	long double deno_exp_num_of_times_in_sj =0;

	for(int sj=0; sj<N; sj++)
	{
		for(int k=0; k < M; k++)  // obs seq value are from 1 to M
		{
			nume_exp_num_of_times_in_sj_observing_vk=0;
			deno_exp_num_of_times_in_sj=0;
			for(int t=0; t<T; t++)
			{
				if(Oi[t]-1==k)
					nume_exp_num_of_times_in_sj_observing_vk += gamma[t][sj];

				deno_exp_num_of_times_in_sj+=gamma[t][sj];
			}
			B_BAR[sj][k]=nume_exp_num_of_times_in_sj_observing_vk/deno_exp_num_of_times_in_sj;
		}
	}// B as B_bar


}//P3_Reestimation_Procedure

/**************************************************************************************************
	Print Model Lamda BAR in the File. (or Display on console)
**************************************************************************************************/
void print_model_lambda_bar(){

	fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Lambda(A):\n");
		for (int si = 0; si < N; ++si){
			fprintf(fp_console, "S[%d]\t", si+1);
			for (int sj = 0; sj < N; ++sj){
					//printf("%0.16g\t", A[si][sj]);
					fprintf(fp_console, "%.16g (%d)\t",A_BAR[si][sj], sj+1);
				}
			fprintf(fp_console, "\n");
		}

	fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Lambda(B):\n");
		for (int si = 0; si < N; ++si){
			fprintf(fp_console, "S[%d]\t", si+1);
			for (int m = 0; m < M; ++m){
					//printf("%.16g\t", B[si][m]);
					fprintf(fp_console, "%.16g (%d)\t",B_BAR[si][m], m+1);
				}
			fprintf(fp_console, "\n");
		}
		
				
	fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Lambda(PI):\n");
			for (int si = 0; si < N; ++si){
				fprintf(fp_console, "%.16g (S%d)\t",PI_BAR[si], si+1);
			}
			fprintf(fp_console, "\n");

}

/**************************************************************************************************
	Making Lambda Stochastic by adjusting values in max value in row.
**************************************************************************************************/
void make_lambda_stochastic(){
	long double row_sum = 0;
	long double row_max =0;
	int max_index = 0;

	// making A Stochastic
	for(int si=0; si<N; si++){
		row_sum=0;
		row_max = 0;
		max_index =-1;		

		for(int sj=0; sj<N; sj++){
			if( A_BAR[si][sj] > row_max )
			{	max_index = sj;
				row_max = A_BAR[si][sj];
			}

			row_sum += A_BAR[si][sj];
		}//row sum
		//if(row_sum != 1){
		//	A_BAR[si][max_index] = (row_sum > 1) ? (A_BAR[si][max_index] - (row_sum-1) ):(A_BAR[si][max_index] + (1-row_sum));
		//}
		A_BAR[si][max_index] -= (row_sum-1); 
		
	}// making A Stochastic

	// making B Stochastic
	for(int sj=0; sj<N; sj++){
		row_sum=0;
		row_max = 0;
		max_index =0;	

		for(int k=0; k<M; k++){
		
			if(B_BAR[sj][k] == 0)
				B_BAR[sj][k] = 1e-030;

			if( B_BAR[sj][k] > row_max )
			{	max_index = k;	
				row_max =  B_BAR[sj][k];
			}

			row_sum += B_BAR[sj][k];

		}
		//if(row_sum != 1){
		//	B_BAR[sj][max_index] = (row_sum > 1) ? (B_BAR[sj][max_index] - (row_sum-1) ):(B_BAR[sj][max_index] + (1-row_sum));
		//}
		B_BAR[sj][max_index] -= (row_sum-1);
	}// making B Stochastic

}//make_lambda_stochastic

/**************************************************************************************************
	Replace Old Model Lambda with New Model LambdaBar after Re-estimation
**************************************************************************************************/
void replace_old_model(){

	for(int si=0; si<N; si++) 
        PI[si]=PI_BAR[si];

    for(int si=0; si < N; si++) //assign the given values for transition probability distribution
        for(int sj=0; sj<N; sj++)
            A[si][sj]=A_BAR[si][sj];

    for(int sj=0; sj<N; sj++) //assign the given values for observation symbol probability distribution
        for(int k=0; k<M; k++)
            B[sj][k]=B_BAR[sj][k];


}//replace_old_model

/**************************************************************************************************
	Print Model Lamda in the File. (or Display on console)
**************************************************************************************************/
void print_model_lambda(){

	fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Lambda(A):\n");
		for (int si = 0; si < N; ++si){
			fprintf(fp_console, "S[%d]\t", si+1);
			for (int sj = 0; sj < N; ++sj){
					//printf("%0.16g\t", A[si][sj]);
					fprintf(fp_console, "%.16g (%d)\t",A[si][sj], sj+1);
				}
			fprintf(fp_console, "\n");
		}

	fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Lambda(B):\n");
		for (int si = 0; si < N; ++si){
			fprintf(fp_console, "S[%d]\t", si+1);
			for (int m = 0; m < M; ++m){
					//printf("%.16g\t", B[si][m]);
					fprintf(fp_console, "%.16g (%d)\t",B[si][m], m+1);
				}
			fprintf(fp_console, "\n");
		}
		
				
	fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Lambda(PI):\n");
			for (int si = 0; si < N; ++si){
				fprintf(fp_console, "%.16g (S%d)\t",PI[si], si+1);
			}
			fprintf(fp_console, "\n");

}

/**************************************************************************************************
	Feed Forward Model | Bakes Model
**************************************************************************************************/
void initialize_feed_forward_model()
{
	/**************** PI ****************/
	PI[0]=1.0; // first state
	for(int si=1; si<N; si++){
		PI[si]=0;	// rest of the state
	}

	/**************** A ****************/
	for(int si=0; si<N; si++)
	{
		for(int sj=0; sj<N; sj++)
		{
			if(si==sj)
				A[si][sj] = 0.8; // Prob of being in same state
			else if(si+1==sj)
				A[si][sj] = 0.2; // Prob to shift in next state
			else A[si][sj] = 0;
		}
	}
	A[N-1][N-1] = 1;		// Forcing to be in same state reaching final state

	/**************** B ****************/
	for(int sj=0; sj<N; sj++)
		for(int k=0; k<M; k++)
			B[sj][k] = 1.0/M;


}//initialize_feed_forward_model

/**************************************************************************************************
	Using Previous Converged Model as starting point
**************************************************************************************************/
void initialize_converged_model()
{
	/**************** PI ****************/
	PI[0]=1.0; // first state
	for(int si=1; si<N; si++){
		PI[si]=0;	// rest of the state
	}

	/**************** A ****************/
	for(int si=0; si<N; si++)
	{
		for(int sj=0; sj<N; sj++)
		{
			A[si][sj] = A_Prev[si][sj];
		}
	}

	/**************** B ****************/
	for(int sj=0; sj<N; sj++)
		for(int k=0; k<M; k++)
			B[sj][k] = B_Prev[sj][k];


}//initialize_feed_forward_model

/**************************************************************************************************
	Saving Converged Model to their Individual files
	Input: digit value
**************************************************************************************************/
void output_lambdaABPi_to_each_file(int d){
	FILE *fp_ind;
	char completePathInd[200];

// Save PI
	sprintf(completePathInd, "%s%s/%s/%s%s.txt", output_folder, output_folder_Model_name, WordNames[d], LambdaFileNames[2], WordNames[d]);  
	fp_ind = fopen(completePathInd, "w");				//to read input codebook
	if(fp_ind == NULL){ 
		perror("\n Error: ");
		printf("\n File Name is: %s", completePathInd);
		getch();
		return;
	}

	for (int si = 0; si < N; ++si){
		fprintf(fp_ind, "%.16g\t",PI[si]);
	}

	fflush(fp_ind); fclose(fp_ind); 

//Save A
	sprintf(completePathInd, "%s%s/%s/%s%s.txt", output_folder, output_folder_Model_name, WordNames[d], LambdaFileNames[0], WordNames[d]);  
	fp_ind = fopen(completePathInd, "w");				//to read input codebook
	if(fp_ind == NULL){ 
		perror("\n Error: ");
		printf("\n File Name is: %s", completePathInd);
		getch();
		return;
	}
	for (int si = 0; si < N; ++si){
		for (int sj = 0; sj < N; ++sj){
			fprintf(fp_ind, "%.16g\t",A[si][sj]);
		}
		fprintf(fp_ind, "\n");
	}
	fflush(fp_ind); fclose(fp_ind); 


//Save B
	sprintf(completePathInd, "%s%s/%s/%s%s.txt", output_folder, output_folder_Model_name, WordNames[d], LambdaFileNames[1], WordNames[d]);  
	fp_ind = fopen(completePathInd, "w");				//to read input codebook
	if(fp_ind == NULL){ 
		perror("\n Error: ");
		printf("\n File Name is: %s", completePathInd);
		getch();
		return;
	}
	for (int si = 0; si < N; ++si){
		for (int m = 0; m < M; ++m){
		
			fprintf(fp_ind, "%.16g\t",B[si][m]);
		}
		fprintf(fp_ind, "\n");
	}
	fflush(fp_ind); fclose(fp_ind); 

}

/**************************************************************************************************
	Convergence Procedure: use each utterance and then converging until Pstar improved
	then average model of all the converged models
**************************************************************************************************/
void covergence_procedure(){

	int skip=1;				// number of times ignore (Pstar_new > Pstar_old) and check once again for more Pstar_new

	for(int rc=1; rc<=repeatConvergence;rc++)//repeat using converged model second time
		{ 
			for(int o=0; o<On; ++o){
				
				if(rc==1){
					printf("\n  -###-###-###-###-###-###->>> Bakis Model | Observation O[%d] < ----------------------- ----------------------- \n", o+1);	
					fprintf(fp_console, "\n  -###-###-###-###-###-###->>> Bakis Model | Observation O[%d] < ----------------------- ----------------------- \n", o+1);	
					initialize_feed_forward_model();
					skip=1;
				}
				else {
					printf("\n -###-###-###-###-###-###->>> Converged Model | Observation O[%d] < ----------------------- ----------------------- \n", o+1);	
					fprintf(fp_console, "\n -###-###-###-###-###-###->>> Converged Model | Observation O[%d] < ----------------------- ----------------------- \n", o+1);	
					initialize_converged_model();
					skip=0;
				}

				printf("\n");
			
				int itr=1;
				do
				{
					P1_Forward_Procedure(O[o],o);
					P1_Backward_Procedure(O[o],o);
					P2_Gamma_Procedure(o);
					P2_Viterbi_Algo(O[o],o);

					if(showAlphaBetaPstarInConsole){
						printf("\n[%d]\t", itr); 
						printf("Alpha P = %g \t", probability_alpha);
						printf("Beta P = %g \t", probability_beta);
						printf("Pstar P = %g \t",Pstar);
					}
					P3_Baum_Welch_Procedure(O[o],o);
					P3_Reestimation_Procedure(O[o],o);
					make_lambda_stochastic();
					replace_old_model();

					itr++;
					Pstar_old = Pstar;
					P2_Viterbi_Algo(O[o],o);
					if(showAlphaBetaPstarInConsole){
						printf("New_Pstar P = %g\t",Pstar);
					}
				}while( (Pstar > Pstar_old && itr <= model_iterations) || (skip--));// while Pstar> Pstar_old loop

				//if(rc==3)
				//for(; itr<=model_iterations; itr++)
				//{
				//	P1_Forward_Procedure(O[o],o);
				//	P1_Backward_Procedure(O[o],o);
				//	P2_Gamma_Procedure(o);
				//	P2_Viterbi_Algo(O[o],o);

				//	if(showAlphaBetaPstarInConsole){
				//		printf("\n[%d]\t", itr); 
				//		printf("Alpha P = %g \t", probability_alpha);
				//		printf("Beta P = %g \t", probability_beta);
				//		printf("Pstar P = %g \t",Pstar);
				//	}
				//	P3_Baum_Welch_Procedure(O[o],o);
				//	P3_Reestimation_Procedure(O[o],o);
				//	make_lambda_stochastic();
				//	replace_old_model();

				//}//for each iterations itr

				// save converged lambda
				for (int si = 0; si < N; ++si){
					for (int sj = 0; sj < N; ++sj){
						converged_A[o][si][sj] = A[si][sj];
					}
				}

				for (int si = 0; si < N; ++si){
					for (int m = 0; m < M; ++m){
						converged_B[o][si][m] = B[si][m];
					}
				}

				printf(" Total Iterations: %d\t", itr-1); 
				printf(" Alpha P = %g \t", probability_alpha);
				printf( "Beta P = %g \t", probability_beta);
				printf( "Pstar P = %g \n",Pstar);
				
				fprintf(fp_console, "Total Iterations: %d\t", itr-1); 
				fprintf(fp_console, "Alpha P = %g \t", probability_alpha);
				fprintf(fp_console, "Beta P = %g \t", probability_beta);
				fprintf(fp_console, "Pstar P = %g \n",Pstar);
			}//for each observation seq 'o'

			// take averaged lambda
			long double lambda_sum=0;
			for (int si = 0; si < N; ++si){
				for (int sj = 0; sj < N; ++sj){
					lambda_sum=0;

					for(int u=0; u<On; u++)
					{
						lambda_sum += converged_A[u][si][sj];
					}

					A[si][sj] = lambda_sum/On ;
					A_Prev[si][sj] = A[si][sj] ;
				}
			}

			for (int si = 0; si < N; ++si){
				for (int m = 0; m < M; ++m){
					lambda_sum=0;
					for(int u=0; u<On; u++)
					{
						lambda_sum += converged_B[u][si][m];
					}
					B[si][m] = lambda_sum/On ;
					B_Prev[si][m] = B[si][m];
				}
			}
		
			fprintf(fp_console, "\n ----------------------- After Taking Average of Converged Lambdas of %d Training Sequence  \n", On);	
			
			if(showStateSeqAlphaBetaInFileForEachObsAfterConverge){
				for(int o=0; o<On; ++o){
					P2_Viterbi_Algo(O[o],o);
					PRINT_P2_PStar_State_Sequence(O[o],o);

						fprintf(fp_console, "\nO[%d]:\t",o+1);
						fprintf(fp_console, "Alpha P = %g\t",P1_Forward_Procedure(O[o],o));
						fprintf(fp_console, "Beta P = %g\t",P1_Backward_Procedure(O[o],o));
						fprintf(fp_console, "Pstar P = %g\t\n",Pstar);
				}
			}
			fprintf(fp_console, "\n ----------------------- -----------------------> Converged Lamda <----------------------- ----------------------- \n");	
			
			print_model_lambda();

		}// for each rc

}//covergence_procedure

/**************************************************************************************************
	Replace Model in Input Folder with Output Folder.
**************************************************************************************************/
void replace_old_models_files(){

	char from_lambda_location[300], to_lambda_location[300];
	FILE *fp_infile, *fp_outfile;
	long double temp;

	for(int d=0; d<W; d++){

		printf("\n ---#---#---#---#---#---#---#--->>> Replacing Model: %s %s <---#---#---#---#---#---#---#---\n", WordFolderName, WordNames[d]); 

// replace A
		sprintf(from_lambda_location, "%s%s/%s/%s%s.txt", output_folder, output_folder_Model_name, WordNames[d], LambdaFileNames[0], WordNames[d]);  // {output_folder/Models}+{0/}+{A_}+{1}+".txt 
		sprintf(to_lambda_location, "%s%s %s/%s%s.txt", input_folder, WordFolderName, WordNames[d], LambdaFileNames[0], WordNames[d]);  // {output_folder/Models}+{0/}+{A_}+{1}+".txt 
	
		fp_infile = fopen(from_lambda_location, "r");
		fp_outfile = fopen(to_lambda_location, "w");
			if(fp_infile == NULL || fp_outfile == NULL ){ 
				perror("\n Error: ");
				printf("\n File Names are: \nSrc: %s \nDest: %s", from_lambda_location, to_lambda_location);
				getch();
				return;
			}

		for (int si = 0; si < N; ++si){
			for (int sj = 0; sj < N; ++sj){
				fscanf(fp_infile, "%Lf",&temp);
				fprintf(fp_outfile, "%.16g\t",temp);
			}
			fprintf(fp_outfile, "\n");
		}

		fflush(fp_infile); fclose(fp_infile); 
		fflush(fp_outfile); fclose(fp_outfile); 

// replace B
		sprintf(from_lambda_location, "%s%s/%s/%s%s.txt", output_folder, output_folder_Model_name, WordNames[d], LambdaFileNames[1], WordNames[d]);  // {output_folder/Models}+{0/}+{B_}+{1}+".txt 
		sprintf(to_lambda_location, "%s%s %s/%s%s.txt", input_folder, WordFolderName, WordNames[d], LambdaFileNames[1], WordNames[d]);  // {output_folder/Models}+{0/}+{B_}+{1}+".txt 
		
		fp_infile = fopen(from_lambda_location, "r");
		fp_outfile = fopen(to_lambda_location, "w");
			if(fp_infile == NULL || fp_outfile == NULL ){ 
				perror("\n Error: ");
				printf("\n File Names are: \nSrc: %s \nDest: %s", from_lambda_location, to_lambda_location);
				getch();
				return;
			}

		for (int si = 0; si < N; ++si){
			for (int m = 0; m < M; ++m){
				fscanf(fp_infile, "%Lf",&temp);
				fprintf(fp_outfile, "%.16g\t",temp);
			}
			fprintf(fp_outfile, "\n");
		}
		fflush(fp_infile); fclose(fp_infile); 
		fflush(fp_outfile); fclose(fp_outfile); 

//replace pi
		sprintf(from_lambda_location, "%s%s/%s/%s%s.txt", output_folder, output_folder_Model_name, WordNames[d], LambdaFileNames[2], WordNames[d]);  // {output_folder/Models}+{0/}+{Pi_}+{1}+".txt 
		sprintf(to_lambda_location, "%s%s %s/%s%s.txt", input_folder, WordFolderName, WordNames[d], LambdaFileNames[2], WordNames[d]);  // {output_folder/Models}+{0/}+{Pi_}+{1}+".txt 
		fp_infile = fopen(from_lambda_location, "r");
		fp_outfile = fopen(to_lambda_location, "w");
			if(fp_infile == NULL || fp_outfile == NULL ){ 
				perror("\n Error: ");
				printf("\n File Names are: \nSrc: %s \nDest: %s", from_lambda_location, to_lambda_location);
				getch();
				return;
			}
		for (int si = 0; si < N; ++si){
			fscanf(fp_infile, "%Lf",&temp);
			fprintf(fp_outfile, "%.16g\t",temp);
		}

		fflush(fp_infile); fclose(fp_infile); 
		fflush(fp_outfile); fclose(fp_outfile);
	}// for each digit d<W

}

/**************************************************************************************************
	Test Offline Utterance of the Digits
**************************************************************************************************/
void offline_testing(int test_word, unsigned int model_type_to_use){

	long double cur_alpha_probability=0;
	long double max_probability=0;
	int word_index=-1;

	offline_correct_count=0;

	for(int u=0; u<Ot; u++){
					printf("\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Utterance: %s_%s/%s_O[%d] \n", WordFolderName, WordNames[test_word], LambdaFileNames[4], u+1);
						fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Utterance: %s_%s/%s_O[%d] \n", WordFolderName, WordNames[test_word], LambdaFileNames[4], u+1); 	 
			max_probability=0;
			cur_alpha_probability=0;
			word_index=-1;
		for(int w=0; w<W; w++){
			readLambdaABPi(w, model_type_to_use);
			cur_alpha_probability=P1_Forward_Procedure(O[u],u);

					printf("O[%d]:W[%s]::  Alpha P = %g\n", u+1, WordNames[w], cur_alpha_probability); 	 
					fprintf(fp_console, "O[%d]:W[%s]::  Alpha P = %g\n", u+1, WordNames[w], cur_alpha_probability); 	 

			if(cur_alpha_probability > max_probability){
				max_probability = cur_alpha_probability;
				word_index=w;
			}
		}//for each w word

		if(word_index==-1)
		{
			printf("\n -----------------> Actual digit: %s  ", WordNames[test_word]);
				printf("\n -----------------> Digit Recognized: %s\n", "NOT RECOGNIZED");
			fprintf(fp_console, "\n -----------------> Actual digit: %s  ", WordNames[test_word]);
			fprintf(fp_console, "\n -----------------> Digit Recognized: %s\n", "NOT RECOGNIZED");
		}
		else
		{
			printf("\n -----------------> Actual digit: %s  ", WordNames[test_word]);
				printf("\n -----------------> Digit Recognized: %s\n", WordNames[word_index]);
			fprintf(fp_console, "\n -----------------> Actual digit: %s  ", WordNames[test_word]);
			fprintf(fp_console, "\n -----------------> Digit Recognized: %s\n", WordNames[word_index]);
		}	
		
		if(strcmp(WordNames[test_word], WordNames[word_index])==0)offline_correct_count++;
	}// for each u utterance

	offline_overall_count +=offline_correct_count;

	fprintf(fp_console, "\n ------------------------------------------------------------------------\n"); 
}//offline_testing

/**************************************************************************************************
	Test Live Utterance of the Digits
**************************************************************************************************/
void live_testing(unsigned int model_type_to_use){

	
	/**************** Observation Sequence Generation ****************/
	live_sequence_generation();
	/**************** Read Observation Sequence ****************/
	readObsSeqOfLiveRecording();

	/**************** Testing ****************/
	int NumOfLiveUtterance=1;

	long double cur_alpha_probability=0;
	long double max_probability=0;
	int word_index=-1;

	offline_correct_count=0;

	for(int u=0; u<NumOfLiveUtterance; u++){
					printf("\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Live Utterance: %s_obs_seq_O[%d] \n", liveRecordingFileName, u+1);
						fprintf(fp_console, "\n ~~~~~~~~~~~~~~~~~~~~~~~~~> Live Utterance: %s_obs_seq_O[%d] \n", liveRecordingFileName, u+1); 	 
			max_probability=0;
			cur_alpha_probability=0;
			word_index=-1;
		for(int w=0; w<W; w++){
			readLambdaABPi(w, model_type_to_use);
			cur_alpha_probability=P1_Forward_Procedure(O[u],u);

					printf("O[%d]:W[%s]::  Alpha P = %g\n", u+1, WordNames[w], cur_alpha_probability); 	 
					fprintf(fp_console, "O[%d]:W[%s]::  Alpha P = %g\n", u+1, WordNames[w], cur_alpha_probability); 	 

			if(cur_alpha_probability > max_probability){
				max_probability = cur_alpha_probability;
				word_index=w;
			}
		}//for each w word
		if(word_index==-1)
		{
			printf("\n -----------------> Digit Recognized: %s\n", "NOT RECOGNIZED");
			fprintf(fp_console, "\n -----------------> Digit Recognized: %s\n", "NOT RECOGNIZED" );
		}
		else
		{
			printf("\n -----------------> Digit Recognized: %s\n", WordNames[word_index]);
			fprintf(fp_console, "\n -----------------> Digit Recognized: %s\n", WordNames[word_index]);
		}
		
	}// for each u utterance
	fprintf(fp_console, "\n ------------------------------------------------------------------------\n"); 
}//live_testing