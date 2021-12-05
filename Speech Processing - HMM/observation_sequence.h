
/******************--------Common Settings--------******************/
//Speech Samples
#define scaleAmp 5000							  // Max Amplitutde Value to Scale.
#define initIgnoreHeaderLines 5			//change  // Number of Initial Header Lines in Txt Files to Ignore.
#define initIgnoreSamples 6400				//change(30frames)  // Number of Initial samples to ignore.
#define initNoiseFrames 10				//change  // Number of Initial Noise Frames to Consider.
#define thresholdNoiseToEnergyFactor 3	//change  // Noise to Energy Factor.
#define samplingRate 16000						  // Sampling Rate of Recording.
#define sizeFrame 320		//change // Number of Samples per Frame.
#define initIgnoreEndFrames 30

const unsigned short totDigits = W;
const unsigned short MaxFrames = Tmax+10;				//change	// max number of frames to consider overall.
const unsigned short q=p;				//change	//#Coefficients (c_i's) that need to be found
const unsigned short NumOfTrainingFiles = On;
const unsigned short NumOfTestingFiles = Ot;

/******************--------Variable Decalaration--------******************/
// DISTANCE, WEIGHTS
const double w_tkh[]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};		//tokhura weigths provided.

// COEFFICIENTS
double w_rsw[q+1];						// weights for raised sine window.
double C_rsw[Tmax][q+1]; // calculating Raised Sine Window Cepstral coefficient and storing for each frame
int OBS_SEQ[Tmax];

// Speech
int samples[500000], maxAmp = 0;							//max number of samples in recording to consider, maximum amplitute for rescaling.
double normSamples[500000];									//normalised samples.
long cntTotSamples = 0, cntTotFrames = 0;					// Total Samples in recording, Total Frames in Recording.
float DCshift =0;
long start=10, end=Tmax-10 ;										//start and end marker, Frames

// ZCR and ENERGY
float cntZCR[MaxFrames], avgZCR[MaxFrames], avgEnergy[MaxFrames];  
float totEnergy=0, noiseEnergy=0, initNoiseZCR=0, initNoiseAvgZCR=0;;
float thresholdZCR=0;
double thresholdEnergy=0;
double maxEnergy = 0;										// Max Energy of the Frame.
int STE_Marker;												// Max STE Marker For the Frame

// File Stream
FILE *fp_obs_seq_ip, *fp_obs_seq_norm, *fp_obsseq_norm_seg, *fp_obsseq_final_op, *fp_obsseq_console;							// file pointer for input stream and output stream.
char OSeqfileNameIp[300], OSeqcompletePathIp[300], OSeqcompletePathNorm[300], OSeqcompletePathNormSeg[300], OSeqcompletePathFinOp[300], OSeqcompletePathConsole[300];
char OSeqcharsLineToRead[50];															// number of characters to read per line in file.
const char filePathInputVoiceTraining[] = "input_voice_training_data/";		 							// Folder name where all Vowels recordings are placed 1 to 100.	input_digits_data/
const char filePathInputVoiceTesting[] = "input_voice_testing_data/";	
const char fileOutputRecordingNorSeg[] = "output_voice_recordings_normalised_segregated/";							// Folder name where output of segregated digits recordings are placed and its normalised files. output_recordings/
const char fileOutputRecordingAnalysis[] = "output_voice_recordings_analysis_files/";								// Folder name where test file analysis are saved. output_test_result/


/**************************************************************************************
	Calculating w[m] weights for raised sine window and storing in array for future use.
**************************************************************************************/
void CalculateWeightsForRaisedSineWindow(){
	for(int i=1;i<=q;i++){
		w_rsw[i] = 1 + (q/2)*sin(M_PI*i/q);
	}
}

/**************************************************************************************
	Calculating R[i] Values using Auto Correlation Method
	Input:	   *s is the pointer to our samples array.
	Output:	   *Rn is the pointer to our R array for storing final Ri values (0 to p).
**************************************************************************************/
void calculateRi_AutoCorrelation(double *s, double *Rn){
	for(int k=0;k<=p;k++) 
	{ 
		Rn[k]=0;
		for(int m=0; m < sizeFrame-k; m++)			
		{
			Rn[k] += s[m]*s[m+k];				// R[0] R[1] .... R[p]
		}	
	}
}

/**************************************************************************************
	Calculating A[i] Values using Levinson Durbin Process
	Input:	   *R is the pointer to our Auto Correlated Energy part.
	Output:	   *A is the pointer to our A array for storing final Ai values (1 to p).
**************************************************************************************/
void calculateAi_Durbin(double *R, double *A){
	 double E[p+1] = {0};
	 double k[p+1] = {0};
	 double alpha[p+1][p+1] = {0};
	 double sum=0;
	 int i, j;

	E[0]=R[0]; //Initialize E[0] with energy

	for( i=1;i<=p;i++)
	{
		if(i==1)
			k[1]=R[1]/R[0]; //special case i=1
		else //find k(i) for all other remaining values of i
		{
			 sum=0;
			for( j=1;j<=i-1;j++)
			{
				sum+=alpha[i-1][j]*R[i-j];
			}
			k[i]=(R[i]-sum)/E[i-1];
		}
		alpha[i][i]=k[i];

		for( j=1;j<=i-1;j++)
			alpha[i][j]=alpha[i-1][j]-k[i]*alpha[i-1][i-j]; //update coefficients from previous values
		
		E[i]=(1-k[i]*k[i])*E[i-1]; //update E[i]
	}
	for( i=1;i<=p;i++){	
		//A[i] = -1*alpha[p][i];
		A[i] = alpha[p][i];				// A[0] A[1] .... A[p]
	}
}

/**************************************************************************************
	Calculating C[i] Cepstral Coefficient Values.
	Input:	   *A is the pointer to our Ai Array.
				sigma is R[0] value.
	Output:	   *A is the pointer to our A array for storing final Ai values (1 to p).
**************************************************************************************/
void calculateCi_Cepstral(double sigma, double *A, double *c){
	int k,m;
	double sum=0;

	c[0]= logl(sigma*sigma); 
	for(m=1; m<=p; m++)
	{
		sum=0;
		for(k=1;k<=m-1;k++) //sum from older cepstral coefficents to compute new cepstral coefficients
			sum+=k*c[k]*A[m-k]/m;
		c[m]=A[m]+sum;		//new cepstral coefficients
	}
	/*
		if(m>p) // For This Assignment this never get executed as we assume q=p
		{
			for(;m<=q;m++)
			{
				sum=0;
				for(k=1;k<=m-1;k++) //sum from older cepstral coefficents to compute new cepstral coefficients
					sum+=k*c[k]*A[m-k]/m;
				c[m]=sum;		//new cepstral coefficients
			}
		}
	*/
}


/**************************************************************************************
	To Display Common Settings used in our System
	Input: File Pointer in case things needed to be written on file.
**************************************************************************************/
void DisplayCommonSettingsObsSeq(FILE *fp_op=NULL){
	// General Information to Display
	if(fp_op==NULL){
		printf("****-------- --------****\n");		
		printf("-Common Settings are : -\n");	
		printf(" P (=Q) = : %d\n", p);
		printf(" Frame Size : %d\n", sizeFrame);	
		printf(" Num of Files For Training : %d\n", NumOfTrainingFiles);	
		printf(" Num of Files For Testing : %d\n", NumOfTrainingFiles);	

		printf(" Tokhura Weights : ");
		for(int i=0; i<q; i++){
			printf("%0.1f(%d) ", w_tkh[i],i+1);
		}
		printf("\n");
		printf(" Amplitutde Value to Scale : %d\n", scaleAmp);			
		printf(" Intital Header Lines Ignore Count : %d\n", initIgnoreHeaderLines); 
		printf(" Intital Samples to Ignore : %d\n",initIgnoreSamples);	
		printf(" Intital Noise Frames Count : %d\n",initNoiseFrames);	
		printf(" Noise to Energy Factor : %d\n",thresholdNoiseToEnergyFactor); 
		printf(" Sampling Rate of Recording: %d\n",samplingRate); 
		printf("----------------------------------------------------------------\n\n");		
	}
	else if(fp_op!=NULL){
		fprintf(fp_op, "****----------------****\n");		
		fprintf(fp_op, "-Common Settings are : -\n");	
		fprintf(fp_op, " P (=Q) = : %d\n", p);
		fprintf(fp_op, " Frame Size : %d\n", sizeFrame);	
		fprintf(fp_op, " Num of Files For Training : %d\n", NumOfTrainingFiles);	
		fprintf(fp_op, " Num of Files For Testing : %d\n", NumOfTrainingFiles);	

		fprintf(fp_op, " Tokhura Weights : ");
		for(int i=0; i<q; i++){
			fprintf(fp_op, "%0.1f(%d) ", w_tkh[i],i+1);
		}
		fprintf(fp_op, "\n");
		fprintf(fp_op, " Amplitutde Value to Scale : %d\n", scaleAmp);			
		fprintf(fp_op, " Intital Header Lines Ignore Count : %d\n", initIgnoreHeaderLines); 
		fprintf(fp_op, " Intital Samples to Ignore : %d\n",initIgnoreSamples);	
		fprintf(fp_op, " Intital Noise Frames Count : %d\n",initNoiseFrames);	
		fprintf(fp_op, " Noise to Energy Factor : %d\n",thresholdNoiseToEnergyFactor); 
		fprintf(fp_op, " Sampling Rate of Recording: %d\n",samplingRate); 
		fprintf(fp_op, "----------------------------------------------------------------\n\n");	
	}
}

/**************************************************************************************
	Normalising and DC Shift of Samples.
	Input (global):		fp_obs_seq_ip is to read from recording text file.
						fp_obs_seq_norm is to save normalised samples into another file.
						fp_obsseq_console is to save analysed values in common output file.
**************************************************************************************/
void normalize_dcshift_samples(){
			cntTotSamples=0; maxAmp = 0;

			int totIgnore;
				 totIgnore=initIgnoreHeaderLines + initIgnoreSamples; //5 + 2 = 7
			// totIgnore=0;
			long sample_index = totIgnore+1; // till 7 samples ignored, sample count is 8, so to make array index 0 there is +1 
			long sample_index_norm = 0;
			double normFactor = 0;
			double normOutput = 0;

				while(!feof(fp_obs_seq_ip)){
					fgets(OSeqcharsLineToRead, 50, fp_obs_seq_ip);
					cntTotSamples += 1 ;  

					if(cntTotSamples > totIgnore){
						sample_index_norm = cntTotSamples - sample_index;
						samples[sample_index_norm] = (int)atoi(OSeqcharsLineToRead);
						DCshift += samples[sample_index_norm];
								if(abs(samples[sample_index_norm]) > maxAmp)
									maxAmp = abs(samples[sample_index_norm]);
					}
				}

			cntTotSamples = cntTotSamples - totIgnore;		// total number of samples stored in array
			DCshift = DCshift/cntTotSamples;				//average DC Shift
			cntTotFrames = (int)floor((1.0*cntTotSamples)/sizeFrame);			//total number of frames.

			start=10;
			end = cntTotFrames-10;

			normFactor = (double)scaleAmp/maxAmp;			//normalising factor

			// saving in normalised file
			for(long i=0; i<cntTotSamples; i++){
				normOutput = (double)(samples[i] - DCshift)*normFactor;
				normSamples[i]=normOutput;
				fprintf(fp_obs_seq_norm, "%lf\n", normSamples[i]);
			}
			
			
	//printf(" TOTAL SAMPLES : %d\n TOTAL FRAMES : %d\n DC SHIFT needed : %lf\n Maximum Amplitude : %d\n Normalization Factor : %lf\n ", cntTotSamples, cntTotFrames, DCshift, maxAmp, normFactor);
	if(fp_obsseq_console!=NULL){
		fprintf(fp_obsseq_console, " TOTAL SAMPLES : %d\n TOTAL FRAMES : %d\n DC SHIFT needed : %lf\n Maximum Amplitude : %d\n Normalization Factor : %lf\n ", cntTotSamples, cntTotFrames, DCshift, maxAmp, normFactor);
	}
}

/**************************************************************************************
	Calculating ZCR and Energy of Frames
	Input (global): fp_obs_seq_norm is to read from normalised file of the samples.
**************************************************************************************/
void zcr_energy_frames(){
	rewind(fp_obs_seq_norm);

	long i,j;
	float s_i, s_i_1=1;

	//totEnergy=0;
	 maxEnergy = 0;										// Max Energy of the Frame.
	 STE_Marker = 0;	
	
	for(i=0;i < cntTotFrames;i++)
		{
			cntZCR[i]=0;
			avgEnergy[i]=0;
			for(j=0;j < sizeFrame ;j++)
			{
				fgets(OSeqcharsLineToRead, 50, fp_obs_seq_norm); // reading from normalised input
				s_i = (float)atof(OSeqcharsLineToRead);
				avgEnergy[i] += (s_i*s_i);
				cntZCR[i] +=  (s_i_1*s_i < 0);
				s_i_1 = s_i;
			}	
			avgEnergy[i]/=sizeFrame;
			avgZCR[i] = cntZCR[i]/sizeFrame;
			//totEnergy+=avgEnergy[i];
			//fprintf(fp_obs_seq_norm, "%f %0.1f \n", avgEnergy[i], cntZCR[i]);	//dumping the features of frames.
			// calculation for detecting STE Frame
			if(avgEnergy[i] > maxEnergy){
				maxEnergy = avgEnergy[i];
				STE_Marker=i;
			}
		}

}


/**************************************************************************************
	Calculating ZCR and Energy of Noise Frames, and Finally Thresholds for ZCR and Energy
	Input (global): fp_obsseq_console is to save analysed values in common output file.
**************************************************************************************/
void noiseEnergy_thresholds_frames(){
	noiseEnergy=0; initNoiseZCR=0; initNoiseAvgZCR=0;
	int i;
	for(i=0;i < initNoiseFrames;i++){
			initNoiseZCR+=cntZCR[i];
			initNoiseAvgZCR+=avgZCR[i];
			noiseEnergy+=avgEnergy[i];
	}
		thresholdZCR=initNoiseZCR/initNoiseFrames;
		noiseEnergy/=initNoiseFrames;
		thresholdEnergy=noiseEnergy*thresholdNoiseToEnergyFactor;

	//printf( "\n---- Initial Noise Frames : %d ----\n\n", initNoiseFrames);
	//printf(" Avg Noise Energy : %lf\n Total Noise ZCR : %0.1f\n Threshold ZCR : %0.1f\n Threshold Energy(Avg Noise*%d) : %0.5lf\n ", noiseEnergy, initNoiseZCR, thresholdZCR, thresholdNoiseToEnergyFactor, thresholdEnergy);
	if(fp_obsseq_console!=NULL){
		fprintf(fp_obsseq_console, "\n---- Initial Noise Frames : %d ----\n\n", initNoiseFrames);
		fprintf(fp_obsseq_console, " Avg Noise Energy : %lf\n Total Noise ZCR : %0.1f\n Threshold ZCR : %0.1f\n Threshold Energy(Avg Noise*%d) : %0.5lf\n ", noiseEnergy, initNoiseZCR, thresholdZCR, thresholdNoiseToEnergyFactor, thresholdEnergy);
	}
	
}


/**************************************************************************************
	Detecting Start and End Marker of the Frame.
	Input (global): fp_obsseq_console is to save analysed values in common output file.
					fp_norm_Seg is to save normalised segregated samples into another file.
**************************************************************************************/
void marker_start_end_segregated(){
	bool flag=false;		//to detect start mark
	// -3 to ignore last 3 frames.
	for(int i=0; i<cntTotFrames-initIgnoreEndFrames; ++i){
			if(!flag && avgEnergy[i+1] > thresholdEnergy && avgEnergy[i+2] > thresholdEnergy && avgEnergy[i+3] > thresholdEnergy && avgEnergy[i+4] > thresholdEnergy){ //&& avgEnergy[i+4] > thresholdEnergy && avgEnergy[i+5] > thresholdEnergy
					start = i;
					flag = 1;
			}
			else if(flag && avgEnergy[i+1] <= thresholdEnergy && avgEnergy[i+2] <= thresholdEnergy && avgEnergy[i+3] <= thresholdEnergy && avgEnergy[i+4] <= thresholdEnergy ){ //&& avgEnergy[i+4] < thresholdEnergy && avgEnergy[i+5] < thresholdEnergy
				end = i;
				flag = 0;
				break;
			}
		}
	
	if(flag == 1) end = cntTotFrames - initIgnoreEndFrames; //if end is not found then making the last frame - 3 as the end marker for the word

	long startSample= (start+1)*sizeFrame;
	long endSample= (end+1)*sizeFrame;
	long totFramesVoice = end-start+1;
	
/****************  saving segregated voice data in different file ****************/
	if(fp_obsseq_norm_seg!=NULL && (segregate_speech || segregate_Live_speech)){
		for(long i=startSample; i<endSample; i++){
			fprintf(fp_obsseq_norm_seg, "%lf\n", normSamples[i]);
		}
	}
	

		//printf("\n---- Segregated Data Saved in File: %s ----\n\n", OSeqcompletePathNormSeg);
		//printf(" START Frame : %ld\t END Frame : %ld\t Total Frames : %ld\n", start+1, end+1, totFramesVoice);
		//printf(" Starting Sample : %ld\t Ending Sample : %ld\n", startSample, endSample);
		//printf("\n--------\n");
		//printf("\n------------------------------------------------------------------------\n");

	double starting_time = (1.0*startSample)/(samplingRate);
	double ending_time = (1.0*endSample)/(samplingRate);
	if(fp_obsseq_console!=NULL){
		fprintf(fp_obsseq_console, "\n---- Segregated Data Saved in File: %s ----\n\n", OSeqcompletePathNormSeg);
		fprintf(fp_obsseq_console, " START Frame : %ld\t END Frame : %ld\t Total Frames : %ld\n", start+1, end+1, totFramesVoice);
		fprintf(fp_obsseq_console, " Starting Sample : %ld\t Ending Sample : %ld\n", startSample, endSample);
		fprintf(fp_obsseq_console, " Starting Time (seconds) : %lf\t Ending Time (seconds) : %lf\n", starting_time, ending_time);
		fprintf(fp_obsseq_console, " Max STE Marker Frame : %d\n", STE_Marker);	
		fprintf(fp_obsseq_console, "\n--------\n");
		fprintf(fp_obsseq_console, "\n------------------------------------------------------------------------\n");
	}
}


/**************************************************************************************************
	Reading CodeBook From the File.
**************************************************************************************************/
void read_codebook_from_file(){
	FILE *fp_cb ;			//to read codebook
	char completePathCB[200];
	sprintf(completePathCB, "%s%s", input_folder, codebook_file_name);  // {input_lamda/}+{codebook_file_name}

	fp_cb = fopen(completePathCB, "r");				//to read input codebook
	if(fp_cb == NULL){ 
				perror("\n Error: ");
				printf("\n codebook File Name is: %s", completePathCB);
				getch();
				return;
	}

	//printf("\n -----------------------> Reading CODEBOOK from File: %s <----------------------- \n", completePathCB);	
		//fprintf(fp_console, "\n -----------------------> Reading CODEBOOK from File: %s <----------------------- \n", completePathCB);	

	for (int i = 0; i < M; ++i){
		//printf("M[%d]\t", i+1);
			//fprintf(fp_console, "M[%d]\t", i+1);
		for (int j = 0; j < p; ++j){
				fscanf(fp_cb,"%Lf",&codebook[i][j]);
				//printf("%Lf (%d)\t",codebook[i][j], j+1);
				//fprintf(fp_console, "%Lf (%d)\t",codebook[i][j], j+1);
			}
		//printf("\n");
		//fprintf(fp_console, "\n");
	}

	fflush(fp_cb); fclose(fp_cb); 
}//read_codebook_from_file

/**************************************************************************************************
	Calculating Coefficients Ai's, Ci's for each frame of the segregated voice data.
**************************************************************************************************/
void calculateCoefficientsForFramesOfSpeech(const long totFramesVoice){

	double markedSamples[sizeFrame]={0};

	//long startSample= (start+1)*sizeFrame;
	//long endSample= (end+1)*sizeFrame;
	//long totFramesVoice = end-start+1;

	long int startOfMarkedFrame = (start+1);		

	//int skipCounter=0;
	for(int ff=0; ff<totFramesVoice; ff++){

		long SampleMarkerToCopy= (startOfMarkedFrame + ff)*sizeFrame;	

		for(int i=0; i<sizeFrame; i++){
			markedSamples[i] = normSamples[i + SampleMarkerToCopy];
		} 


	// Calculating R_i values using AutoCorrelation Method.
		double R[p+1] = {0};
			calculateRi_AutoCorrelation(markedSamples, R);
			///* Not Applicable for this assignment as we are taking Max STE Frame.
			if(R[0]==0){
				printf("\n R[0] Energy should not be ZERO, Skipping frame %d, \n", (startOfMarkedFrame + ff));
				system("pause");
				//skipCounter++;
				continue;
			}
			//*/
		// calculating A_i using Durbin algo.
		double A[p+1] = {0};
			calculateAi_Durbin(R, A);

		// calculating Cepstral coefficient.
		double C[q+1] = {0};
			calculateCi_Cepstral(R[0], A, C);

		if(fp_obsseq_console!=NULL && showCoefficientsInFile)
		{
			fprintf(fp_obsseq_console, "\n\t\t\tFile Name: %s, Frame: %d", OSeqfileNameIp, startOfMarkedFrame +ff);
			fprintf(fp_obsseq_console, "\n\nR Coefficient Values \t LPC Coefficient values \t Cepstral Coefficient values \t Raised Sine Window\n");
			fprintf(fp_obsseq_console, "R[%d] = %lf \n",0, R[0]);
		}
		/****************  Final Values of coefficients ****************/ 	
		
		for(int i=1;i<=p;i++){	
				C_rsw[ff][i] = C[i]*w_rsw[i];   
			if(fp_obsseq_console!=NULL && showCoefficientsInFile){	
				fprintf(fp_obsseq_console, "R[%d] = %lf \t ",i, R[i]);
				fprintf(fp_obsseq_console, "A[%d] = %lf \t ", i, A[i]);
				fprintf(fp_obsseq_console, "C[%d] = %lf \t", i, C[i]);	
				fprintf(fp_obsseq_console, "*(%lf)=> C[%d] = %lf \n", w_rsw[i], i, C_rsw[ff][i]);
			}
		}
		/* For This Assignment this never get executed as we assume q=p
		for( i=p+1;i<=q;i++){
			C_rsw[ff][i] = C[i]*w_rsw[i];  
		}
		*/
		
	}//end of ff loop
}//calculateCoefficientsForFramesOfSpeech

/**************************************************************************************************
	Generate Observation Sequence Based on the Ci's of the frame and Minimum Tokhura Distance
**************************************************************************************************/
void generateObeservationSequence(const long totFramesVoice)
{
	
	double temp;
	double td_cb[M] = {0};  //save min distance of Frame Ci's with each Codebook Index;
	int min_index=0;

	for(int ff=0; ff<totFramesVoice; ff++){

		for (int m= 0; m < M; ++m){
				td_cb[m]=0;
			for (int k = 1; k <= p; ++k){
				temp = C_rsw[ff][k] - codebook[m][k-1];
				td_cb[m] += w_tkh[k-1]*temp*temp;
			}// for each k col of codebook
		}// for each m row index of codebook

		min_index=0;
		for (int m= 0; m < M; ++m){
			if(td_cb[m] < td_cb[min_index]){
				min_index = m;
			}
		}// for each m row index of codebook

		OBS_SEQ[ff] = min_index+1;

	}//end of ff loop

}//generateObeservationSequence

/**************************************************************************************************
	Print Observation Sequence in the file
**************************************************************************************************/
void print_observation_sequence(const long totFramesVoice, const int digit){

	fprintf(fp_obsseq_final_op, "-------------------------------- %d:%d -------------------------\n", digit, totFramesVoice);
	for(int ff=0; ff<totFramesVoice; ff++){
		fprintf(fp_obsseq_final_op, "%d ", OBS_SEQ[ff]);
	}
	fprintf(fp_obsseq_final_op, "\n");

}//print_observation_sequence

/**************************************************************************************************
	Generate Observation Sequence For Recordings present in the disk
**************************************************************************************************/
int sequence_generation(unsigned short int seq_type)
{

	unsigned short NumOfFiles;
	char OSeqfilenameIpformatSpecifier[50];
	char filePathInputVoice[50];
	char lambda_obs_seq_file_name[50];
	char seq_type_name[20];
	char obs_seq_name[50]="OBSERVATION SEQUENCE FOR";

	if(seq_type == 1)
	{	
		NumOfFiles = On;
		strcpy(filePathInputVoice, filePathInputVoiceTraining);
		strcpy(OSeqfilenameIpformatSpecifier, "training_%s_%s_%s%d");
		strcpy(lambda_obs_seq_file_name, LambdaFileNames[3]);
		strcpy(seq_type_name, "TRAINING");
	}
	else if(seq_type == 2)
	{	
		NumOfFiles = Ot;
		strcpy(filePathInputVoice, filePathInputVoiceTesting);
		strcpy(OSeqfilenameIpformatSpecifier, "testing_%s_%s_%s%d");
		strcpy(lambda_obs_seq_file_name, LambdaFileNames[4]);
		strcpy(seq_type_name, "TESTING");
	}

	if(codebook_universe_generation)
	{
		//to save Cepstral Coefficients 
		sprintf(OSeqcompletePathFinOp, "%sUniverse.csv", input_folder);  		
		fp_obsseq_final_op = fopen(OSeqcompletePathFinOp, "w"); 
		if(fp_obsseq_final_op==NULL){ 
			perror("\n Error: ");
			printf("\n File Name : \n  %s\n", OSeqcompletePathFinOp);
			getch();
			return 1;
		}
		strcpy(obs_seq_name, "Cepstral Coefficients of Frames in");
	}

	for(int d = 0 ; d<totDigits ; d++) //iterating through all digits. totDigits
	{			
		
		//if(codebook_universe_generation)
		//{
			printf("\n\n\t ---#---#---#---#---#--- GENERATING %s (%s): %s %s ---#---#---#---#---#---\n",obs_seq_name, seq_type_name, WordFolderName,  WordNames[d]);
		//}
		//else
		//{
		//	printf("\n\t ---#---#---#---#---#--- GENERATING OBSERVATION SEQUENCE FOR (%s): %s %s ---#---#---#---#---#---\n", seq_type_name, WordFolderName,  WordNames[d]);
		//}
		for(int fileCounter=1; fileCounter <= NumOfFiles ; ++fileCounter)//iterating through all files of given digits (1 to X).
		{
		/**************** Creating necessary file Path for data. ****************/
			
			//input file name
			sprintf(OSeqcompletePathIp, "%s%s/%s%d.txt", filePathInputVoice, WordNames[d], voice_data_prefix, fileCounter); // filePathInputVoiceTraining/ {0} + / + {obs_} + {1}.txt
			//segregated file data name
			sprintf(OSeqfileNameIp, OSeqfilenameIpformatSpecifier, WordFolderName, WordNames[d], voice_data_prefix, fileCounter); //OSeqfileNameIp = {Digit} +_+ {0} +_+ {obs_} + {1}
			sprintf(OSeqcompletePathNorm, "%s%s_normalized_samples.txt", fileOutputRecordingNorSeg, OSeqfileNameIp); // fileOutputRecordingNorSeg/ {OSeqfileNameIp}_normalized_samples.txt
			sprintf(OSeqcompletePathNormSeg, "%s%s_normalized_segregated_data.txt", fileOutputRecordingNorSeg, OSeqfileNameIp); // fileOutputRecordingNorSeg/ {OSeqfileNameIp} + _normalized_segregated_data.txt
			//to save analysis file
			sprintf(OSeqcompletePathConsole, "%s%s_analysis.txt", fileOutputRecordingAnalysis, OSeqfileNameIp);  // {fileOutputRecordingAnalysis/}+{OSeqfileNameIp}+"_analysis.txt 
			/**************** Opening respective files. ****************/
			fp_obs_seq_ip = fopen(OSeqcompletePathIp, "r");				//to read input file
			fp_obs_seq_norm = fopen(OSeqcompletePathNorm, "w+");		//to save normalised samples
			fp_obsseq_norm_seg = fopen(OSeqcompletePathNormSeg, "w");  //to save segregated recording from start to end
			fp_obsseq_console = fopen(OSeqcompletePathConsole, "w");	// to save analysis data of each file
			
			if(fileCounter==1){
				DisplayCommonSettingsObsSeq(fp_obsseq_console);
				if(!codebook_universe_generation)
				{
					//to save observation sequence 
					sprintf(OSeqcompletePathFinOp, "%s%s %s/%s%s.txt", input_folder, WordFolderName, WordNames[d], lambda_obs_seq_file_name, WordNames[d]);  // {input_lamda/}+{WordFolderName}+" "+{1}+"/"+{obs_seq__}+{1}+".txt 
									
					fp_obsseq_final_op = fopen(OSeqcompletePathFinOp, "w"); //to save compelete observation sequence in one file
				}
			}
			if(fp_obs_seq_ip == NULL || fp_obs_seq_norm == NULL || fp_obsseq_norm_seg == NULL ||  fp_obsseq_console==NULL ){ 
						perror("\n Error: ");
						printf("\n File Names are : \n  %s, \n  %s, \n  %s, \n %s \n", OSeqcompletePathIp, OSeqcompletePathNorm, OSeqcompletePathNormSeg, OSeqcompletePathConsole  );
						getch();
						return 1;
			}
			
			
			
		if(fileCounter==1){
			printf("  ---->  FILE: %s,\n", OSeqcompletePathIp);  
		}
		else
		{
			printf("\t %s%d.txt,", voice_data_prefix, fileCounter);  
		}


		fprintf(fp_obsseq_console, "\n ----------------------- START - ANALYZING OF FILE: %s ----------------------- \n", OSeqcompletePathIp);

		/**************** DC Shift and Normalizing ****************/
			normalize_dcshift_samples();

		/**************** Frames ZCR and Energy. STE Marker ****************/
			zcr_energy_frames();

		   //if(segregate_speech){						//only if you want to segregate speech into separate file.
			/****************  calculating noise energy and threshold. ****************/
				noiseEnergy_thresholds_frames();						// if you want to calculate thresholds for zcr and energy
					
			/**************** start and end marker of speech ****************/
				marker_start_end_segregated();							//this and above func, if you want to detect start, end marker of speech, and to save it in separate file.
				fclose(fp_obsseq_norm_seg);	// closing file stream
			//}
		 //  else
		 //  {
			//   fclose(fp_obsseq_norm_seg);		// closing file stream
			//   remove(OSeqcompletePathNormSeg);		//removing unnecessory file created.
		 //  }
			if(!segregate_speech)
			{
				remove(OSeqcompletePathNormSeg);		//removing unnecessory
			}

		  // closing file stream, as no longer needed.
		   fflush(fp_obs_seq_ip); fclose(fp_obs_seq_ip); 
		   fflush(fp_obs_seq_norm); fclose(fp_obs_seq_norm);
		   remove(OSeqcompletePathNorm);	//comment it if you want to keep normalised data file.

		/****************  Calculating Coefficients for Voiced Frames of File ****************/
			long totFramesVoice = end-start+1;
			calculateCoefficientsForFramesOfSpeech(totFramesVoice); //for each frame calculate coefficients
			
			if(codebook_universe_generation)
			{
				for(int ff=0; ff<totFramesVoice; ff++){
					for(int i=1;i<=p;i++){
						fprintf(fp_obsseq_final_op, "%lf,", C_rsw[ff][i]);
					}
					fprintf(fp_obsseq_final_op, "\n");
				}
			}
			else
			{
				generateObeservationSequence(totFramesVoice);	//for each frame calculate codebook index
				print_observation_sequence(totFramesVoice,  d);		//print observation seq in file
			}

				//printf("\n ----------------------- END Analyzing OF File: %s ----------------------- \n", OSeqfileNameIp);  
				fprintf(fp_obsseq_console, "\n ----------------------- END - ANALYZING OF FILE: %s ----------------------- \n", OSeqfileNameIp);
			
			// closing  stream, as no longer needed.
			fflush(fp_obsseq_console); fclose(fp_obsseq_console);
		}//end of filecounter loop -------------------------------------------------------------------------------------------------------------------
	
	if(!codebook_universe_generation)
	{
		printf("\n\n  ----> Digit %s:: \n\t Observation Sequence File Generated: %s\n\n", WordNames[d], OSeqcompletePathFinOp); 
		printf("\n-----------------------------------------------------\n");
			
		// closing file stream, as no longer needed.
		fflush(fp_obsseq_final_op); fclose(fp_obsseq_final_op);
	}
		//system("pause");
	}//end of digit loop ------------------------------------------------------------------------------------------------------------------------------
	
	if(codebook_universe_generation)
	{
		printf("\n\n\n  ----> CodeBook Universe File Generated: %s\n\n", OSeqcompletePathFinOp); 
		printf("\n-----------------------------------------------------\n");
		fflush(fp_obsseq_final_op); fclose(fp_obsseq_final_op);
	}
	
	return 0;
}//sequence_generation

/**************************************************************************************************
	Generate Observation Sequence For one Live Recording
**************************************************************************************************/
int live_sequence_generation()
{	
	printf("\n-----------------------------------------------------\n");
	printf("\n\t ---~---~---~---~---~--- GENERATING OBSERVATION SEQUENCE FOR (%s) ---~---~---~---~---~---\n", "LIVE VOICE RECORDING");
	/**************** Creating necessary file Path for data. ****************/
			
	//input file name
	sprintf(OSeqcompletePathIp, "%s%s.txt", input_live_voice_folder, liveRecordingFileName); 
	//segregated file data name 
	sprintf(OSeqcompletePathNorm, "%s%s_normalized_samples.txt", fileOutputRecordingNorSeg, liveRecordingFileName); 
	sprintf(OSeqcompletePathNormSeg, "%s%s_normalized_segregated_data.txt", fileOutputRecordingNorSeg, liveRecordingFileName); 
	//to save analysis file
	sprintf(OSeqcompletePathConsole, "%s%s_analysis.txt", fileOutputRecordingAnalysis, liveRecordingFileName);
	
	/**************** Opening respective files. ****************/
	fp_obs_seq_ip = fopen(OSeqcompletePathIp, "r");				//to read input file
	fp_obs_seq_norm = fopen(OSeqcompletePathNorm, "w+");		//to save normalised samples
	fp_obsseq_norm_seg = fopen(OSeqcompletePathNormSeg, "w");  //to save segregated recording from start to end
	fp_obsseq_console = fopen(OSeqcompletePathConsole, "w");	// to save analysis data of each file

	DisplayCommonSettingsObsSeq(fp_obsseq_console);
	//to save observation sequence 
	sprintf(OSeqcompletePathFinOp, "%s%s_obs_seq_.txt", input_live_voice_folder, liveRecordingFileName);  //		
	fp_obsseq_final_op = fopen(OSeqcompletePathFinOp, "w"); //to save compelete observation sequence in one file

	if(fp_obs_seq_ip == NULL || fp_obs_seq_norm == NULL ||  fp_obsseq_final_op==NULL || fp_obsseq_norm_seg == NULL ||  fp_obsseq_console==NULL ){ 
			perror("\n Error: ");
			printf("\n File Names are : \n  %s, \n  %s, \n  %s, \n  %s \n %s \n", OSeqcompletePathIp, OSeqcompletePathNorm, OSeqcompletePathNormSeg,  OSeqcompletePathFinOp, OSeqcompletePathConsole  );
			getch();
			return EXIT_FAILURE;
	}
			
	printf("\n  ----> ANALYZING OF FILE: %s", OSeqcompletePathIp);   
	
	fprintf(fp_obsseq_console, "\n ----------------------- START - ANALYZING OF FILE: %s ----------------------- \n", OSeqcompletePathIp);

	/**************** DC Shift and Normalizing ****************/
		normalize_dcshift_samples();

	/**************** Frames ZCR and Energy. STE Marker ****************/
		zcr_energy_frames();

		//if(segregate_Live_speech){						//only if you want to segregate speech into separate file.
		/****************  calculating noise energy and threshold. ****************/
			noiseEnergy_thresholds_frames();						// if you want to calculate thresholds for zcr and energy
					
		/**************** start and end marker of speech ****************/
			marker_start_end_segregated();							//this and above func, if you want to detect start, end marker of speech, and to save it in separate file.
			fclose(fp_obsseq_norm_seg);	// closing file stream
		//}
		//else
		//{
		//	fclose(fp_obsseq_norm_seg);		// closing file stream
		//	remove(OSeqcompletePathNormSeg);		//removing unnecessory file created.
		//}
		if(!segregate_Live_speech)
		{
			remove(OSeqcompletePathNormSeg);		//removing unnecessory
		}

	 // closing file stream, as no longer needed.
		fflush(fp_obs_seq_ip); fclose(fp_obs_seq_ip); 
		fflush(fp_obs_seq_norm); fclose(fp_obs_seq_norm);
		remove(OSeqcompletePathNorm);	//comment it if you want to keep normalised data file.

		/****************  Calculating Coefficients for Voiced Frames of File ****************/
			long totFramesVoice = end-start+1;
			calculateCoefficientsForFramesOfSpeech(totFramesVoice); //for each frame calculate coefficients
			generateObeservationSequence(totFramesVoice);	//for each frame calculate codebook index
			print_observation_sequence(totFramesVoice,  -1);		//print observation seq in file
				
				//printf("\n ----------------------- END Analyzing OF File: %s ----------------------- \n", OSeqfileNameIp);  
				fprintf(fp_obsseq_console, "\n ----------------------- END - ANALYZING OF FILE: %s ----------------------- \n", OSeqfileNameIp);
			
			// closing  stream, as no longer needed.
			fflush(fp_obsseq_console); fclose(fp_obsseq_console);

	
	
		printf("\n\n  ----> Live Observation Sequence File Generated: %s\n\n", OSeqcompletePathFinOp); 
		printf("\n-----------------------------------------------------\n");
			
		// closing file stream, as no longer needed.
		fflush(fp_obsseq_final_op); fclose(fp_obsseq_final_op);
		//system("pause");

	return 0;
}//LiveSequence
