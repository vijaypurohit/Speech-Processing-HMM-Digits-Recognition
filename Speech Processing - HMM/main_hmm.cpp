// Speech Processing - HMM.cpp : Defines the entry point for the console application.
/******************************************************************************************************************************************
	Speech Processing CS 566: Assignment 05 (HMM)
	Roll No: 214101058 MTech CSE'23 IITG
*********************************************************************************************************************************************/

#include "stdafx.h"
#pragma warning (disable : 4996)		//to disable fopen() warnings
#include <stdio.h>
#include <stdlib.h>		//atoi, atof
#include <string.h>		//strcpy
#include <conio.h>		//getch,
#include <time.h>		// time(NULL)
#define _USE_MATH_DEFINES	// for pi
#include <math.h>       // abs, log, sin, floor
#include <windows.h>
#pragma comment(lib,"Winmm.lib")	// lib for listenning to sound

/******************--------Common Settings--------******************/
#define p 12						 //	Number of Cepstral Coefficients, dimensionality of vector x
//#define codeBkSize 32		//change // CodeBook Y Size = K = M		
//#define distr_delta 0.000001			 // 0.00001 or 0.000001 abs(old_distortion - new_dist) > delta. 
#define W 10				//change // Number of Words/Digits//HMM
#define N 5					//change // Number of States per HMM
#define M 32				//change // Number of Distinct Observation Symbols
#define Tmax 150			//change // Max Length of Observation Sequence
#define On 20				//change // Number of Training Observations Given in one file
#define Ot 20				//change // Number of Testing Observations Given in one file
#define model_iterations 200				//change // Number of times Re-estimation has to be done
#define repeatConvergence 2

#define TRAINING 1
#define TESTING 2
#define InputFolderModel 1
#define OutputFolderModel 2

//------------------------------------------------- DEBUG PURPOSE ----------------------------------------------------------------
bool segregate_speech = true;			// True: to segreagate speech data with respect to start and end marker in its output folder. 
bool segregate_Live_speech = true;
bool showCoefficientsInFile = false;		// True: Show Coefficients Values R, A, C's in each voice analysised files
bool showAlphaBetaPstarInConsole = false;	//True: to show alpa, beta probabilities in the console for each observation sequence
bool showStateSeqAlphaBetaInFileForEachObsAfterConverge = true;

//----------------------------------------------------------------  ----------------------------------------------------------------
/******************--------Variable Decalaration--------******************/
// HMM Lamda
long double PI[N];			// Initial State Distribution
long double A[N][N];		// State Transition Probability Distribution
long double B[N][M];		// Observation Symbol Probability Distribution						// 
int O[On>Ot?On:Ot][Tmax];			// N observation Sequence Length of T (for training, testing), Observation Sequence Given, Codebook Index Entry
int OFmax[On>Ot?On:Ot];				// max frames in each observation sequence

int offline_correct_count=0;	// to detect number of utterance of digit correctly recognized
int offline_overall_count=0;	// to detect number of utterance of digit correctly recognized in all files

// P1 Evaluation Problem (Scoring Problem)
long double alpha[Tmax][N]; // T*N, forward procedure variable, Prob of Partial Observation Seq O1 to Ot until time t and state Si  at time t given model lambda
long double beta[Tmax][N];  // T*N, backward procedure variable, Prob of Partial Obser Seq O(t+1) to Ot given state Si at time t and mdodel lamda
long double probability_alpha, probability_beta;

// P2 Uncovering the Problem
long double gamma[Tmax][N]; // T*N, Prob of being in state Si at time t given the Obser Seq and model lamda
long double delta[Tmax][N]; // T*N, Viterbi Algo variable, best score along a single path, at time t which accounts for the first t obser and ends in state Si
int psi[Tmax][N];   // T, Viterbi Algo variable, Book Keeping Var to keep track of argument that maximized the delta
long double Pstar, Pstar_old;	 // Viterbi Algo variable, max probability of delta_T_i 
int Qstar[Tmax];	 // Viterbi Algo variable, State Sequence path

// P3 Reestimation Problem
long double XI[Tmax][N][N];		// Xai Matrix
long double PI_BAR[N];			// Re-estimated Initial State Distribution
long double A_BAR[N][N];		// Re-estimated State Transition Probability Distribution
long double B_BAR[N][M];		// Re-estimated Observation Symbol Probability Distribution
long double converged_A[On][N][N];		// save all the converged A matrix of training Sequence
long double converged_B[On][N][M];		// save all the converged B matrix of training sequence			
long double A_Prev[N][N];		// Previous Averaged A Matrix
long double B_Prev[N][M];		// Previous Averaged B Matrix

// CodeBook
long double codebook[M][p];
const char codebook_file_name[] = "codebook_fb.txt";
bool codebook_universe_generation = false;
// Files								
const char input_folder[] = "input_lamda/";
const char output_folder[] = "output/";
const char WordFolderName[] = "Digit";
const char *LambdaFileNames[] = {"A_","B_","Pi_","obs_seq_training_", "obs_seq_testing_"};
const char *WordNames[]={"0","1","2","3","4","5","6","7","8","9"};
const char voice_data_prefix[] = "obs_";
const char output_folder_Model_name[] = "Models";

FILE *fp_console ;		//to write output file		
char completePathOuput[200];

//Live Voice
time_t timestamp;				//timestamp for live voice filename
char liveRecordingCommand[300], liveRecordingFileName[100];
const char recording_module_exe_path[] = "RecordingModule\\Recording_Module.exe";
const char input_live_voice_folder[] = "input_live_voice_data/";
#define liveRecordingWAV "input_live_voice_data/live_input.wav"


// Observation Sequence Generations Functions	
#include "observation_sequence.h"

// hmm functions	
#include "hmm_functions.h"

/**************************************************************************************
	To Display Common Settings used in our System
	Input: File Pointer in case things needed to be written on file.
**************************************************************************************/
void DisplayCommonSettings(FILE *fp_set=NULL){
	// General Information to Display
	if(fp_set==NULL){
		printf("****-------- WELCOME TO HMM --------****\n");		
		printf("-Common Settings are : -\n");	
		printf(" P (=Q)(#of Cepstral Coefficients) : %d\n", p);
		printf(" Number of Words/Digits/HMM (W) : %d\n", W);	
		printf(" Number of States per HMM (N) : %d\n", N);	
		printf(" Number of Distinct Observation Symbols (M) or CodeBook Size (Y) : %d\n", M); 	
		printf(" Max Length of Observation Sequence (T) : %d\n", Tmax);			
		printf(" Number of Training Observations : %d\n", On);	
		printf(" Number of Testing Observations : %d\n", Ot);

		printf("\n");
		printf(" Frame Size : %d\n", sizeFrame);	
		printf(" Tokhura Weights : ");
		for(int i=0; i<q; i++){
			printf("%0.1f(%d) ", w_tkh[i],i+1);
		}
		printf("\n Amplitutde Value to Scale : %d\n", scaleAmp);			
		printf(" Intital Header Lines Ignore Count : %d\n", initIgnoreHeaderLines); 
		printf(" Intital Samples to Ignore : %d\n",initIgnoreSamples);	
		printf(" Intital Noise Frames Count : %d\n",initNoiseFrames);	
		printf(" Noise to Energy Factor : %d\n",thresholdNoiseToEnergyFactor); 
		printf(" Sampling Rate of Recording: %d\n",samplingRate); 
		printf("----------------------------------------------------------------\n\n");	
	}
	else{
		//printing in file
		fprintf(fp_console,"****-------- WELCOME TO HMM --------****\n");		
		fprintf(fp_console,"-Common Settings are : -\n");	
		fprintf(fp_console," P (=Q)(#of Cepstral Coefficients) : %d\n", p);
		fprintf(fp_console," Number of Words/Digits/HMM (W) : %d\n", W);	
		fprintf(fp_console," Number of States per HMM (N) : %d\n", N);	
		fprintf(fp_console," Number of Distinct Observation Symbols (M) or CodeBook Size (Y) : %d\n", M); 	
		fprintf(fp_console," Max Length of Observation Sequence (T) : %d\n", Tmax);			
		fprintf(fp_console," Number of Training Observations : %d\n", On);	
		fprintf(fp_console," Number of Testing Observations : %d\n", Ot);

		fprintf(fp_console,"\n");
		fprintf(fp_console," Frame Size : %d\n", sizeFrame);	
		fprintf(fp_console," Tokhura Weights : ");
		for(int i=0; i<q; i++){
			fprintf(fp_console,"%0.1f(%d) ", w_tkh[i],i+1);
		}
		fprintf(fp_console,"\n Amplitutde Value to Scale : %d\n", scaleAmp);			
		fprintf(fp_console," Intital Header Lines Ignore Count : %d\n", initIgnoreHeaderLines); 
		fprintf(fp_console," Intital Samples to Ignore : %d\n",initIgnoreSamples);	
		fprintf(fp_console," Intital Noise Frames Count : %d\n",initNoiseFrames);	
		fprintf(fp_console," Noise to Energy Factor : %d\n",thresholdNoiseToEnergyFactor); 
		fprintf(fp_console," Sampling Rate of Recording: %d\n",samplingRate); 
		fprintf(fp_console,"----------------------------------------------------------------\n\n");	

	}
}

/**************************************************************************************
	Main Function
**************************************************************************************/
int _tmain(int argc, _TCHAR* argv[])
{

	 /*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
																		 Intialization
	-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
	CalculateWeightsForRaisedSineWindow();	// calculating weights for Raised sine window before hand using in program.
	read_codebook_from_file();				// read codebook from file.
	//RecordMyVoice();
	//
	//getch();	return 0;
	char choice;		// choice exercised.	 
  do{
		char ch;
		double accuracy=0, final_accuracy=0;
		unsigned int model_type_use = OutputFolderModel;
		char model_type_string[20] = "new_model";
		bool temp=false;
		char recognised, correct_voice;

		printf("\n\n\n");
		system("pause");
		system("cls");

		printf("\n\n ------- -------~~ HMM MENU ~~------- -------");

		printf("\n d. DISPLAY Common System Settings Used.");
		printf("\n 1. OBS SEQ Generation: All Training Files.");
		printf("\n 2. OBS SEQ Generation: All Testing Files.");

		printf("\n\n 3. CONVERGE: Converge Model Lambda for Each Digit. ");

		printf("\n\n 4. TESTING: Offline Testing of Digit: Using Old Input Folder Model");

		printf("\n 5. TESTING: Offline Testing of Digit: Using New Converged Output Folder Model");
		printf("\n 6. REPLACE: OLD MODEL In Default Input Folder with NEW CONVERGED MODEL in Output Folder.");

		printf("\n\n 7. TESTING: Live Testing: Using Old Input Folder Model.");
		printf("\n 8. TESTING: Live Testing: Using New Converged Output Folder Model.");

		printf("\n\n u. CodeBook: Generate Cepstral Coeff Universe File From Training Files");

		printf("\n\n n. Exit - Bye		\n\n  --Choice : ");
		scanf("%c%*c", &ch);
		printf("\n <-------->");
        
		switch (ch) {
			case 'd' : DisplayCommonSettings();
				break;

/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
											"Training Observation" Sequence Generation for Each Digit
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			case '1' : 
				// generate observation sequence from the training files. each digit observation sequence is saved in the input lambda folder.
				// also generate normalised files and segrageted data of voice, and analysis of the voice data.
				sequence_generation(TRAINING);
				break;
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
											"Testing Observation" Sequence Generation for Each Digit
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			case '2' : 
				// generate observation sequence from the testing files. each digit observation sequence is saved in the input lambda folder.
				// also generate normalised files and segrageted data of voice, and analysis of the voice data.
				sequence_generation(TESTING);
				break;
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
											Converge Model Lambda For Each Digit
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			case '3' : 
					for(int d=0; d<W; d++){
						/**************** Creating necessary file Path for data. ****************/
						sprintf(completePathOuput, "%s%s_%s_HMM_Converged_log.txt", output_folder, WordFolderName, WordNames[d]);   
						/**************** Opening respective files. ****************/
						fp_console = fopen(completePathOuput, "w");					//to read input observation sequence
						if(fp_console == NULL){ 
								perror("\n Error: ");
								printf("\n File Names is: %s \n ", completePathOuput);
								getch();
								return EXIT_FAILURE;
						}
							fprintf(fp_console, "\n ----------------------- -----------------------> CONVERGING LAMDA : %s %s <----------------------- -----------------------\n", WordFolderName, WordNames[d]); 	 
							fprintf(fp_console, "\n ------------------------------------------------------------------------\n"); 	 
						/**************** Reading  Obs Seq from File ****************/
							readObsSeq(d, TRAINING);		// for each digit read their training sequence observations
							fprintf(fp_console, "\n ------------------------------------------------------------------------\n"); 	 

						/**************** Making Converged Lambda From Bakis Model ****************/
							covergence_procedure();				// FOR EACH observation seq generate their model and then converge finally by taking average of all.
							output_lambdaABPi_to_each_file(d);
								printf("\n\n\t -------->> New Lambda Files Saved: %s%s/%s/", output_folder, output_folder_Model_name, WordNames[d]);
								printf("\n\n\t -------->> Convergence Done, Log File Generated: %s\n\n", completePathOuput);
									
								fprintf(fp_console, "\n ---------------------------------- --------------------------------------"); 
								fprintf(fp_console, "\n <---------------------------------- ----------------------------------> END <---------------------------------- -------------------------------------->"); 
						 fflush(fp_console); fclose(fp_console); 
					}// for each digit d<W Converge Model Lambda For Each Digit
				break;
			
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
											OFFLINE TESTING USING INPUT FOLDER MODEL
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			case '4' : 

				model_type_use = InputFolderModel;
				strcpy(model_type_string,"old_model");
			
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
											OFFLINE TESTING USING OUTPUT FOLDER MODEL
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			case '5' : 
				offline_overall_count=0;
				for(int d=0; d<W; d++){
					system("cls");
					/**************** Creating necessary file Path for data. ****************/
					sprintf(completePathOuput, "%s%s_%s_HMM_offline_test_result_%s.txt", output_folder, WordFolderName, WordNames[d],model_type_string);  // {output/}+{WordFolderName}+"_"+{1}+".txt"
					/**************** Opening respective files. ****************/
					fp_console = fopen(completePathOuput, "w");					//to read input observation sequence
					if(fp_console == NULL){ 
							perror("\n Error: ");
							printf("\n File Names is: %s \n ", completePathOuput);
							getch();
							return EXIT_FAILURE;
					}
						fprintf(fp_console, "\n ----------------------- ----------------------- > OFFLINE TESTING : %s %s < ----------------------- -----------------------\n", WordFolderName, WordNames[d]); 	 
						fprintf(fp_console, "\n ------------------------------------------------------------------------\n"); 	 
					/**************** Reading  Obs Seq from File ****************/
						readObsSeq(d, TESTING);
						fprintf(fp_console, "\n ------------------------------------------------------------------------\n"); 	

					/**************** Offline Testing ****************/
						offline_testing(d, model_type_use);

						 accuracy = (double)(offline_correct_count*1.0/Ot)*100;
						printf("\n\t FOR %s %s | Accuracy:  %0.2f %%\n\n", WordFolderName, WordNames[d], accuracy); 
							fprintf(fp_console, "\n\t FOR %s %s | Accuracy:  %0.2f %%\n\n", WordFolderName, WordNames[d], accuracy);

						//printf("\n\n\t -------->> New Lambda Files Saved: %s%s/%s/", output_folder, output_folder_Model_name, WordNames[d]);
						printf("\n\n\t -------->> Offline Testing Done, Log File Generated: %s\n\n", completePathOuput);
							fprintf(fp_console, "\n ---------------------------------- --------------------------------------"); 
							fprintf(fp_console, "\n <---------------------------------- ----------------------------------> END <---------------------------------- -------------------------------------->"); 
					fflush(fp_console); fclose(fp_console); 
					system("pause");
					
				}// for each digit d<W  OFFLINE TESTING
					 final_accuracy = (double)(offline_overall_count*1.0/(Ot*W))*100;
						printf("\n\t Overall Accuracy:  %0.2f %%\n\n", final_accuracy); 
						
				break;	
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
										Replace OLD MODEL In Default Input Folder with NEW CONVERGED MODEL in Output Folder. 
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			case '6' : replace_old_models_files(); 
				break;
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
										TESTING: Live Testing: Using Old Input Folder Model.
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			
			case '7' : model_type_use = InputFolderModel;
					   strcpy(model_type_string,"old_model");
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
										TESTING: Live Testing: Using New Converged Output Folder Model
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			
			case '8' : 
					printf("\n Duration 3 sec: \n"); 
					timestamp = time(NULL);
					//timestamp = 30;

					/**************** Creating necessary file Path for data. ****************/
					sprintf(liveRecordingFileName, "live_%ld", timestamp);  //file_name
					sprintf(liveRecordingCommand, "%s 3 %s %s%s.txt", recording_module_exe_path, liveRecordingWAV, input_live_voice_folder, liveRecordingFileName);  
					//printf("\n path: %s\n",liveRecordingCommand );
					/**************** Creating necessary file Path for data. ****************/
						sprintf(completePathOuput, "%s%s_test_result_%s.txt", input_live_voice_folder, liveRecordingFileName ,model_type_string);  
						/**************** Opening respective files. ****************/
						fp_console = fopen(completePathOuput, "w");					//to read input observation sequence
						if(fp_console == NULL){ 
								perror("\n Error: ");
								printf("\n File Names is: %s \n ", completePathOuput);
								getch();
								return EXIT_FAILURE;
						}
					// Right Click Project Name in Solution Explorer
					// Select Propertes --> Linker --> Input
					// Select Additional Dependencies --> Edit
					// Add winmm.lib
					do{
						
						
						do
						{
							/**************** Execute the Live Recording Module ****************/
							system(liveRecordingCommand);		//execute the command
							//USAGE : "Recording_Module.exe" <duration_in_seconds> <output_mono_wav_file_path> <output_text_file_path>
						
							printf("\n Playing Sound: "); 
							PlaySound(TEXT(liveRecordingWAV), NULL, SND_SYNC );
							printf("\n\n Is Word Correctly Spoken (y/n) ?" 
								"\n\t n: will repeat the recording process. "
								"\n\n  --Choice :  ");
							//scanf("%c%*c",&correct_voice);
							scanf("%c",&correct_voice);

							while ((getchar()) != '\n');

						}while(correct_voice!='y');
						
						printf("\n\n");

							fprintf(fp_console, "\n ----------------------- ----------------------- > LIVE TESTING : %s < ----------------------- \n",liveRecordingFileName );
						/**************** Live Testing ****************/
							live_testing(model_type_use);

						printf("\n Is Word Correctly Recognised (y/n/e) ?" 
								"\n n: will repeat the recording process again. "
								"\n e: exit. "
								"\n\n  --Choice :  ");

						scanf("%c%*c",&recognised);
						fprintf(fp_console, "\n Digit Correct Recogntion Status (y/n/e): %c", recognised); 
						fprintf(fp_console, "\n ---------------------------------- --------------------------------------"); 
						fflush(fp_console);
					 }while(recognised == 'n');
					fprintf(fp_console, "\n ---------------------------------- --------------------------------------"); 
							fprintf(fp_console, "\n <---------------------------------- ----------------------------------> END <---------------------------------- -------------------------------------->"); 
						fflush(fp_console); fclose(fp_console); 
				break;
			case 'u' :   
					//generate_codebook_universe(TRAINING);
						temp=segregate_speech;
						segregate_speech=false;
						codebook_universe_generation=true;
							sequence_generation(TRAINING);
						codebook_universe_generation=false;
						segregate_speech=temp;

				break;
			case 'n' :   printf("\n Bye \n");  
				break;
			default  :   printf("\n--Invalid Choice. Enter Again \n");
		}//switch
		choice=ch;
	} while (choice != 'n');

	printf("\n---------------------------------- ENTER TO EXIT --------------------------------------\n");
	getch();
	return 0;
}

//void generate_codebook_universe(unsigned short int seq_type){
//
//	unsigned short NumOfFiles;
//	char OSeqfilenameIpformatSpecifier[50];
//	char filePathInputVoice[50];
//	char lambda_obs_seq_file_name[50];
//	char seq_type_name[20];
//
//	if(seq_type == 1)
//	{	
//		NumOfFiles = On;
//		strcpy(filePathInputVoice, filePathInputVoiceTraining);
//		strcpy(OSeqfilenameIpformatSpecifier, "training_%s_%s_%s%d");
//		strcpy(lambda_obs_seq_file_name, LambdaFileNames[3]);
//		strcpy(seq_type_name, "TRAINING");
//	}
//	else if(seq_type == 2)
//	{	
//		NumOfFiles = Ot;
//		strcpy(filePathInputVoice, filePathInputVoiceTesting);
//		strcpy(OSeqfilenameIpformatSpecifier, "testing_%s_%s_%s%d");
//		strcpy(lambda_obs_seq_file_name, LambdaFileNames[4]);
//		strcpy(seq_type_name, "TESTING");
//	}
//
//	//to save Cepstral Coefficients 
//	sprintf(OSeqcompletePathFinOp, "%sUniverse.csv", input_folder);  		
//	fp_obsseq_final_op = fopen(OSeqcompletePathFinOp, "w"); 
//	if(fp_obsseq_final_op==NULL){ 
//		perror("\n Error: ");
//		printf("\n File Name : \n  %s\n", OSeqcompletePathFinOp);
//		getch();
//		return ;
//	}
//
//	for(int d = 0 ; d<totDigits ; d++) //iterating through all digits. totDigits
//	{			
//		printf("\n\n\t ---#---#---#---#---#--- GENERATING Cepstral Coefficients of Frames in (%s): %s %s ---#---#---#---#---#---\n", seq_type_name, WordFolderName,  WordNames[d]);
//
//		for(int fileCounter=1; fileCounter <= NumOfFiles ; ++fileCounter)//iterating through all files of given digits (1 to X).
//		{
//		/**************** Creating necessary file Path for data. ****************/
//			
//			//input file name
//			sprintf(OSeqcompletePathIp, "%s%s/%s%d.txt", filePathInputVoice, WordNames[d], voice_data_prefix, fileCounter); 
//			//segregated file data name
//			sprintf(OSeqfileNameIp, OSeqfilenameIpformatSpecifier, WordFolderName, WordNames[d], voice_data_prefix, fileCounter); 
//			sprintf(OSeqcompletePathNorm, "%s%s_normalized_samples.txt", fileOutputRecordingNorSeg, OSeqfileNameIp); 
//			sprintf(OSeqcompletePathNormSeg, "%s%s_normalized_segregated_data.txt", fileOutputRecordingNorSeg, OSeqfileNameIp); 
//			//to save analysis file
//			sprintf(OSeqcompletePathConsole, "%s%s_analysis.txt", fileOutputRecordingAnalysis, OSeqfileNameIp);  
//			/**************** Opening respective files. ****************/
//			fp_obs_seq_ip = fopen(OSeqcompletePathIp, "r");				//to read input file
//			fp_obs_seq_norm = fopen(OSeqcompletePathNorm, "w+");		//to save normalised samples
//			fp_obsseq_norm_seg = fopen(OSeqcompletePathNormSeg, "w");  //to save segregated recording from start to end
//			fp_obsseq_console = fopen(OSeqcompletePathConsole, "w");	// to save analysis data of each file
//			if(fileCounter==1){
//				DisplayCommonSettingsObsSeq(fp_obsseq_console);
//			}
//			if(fp_obs_seq_ip == NULL || fp_obs_seq_norm == NULL || fp_obsseq_norm_seg == NULL ||  fp_obsseq_console==NULL ){ 
//					perror("\n Error: ");
//					printf("\n File Names are : \n  %s, \n  %s, \n  %s, \n %s \n", OSeqcompletePathIp, OSeqcompletePathNorm, OSeqcompletePathNormSeg, OSeqcompletePathConsole  );
//					getch();
//					return ;
//			}
//			
//		if(fileCounter==1){
//			printf("  ----> FILE: %s,\n", OSeqcompletePathIp);  
//		}
//		else
//		{
//			printf("\t %s%d.txt,", voice_data_prefix, fileCounter);  
//		}
//
//
//		fprintf(fp_obsseq_console, "\n ----------------------- START - ANALYZING OF FILE: %s ----------------------- \n", OSeqcompletePathIp);
//
//		/**************** DC Shift and Normalizing ****************/
//			normalize_dcshift_samples();
//
//		/**************** Frames ZCR and Energy. STE Marker ****************/
//			zcr_energy_frames();
//
//		   //if(segregate_speech){						//only if you want to segregate speech into separate file.
//			/****************  calculating noise energy and threshold. ****************/
//				noiseEnergy_thresholds_frames();						// if you want to calculate thresholds for zcr and energy
//					
//			/**************** start and end marker of speech ****************/
//				marker_start_end_segregated();							//this and above func, if you want to detect start, end marker of speech, and to save it in separate file.
//				fclose(fp_obsseq_norm_seg);	// closing file stream
//			//}
//		   //else
//		   //{
//			  // fclose(fp_obsseq_norm_seg);		// closing file stream
//			  // remove(OSeqcompletePathNormSeg);		//removing unnecessory file created.
//		   //}
//			if(!segregate_speech)
//			{
//				remove(OSeqcompletePathNormSeg);		//removing unnecessory
//			}
//
//		  // closing file stream, as no longer needed.
//		   fflush(fp_obs_seq_ip); fclose(fp_obs_seq_ip); 
//		   fflush(fp_obs_seq_norm); fclose(fp_obs_seq_norm);
//		   remove(OSeqcompletePathNorm);	//comment it if you want to keep normalised data file.
//
//		/****************  Calculating Coefficients for Voiced Frames of File ****************/
//			long totFramesVoice = end-start+1;
//			calculateCoefficientsForFramesOfSpeech(totFramesVoice); //for each frame calculate coefficients
//
//			for(int ff=0; ff<totFramesVoice; ff++){
//				for(int i=1;i<=p;i++){
//					fprintf(fp_obsseq_final_op, "%lf,", C_rsw[ff][i]);
//				}
//				fprintf(fp_obsseq_final_op, "\n");
//			}
//				
//				//printf("\n ----------------------- END Analyzing OF File: %s ----------------------- \n", OSeqfileNameIp);  
//				fprintf(fp_obsseq_console, "\n ----------------------- END - ANALYZING OF FILE: %s ----------------------- \n", OSeqfileNameIp);
//		
//			fflush(fp_obsseq_console); fclose(fp_obsseq_console);
//		}//end of filecounter loop -------------------------------------------------------------------------------------------------------------------
//		//system("pause");
//	}//end of digit loop ------------------------------------------------------------------------------------------------------------------------------
//		
//	printf("\n\n  ----> CodeBook Universe File Generated: %s\n\n", OSeqcompletePathFinOp); 
//	printf("\n-----------------------------------------------------\n");
//	fflush(fp_obsseq_final_op); fclose(fp_obsseq_final_op);
//}//generate_codebook_universe