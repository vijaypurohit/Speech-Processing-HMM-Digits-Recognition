/******************************************************************************************************************************************
	Speech Processing CS 566: Assignment 05 (HMM)
	Roll No: 214101058 MTech CSE'23 IITG
	Input: 
		*	input_lamda/* = Contains: PRESENTLY USED "Lamda Model for Each Digit"
			**	      = Contains: "Universe.csv" of Cepstral Coefficients of Training Files
			**	      = Contains: "Codebook" file to be read.
		*	input_live_voice_data/*	= Contains: Live Recording Files, 
							Their Observation Sequences, 
							Their Test Result Using Model (input or newly converged)
						 (Can Clean These Files Except the Folders Present)
			**	input_live_voice_data/TRAINING	= Contains Digits Recordings generated Using Application for Training Purpose 
			**	input_live_voice_data/TESTING	= Contains: Digits Recordings generated Using Application for Testing Purpose
			*** 	in the end, replace manually
					input_live_voice_data/TRAINING/* --> input_voice_training_data/*
					input_live_voice_data/TESTING/* --> input_voice_testing_data/*
		*	input_voice_training_data/*	= Contains: Training Utterance Recordings for Input into model
		*	input_voice_testing_data/*	= Contains: Test Utterance Recordings for Input into model
		* 	RecordingModule/*		= Contains: Recording Module Files

	Output:  
		*	output/*	= Contains:	Result Analysis of Converged Model for Each Digit.
			**	output/Models/*	= Contains: Newly Generated Model using Input Trainning Files.
		*	output_voice_recordings_analysis_files	= Contains: Recording Analysis files which 
								shows Frames used, Samples used, STE Marker, Cepstral Coefficients etc
								For Files of Input Training, Input Testing, Live Recordings
		*	output_voice_recordings_normalised_segregated	= Contains: Segragated Speech Part using Start and End Marker
								For Files of Input Training, Input Testing, Live Recordings
	Debug Variables:	
		* 	segregate_speech :	True: to segreagate speech data with respect to start and end marker in its output folder (output_voice_recordings_normalised_segregated). 
		*	segregate_Live_speech :	True: to segreagate Live Recording data with respect to start and end marker in its output folder (output_voice_recordings_normalised_segregated). 
		*	showCoefficientsInFile :	True: show Coefficients Values R, A, C's of each frames in its analysised files (output_voice_recordings_analysis_files).
		*	showAlphaBetaPstarInConsole :	True: to show alpa, beta probabilities in the console for each observation sequence. (also saved in files in (output/) )
		*	showStateSeqAlphaBetaInFileForEachObsAfterConverge :	True: It will save each utterance alpha, beta probabilites and state sequence in the file in (output/).

*********************************************************************************************************************************************/
FILES:
	* 	main_hmm.cpp	= Main File Contains Menu	
	* 	hmm_functions.h	= Contains HMM Functions and Solution
	* 	observation_sequence.h	= Observation Sequence Functions

For Error: PlaySound() is not Identified:
Do:
	1 Right Click Project Name in Solution Explorer
	2 Select Propertes --> Linker --> Input
	3 Select Additional Dependencies --> Edit
	4 Add name " winmm.lib "
-------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------
Instructions to execute Code.
-----------------------------
1. Open it in Visual Studio 2010. 
	Main file: main_hmm.cpp 
2. Compile it and Run. Console window will show up.
	Interact With Menu
		Output will be shown on the Console.
		Detailed Output *.txt will be present in their respective folder.
3. Take Care:
	To generate The Respective Sequence (Training/Testing) before Converging or Testing.	
-------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------
THE END.
-----------------------------