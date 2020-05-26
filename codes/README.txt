				Supplementary Material and Code of 
"Triplet Matching for Estimating Causal Effects with Three Treatment Arms: 
	A Comparative Study of Mortality by Trauma Center Level"


We provide as online supplement:

1) The document "Supplementary Material.pdf", which provides proofs and results that 
	did not fit into the main body of the paper.
	
2) The folder "code", which provides the code used for the paper. 
	The folder has two subfolders:
	
	a) "simulation analysis" contains the code used for the simulation analysis.
		In particular:
		- Functions implementing the matching algorithm and to generate the samples are contained respectively in "lib/functionsMatching_7.R" and "lib/functionsDataGeneration_5.R".
		- The simulations have been run with the code "run simulations.R".
		- The output of the simulations are stored in the folder "result simulations"
		- The result of the simulations are analyzed and interpreted with the code in "analysis result simulations.R"

	b) "analysis NEDS" contains the code used for the application of the 3-way matching to the NEDS dataset.
		NOTE: the code in this folder is not replicable since we are not allowed to distribute the NEDS data. 
		The code is organized as follows:
		- "matching.R" contains the application of the matching algorithm to the NEDS dataset.
		- "standardized differences.R" evaluates the standardized differences of the baseline variables in the matched sample.
		- "outcomeAnalysis.Rmd" and "outcomeAnalysis.pdf" provide the code and the results, respectively, of the outcome analysis in the matched sample, including the sensitivity analysis to hidden bias.