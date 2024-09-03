# pdb_rank
Purpose:
These scripts rank PDB models generated from de novo computational protein structure prediction. Most servers internally rank their structures, but it is difficult to compare results between servers. This script uses two benchmarks to rank PDB inputs, 1) Amino acid intensity comparison between the input FASTA sequence of the dimer/multimer and the amino acid sequences gained from the input PDB files and 2) A ranking of the agreeable phi and psi angles from a Ramachandran plot. The scores are combined to yield a final ranking, with the highest ranked PDB file having the highest score. The scores are output in the output text file along with their corresponding file names.
Command: 
python rank_pdbs_2.py [path_to_fasta_file] [path_to_pdb_directory] [path_to_output_file]
Takes input of a fasta file with a single sequence containing all the proteins involved (concatenated)
  Example: 
  >concatenated_seqs
  (combined protein sequences here containing both monomers or the dimer + monomer, case dependent)
Second input is a directory of all the PDB files being used from de novo protein structure prediction servers 
De novo protein prediction PDB source suggestions(multiple sources is best since most sources internally rank their models):
  ZDock (Pierce BG, Wiehe K, Hwang H, Kim BH, Vreven T, Weng Z. (2014) ZDOCK Server: Interactive Docking Prediction of Protein-Protein Complexes and Symmetric Multimers. Bioinformatics 30(12): 1771-3.)

  ClusPro (Desta IT, Porter KA, Xia B, Kozakov D, Vajda S. Performance and Its Limits in Rigid Body Protein-Protein Docking. Structure. 2020 Sep; 28 (9):1071-1081.; Vajda S, Yueh C, Beglov D, Bohnuud T, Mottarella SE, Xia B, Hall DR, Kozakov D. New additions to the ClusPro server motivated by CAPRI. Proteins: Structure, Function, and Bioinformatics. 2017 Mar; 85(3):435-444.; Kozakov D, Hall DR, Xia B, Porter KA, Padhorny D, Yueh C, Beglov D, Vajda S. The ClusPro web server for protein-protein docking. Nature Protocols. 2017 Feb;12(2):255-278.; Kozakov D, Beglov D, Bohnuud T, Mottarella S, Xia B, Hall DR, Vajda, S. How good is automated protein docking? Proteins: Structure, Function, and Bioinformatics. 2013 Dec; 81(12):2159-66.  pdf)

  HDock (Yan Y, Tao H, He J, Huang S-Y.* The HDOCK server for integrated protein-protein docking. Nature Protocols, 2020; doi: https://doi.org/10.1038/s41596-020-0312-x.; Yan Y, Zhang D, Zhou P, Li B, Huang S-Y. HDOCK: a web server for protein-protein and protein-DNA/RNA docking based on a hybrid strategy. Nucleic Acids Res. 2017;45(W1):W365-W373.; Yan Y, Wen Z, Wang X, Huang S-Y. Addressing recent docking challenges: A hybrid strategy to integrate template-based and free protein-protein docking. Proteins 2017;85:497-512.; Huang S-Y, Zou X. A knowledge-based scoring function for protein-RNA interactions derived from a statistical mechanics-based iterative method. Nucleic Acids Res. 2014;42:e55.; Huang S-Y, Zou X. An iterative knowledge-based scoring function for protein-protein recognition. Proteins 2008;72:557-579.)

  GRAMM Docking (Singh, A., Copeland, M.M., Kundrotas, P.J., Vakser, I.A., 2024, GRAMM Web Server for Protein Docking, Methods in Mol. Biol., 2714:101-112.)

Intensity algorithm:
Takes the average amino acid intensity, which is calculated as a function of the position in the sequence. If the amino acid exists at position i, it gets an intensity value, A, that is a function of the sequence position. If it does not exist at that position, then it gets a value of 0. The intensities from all 20 amino acids are combined in order to yield an overall average intensity value for the sequence. These values are compared pairwise to that of the input fasta sequence to get the percent difference for each PDB file. This is subtracted from 100 to give the highest score to the sequence that has the most similar intensity to the input sequence.

Ramachandran:
Not only are the favorable angles factored into the score, but their distance from the favorable regions is also calculated as part of the score. In the allowed regions, the distance from the favorable region is calculated to give a higher score to residues with angles closer to the favorable regions, and the same is done for the disallowed regions with reference to the allowed regions. This distance-based scoring system allows for more precise scoring, allowing differentiation between multiple models from the same source.

The best PDB file will not be output, but instead must be retrieved from the input pdb folder. The names of the PDB files will be output with their corresponding scores and can be visualized using an outside program such as PyMol.

