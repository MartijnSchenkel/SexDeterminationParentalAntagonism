The data in this folder are used to generate Figure 4 from the manuscript. Below is a brief 
description of each data set. Further parameter details are included in the manuscript. Each
data set represents the summarized data from 5,000 simulations with specific parameter values
being randomized for each simulation. The first line in each file contains the variable names, 
which are as follows (SA = sexually-antagonistic selection only, PA = parentally-antagonistic 
selection only)

FAMat 		= Frequency of s_22 on the maternal copy in females
FAPat		= Frequency of s_22 on the paternal copy in females
FXYMat		= Frequency of s_12 on the maternal copy in females
FXYPat		= Frequency of s_12 on the paternal copy in females
MAMat		= Frequency of s_22 on the maternal copy in males
MAPat		= Frequency of s_22 on the paternal copy in males
MXYMat		= Frequency of s_12 on the maternal copy in males
MXYPat		= Frequency of s_12 on the paternal copy in males
SPY		= Value of u_1 (selective effect of p_12 in homozygotes)
SPA		= Value of u_2 (selective effect of p_22 in homozygotes)
AY		= Value of a_1 (fitness modifier of p_12 in p_12/p_11 heterozygotes) (PA)
AA		= Value of a_2 (fitness modifier of p_22 in p_22/p_21 heterozygotes) (PA)
BY		= Value of b_1 (fitness modifier of p_12 in p_22/p_21 heterozygotes) (PA)
BA		= Value of b_2 (fitness modifier of p_12 in p_21/p_22 heterozygotes) (PA)
RY		= Value of r_1 (recombination between S_1 and P_1)
RA		= Value of r_2 (recombination between S_2 and P_2)
CY		= Value of c_1 (scaling parameter for b_1) (PA)
CA		= Value of c_2 (scaling parameter for b_2) (PA)
HFY		= Value of h_F_1 (dominance of p_12 in female heterozygotes) (SA)
HMY		= Value of h_M_1 (dominance of p_12 in male heterozygotes) (SA)
HFA		= Value of h_F_2 (dominance of p_22 in female heterozygotes) (SA)
HMA		= Value of h_M_2 (dominance of p_22 in male heterozygotes) (SA)
SA_XY		= Logical variable denoting if P1 is under sexual (1) or parental (0) antagonism
SA_A		= Logical variable denoting if P2 is under sexual (1) or parental (0) antagonism
MaleDeterminer	= Logical variable denoting if s_22 is male- (1) or female-determining (0)

./results_Fem_PS_combined
  Data from simulations with parentally-antagonistic selection on P1 and sexually-antagonistic 
  selection on P2; s_22 represents a dominant feminizing allele.

./results_Fem_SP_combined
  Data from simulations with sexually-antagonistic selection on P1 and parentally-antagonistic 
  selection on P2; s_22 represents a dominant feminizing allele.

./results_Mas_PS_combined
  Data from simulations with parentally-antagonistic selection on P1 and sexually-antagonistic 
  selection on P2; s_22 represents a functionally-homologous masculinizing allele.

./results_Mas_SP_combined
  Data from simulations with sexually-antagonistic selection on P1 and parentally-antagonistic 
  selection on P2; s_22 represents a functionally-homologous masculinizing allele.