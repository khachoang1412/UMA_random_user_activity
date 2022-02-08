# Unourced Multiple Access With Random User Activity

This repository contains the matlab numerical routines of the paper:

[1] K.-H. Ngo, A. Lancho, G. Durisi and A. Graell i Amat, "Unsourced Multiple Access With Random User Activity," submitted to IEEE Trans. Inf. Theory, Jan. 2022.

A short version of this work was reported in

[2] K.-H. Ngo, A. Lancho, G. Durisi and A. Graell i Amat, "Massive Uncoordinated Access With Random User Activity," in Proc. IEEE Int. Symp. Inf. Theory (ISIT), Melbourne, Australia, Jul. 2021, pp. 3014-3019. 

Please, cite the aforementioned papers if you use this code.

## Content of the repository

This repository contains the codes to evaluate random-coding bounds on the error probability of an unsourced multiple access (UMA) system, and the required energy per bit (EbN0) for such probability to lie below a certain threshold. There are four folders:

1. `/RCU_KaKnown/`: RCU bound on the per-user probability of error (PUPE) for an UMA system where the number of active users, denoted by Ka, is known.
  - `RCU_KaFixedKnown.m`: RCU bound for the case where Ka is fixed and known. This is a straightforward extension of the bound in the following seminal paper to the complex-valued case: 

  [3] Y. Polyanskiy, "A perspective on massive random-access," in Proc. IEEE Int. Symp. Inf. Theory (ISIT), 2017, Aachen, Germany, pp. 2523-2527.
  
  To evaluate the RCU bound in [3], simply adjust the framelength, i.e., the degrees of freedom.
  
  - `RCU_KaRandomKnown.m`: RCU bound for the case where Ka is random (follows a distribution with specified PMF) and known. This bound is obtained by averaging the bound above for fixed and known Ka over the distribution of Ka. 
  - `RCU_KaPoissonKnown.m`: RCU bound for the case where Ka follows a Poisson distribution and is known.
  - `EbN0_KaPoissonKnown.m`: Find the minimal required EbN0 such that the RCU bound for Ka following a Poission distribution and known lie below a certain threshold.
  - `binary_search.m`, `golden_search.m`: Auxiliary functions.
  
2. `/RCU_KaUnknown/`: RCU bounds proposed in [1], [2] on the misdetection (MD) and false alarm (FA) probabilities for an UMA system where Ka is random and unknown.
  - `RCU_KaRandomUnknown.m`: RCU bounds for the case where Ka follows a distribution with specified PMF ([1, Theorem 1]). 
  - `RCU_KaPoissonUnknown.m`: RCU bounds for the case where Ka follows a Poisson distribution.
  - `RCU_floor_KaRandomUnknown.m`: Error floors characterized in [1, Theorem 3].
  - `EbN0_KaPoissonUnknown.m`: Find the minimal required EbN0 such that the RCU bounds on the MD and FA probabilities lie below certain thresholds.
  - `binary_search_P_MDFA.m`, `golden_search_P1_MDFA`: Auxiliary functions.

3. `/SA-MPR/`: Application of the bounds above for slotted ALOHA (SA) with multi-packet reception (MPR).
  - `EbN0_SAMPR_KaPoissonKnown.m`: Find the minimal required EbN0 for SA-MPR, Ka follows a Poisson distribution and is known.
  - `EbN0_SAMPR_KaPoissonUnknown.m`: Find the minimal required EbN0 for SA-MPR, Ka follows a Poisson distribution and is unknown. See [1, Corollary 2].
 
4. `/TIN/`: RCUs bound for a scheme that treats interference as noise (TIN).
  - `EbN0_TIN_KaPoissonUnknown.m`: Find the minimal required EbN0 for TIN, Ka follows a Poisson distribution and is unknown.

**Note**: the codes for the SPARC and enhance SPARC schemes are provided by the group of Krishna R. Narayanan, and hence not included in this repository. 
