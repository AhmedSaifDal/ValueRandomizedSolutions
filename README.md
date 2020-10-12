# ValueRandomizedSolutions

This is a repositery for the codes used and results presented in the numerical testing section of the paper titled "The Value of Randomized Solutions in Mixed-Integer Distributionally Robust Optimization Problems" by Erick Delage and Ahmed Saif. Below is a description of these codes and results.

1. "DRAPr.m" is a Matlab function that solves the Distributionally Robust Assignment Problem (DRAP) with a Wasserstein distributional ambiguity set and a polyhedral support set, WITH and WITHOUT randomization.

2. "DRUFLPr.m" is a Matlab function that solves the Distributionally Robust Uncapacitated Facility Location Problem (DRUFLP) with a Wasserstein distributional ambiguity set and a polyhedral support set, WITH and WITHOUT randomization.

3. "DRCFLP.m" is a Matlab function that solves the Distributionally Robust Capacitated Facility Location Problem (DRCFLP) with a Wasserstein distributional ambiguity set and a polyhedral support set, WITHOUT randomization.

4. "DRCFLPr.m" is a Matlab function that solves the Distributionally Robust Uncapacitated Facility Location Problem (DRCFLP) with a Wasserstein distributional ambiguity set and a polyhedral support set, WITH randomization.

5. "DRAPr_solver.m" is a Matlab code for conducting the DRAP numerical experiments described in Section 7.1. It performs the following tasks:

        a. Generates random instances of the DRAP.
  
        b. Solves the deterministic and randomized strategy DRAP by calling the function DRAPr.
  
        c. Tests the solutions obrained on out-of-sample data.
  
        d. Writes the results to an Excel file.
  
6. "DRUFLPr_solver.m" is a Matlab code for conducting the DRUFLP numerical experiments described in Section 7.2. It performs the following tasks:

        a. Generates random instances of the DRUFLP.
  
        b. Solves the deterministic and randomized strategy DRUFLP by calling the function DRUFLPr.
  
        c. Tests the solutions obrained on out-of-sample data.
  
        d. Writes the results to an Excel file.
        
        
7. "DRCFLPr_solver.m" is a Matlab code for conducting the DRCFLP numerical experiments described in Section 7.3. It performs the following tasks:

        a. Generates random instances of the DRCFLP.
  
        b. Solves the deterministic strategy DRCFLP, and its linear relaxation, by calling the function DRUFLP.
        
        c. Solves the randomized strategy DRCFLP by calling the function DRUFLPr.
  
        d. Tests the solutions obrained on out-of-sample data.
  
        e. Writes the results to an Excel file.
        
 
8. "Value of Randomization DRAP_100_10.xlsx" is an Excel spreadsheet that contains the results used to produce Figure 2 in the paper (DRAP results).

9. "Value of Randomization UFLP_300_10.xlsx" is an Excel spreadsheet that contains the results used to produce Figure 3 in the paper (DRUFLP results).

10. "Value of Randomization CFLP_20_3_10.xlsx" and "Value of Randomization CFLP_20_5_10.xlsx" are Excel spreadsheets that contains the results used to produce Figure 4 in the paper (DRCFLP results).
