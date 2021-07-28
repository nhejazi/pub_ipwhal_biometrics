# Summary of HAL-IPW experiments

* Experimental parameters:
  * DGPs:
      * 1a: no positivity issues (g0: min = 0.31, max = 0.83)
      * 1b: randomized_(g0: min = max = 0.5)
      * 1c: Q0 indep. of A_ (g0: min = 0.29, max = 0.92)
      * 2a: severe positivity issues (g0: min = 0.018, max = 0.77)
      * 2b: moderate positivity issues (g0: min = 0.047, max = 0.69)
      * 3a: no treatment effect, from Ashkan (g0: min = 0.28, max = 0.83)
      * 3b: no treatment effect, from Ashkan (g0: min = 0.18, max = 0.58)
  * estimator: non-stabilized Horvitz-Thompson IPW
  * truncation:
      * profile: choose optimal regularization parameter then truncate within
      * joint: choose optimal combination of regularization and truncation
  * outcome mechanism:
      * Q0: use the true outcome mechanism in the EIF
      * Qn: use the estimated outcome mechanism in the EIF

* Using D_CAR-based selector with estimated outcome Qn ({Mn_max, delta_trunc}):

  |DGP  | HT, profile  | HT, joint   |
  |-----|--------------|-------------|
  |1a   |  {140, 0.30} | {80, 0.20}  |
  |1b   |  {350, 0.50} | {350, 0.50} |
  |1c   |  {100, 0.25} | {150, 0.25} |
  |2a   |  {80, 0.40}  | {80, 0.40}  |
  |2b   |  {120, 0.30} | {120, 0.30} |
  |3a   |  {80, 0.30}  | {80, 0.30}  |
  |3b   |  {400, 0.25} | {400, 0.25} |

* Using score-based selector with estimated outcome Qn:

  |DGP  |    Mn_max    | delta_trunc |
  |-----|--------------|-------------|
  |1a   |      200     |     0.30    |
  |1b   |      500     |     0.50    |
  |1c   |      200     |     0.25    |
  |2a   |      200     |     0.30    |
  |2b   |      200     |     0.30    |
  |3a   |      200     |     0.30    |
  |3b   |      400     |     0.30    |

* Using plateau-based selector with estimated outcome Qn:

  |DGP  |    Mn_max    | delta_trunc |
  |-----|--------------|-------------|
  |1a   |      100     |     0.30    |
  |1b   |      350     |     0.50    |
  |1c   |      NA      |     NA      |
  |2a   |      100     |     0.30    |
  |2b   |      100     |     0.30    |
  |3a   |      100     |     0.30    |
  |3b   |      200     |     0.30    |
