# Sharp2021
Code for implementing iterative acceleration techniques for general systems and within the forward-backward sweep method (FBSM), featured in Sharp et al. 2021:

JA Sharp, K Burrage, MJ Simpson. 2021 Implementation and acceleration of optimal control for systems biology. _Journal of the Royal Society Interface_ **18**, 20210241. (doi.org/10.1098/rsif.2021.0241). 

This repository contains six subfolders corresponding to each of the linear and acute myeloid leukaemia (AML) control problems considered in the work of Sharp et al. 2021, and one subfolder containing the implementation of the fixed point iteration method and acceleration algorithms considered in the work, for a system of arbitrary size N. Where functions can be called by users, we provide example function calls corresponding to results presented in the paper. 

The following system corresponds to Equation S24 of the supplementary material of Sharp et al. 2021, and can be used for example function calls to test the algorithms in the **Nonlinear_System_Solvers** subfolder.   
F = @(X) [X(1)-(X(1)^2+X(2)^2-5)/4;
 X(2)-(X(1)*X(2)-2)/2;
-(X(1)*X(3)-X(2))/3];  
X0 = [0;0;0];

FixedPoint.m - Example function call:  [root,Fevals] = FixedPoint(F,X0,1e-10,100)   
Aitken.m - Example function call:  [root,Fevals] = Aitken(F,X0,1e-10,100,3)     
Anderson.m - Example function call:  [root,Fevals] = Anderson(F,X0,1e-10,100,3,1e10)   
Steffensen.m - Example function call:  [root,Fevals] = Steffensen(F,X0,1e-10,100,3)   
Wegstein.m - Example function call:  [root,Fevals] = Wegstein(F,X0,1e-10,100,3,0,0,0)    

These scripts solve equations of the form _X = F(X)_, for the user specified 'F', initial guess 'X0', convergence tolerance 'tol', and maximum function evaluations 'MaxFevals'. The acceleration methods each require user specified tuning parameters. For the Aitken and Steffensen methods, input 'm' specifies the desired dimension of the _N x m_ difference matrices for the "partial" implementation discussed in Sharp2021; setting _m = N_ corresponds to the standard implementation of the Aitken and Steffensen methods. For Anderson Acceleration, input 'M' determines the maximum number of previous iterations to incorporate in each iteration, while input 'Droptol' specifies the maximum accptable condition number of the residual difference matrix, _dG_. For Wegstein's method, input 'nth' specifies how frequently to update q; every nth iteration, input 'bounding' acts as a switch to turn on or off bounds on _q_, with bounds applied if _bounding=1_. Inputs 'lower' and 'upper' specify the lower and upper bounds to apply when _bounding = 1_. In addition to these tuning parameters, 

Each of the remaining subfolders correspond to a control problem presented in Sharp et al. 2021. The contents of each subfolder and a description are provided below. Note that control problems with fixed endpoints require solving several FBSM problems, and may run for several minutes on a standard desktop computer.  

**Linear_continuous** - Corresponds to the linear continuous control problem described in Section 3 of Sharp et al. 2021. 

This subfolder contains 5 scripts/functions for implementation:
- Linear_base.m - this script can be run directly to solve the linear continuous control problem using the standard FBSM without acceleration
- Linear_Aitken.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear continuous control problem using the FBSM with the parital Aitken method. Example function call: [Control,Fevals]  = Linear_Aitken(1e-10,100,3);
- Linear_Steffensen.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear continuous control problem using the FBSM with the parital Steffensen method. Example function call: [Control,Fevals]  = Linear_Steffensen(1e-10,100,6); 
- Linear_Wegstein.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower' and 'upper' as described above, to solve the linear continuous control problem using the FBSM with Wegstein's method. Example function call: [Control,Fevals]  = Linear_Wegstein(1e-10,100,4,1,-2,0);
- Linear_Anderson.m - this function can be called with inputs 'tol', 'MaxFevals', 'M' and 'Droptol' as described above, to solve the linear continuous control problem using the FBSM with Anderson Acceleration. Example function call: [Control,Fevals]  = Linear_Anderson(1e-10,100,4,1e10);

The above scripts are dependent on the following (all contained within the Linear_continuous folder):
- Control.m - function for the control     
- State.m - function for the state equations
- Costate.m - function for the costate equations
- FBSM.m - subroutine that performs one iteration of the forward-backward sweep method and updates the control

**Linear_bangbang** - Corresponds to the linear bangbang control problem described in Section 3 of Sharp et al. 2021. 

This subfolder contains 5 scripts/functions for implementation:
- Linear_base_BB.m - this script can be run directly to solve the linear bangbang control problem using the standard FBSM without acceleration
- Linear_Aitken_BB.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear bangbang control problem using the FBSM with the parital Aitken method. Example function call: [Control,Fevals]  = Linear_Aitken_BB(1e-10,100,1);
- Linear_Steffensen_BB.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear bangbang control problem using the FBSM with the parital Steffensen method. Example function call: [Control,Fevals]  = Linear_Steffensen_BB(1e-10,100,7);
- Linear_Wegstein_BB.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower' and 'upper' as described above, to solve the linear bangbang control problem using the FBSM with Wegstein's method. Example function call: [Control,Fevals]  = Linear_Wegstein_BB(1e-10,100,1,0,0,0);
- Linear_Anderson_BB.m - this function can be called with inputs 'tol', 'MaxFevals', 'M' and 'Droptol' as described above, to solve the linear bangbang control problem using the FBSM with Anderson Acceleration. Example function call: [Control,Fevals]  = Linear_Anderson_BB(1e-10,100,1,1e10);

The above scripts are dependent on the following (all contained within the Linear_bangbang folder):
- Control_BB.m - function for the control     
- State.m - function for the state equations
- Costate_BB.m - function for the costate equations
- FBSM_BB.m - subroutine that performs one iteration of the forward-backward sweep method and updates the control

**Linear_fixedfinalstate** - Corresponds to the linear continuous control problem with fixed endpoint described in Section 3 of Sharp et al. 2021.

This subfolder contains 5 scripts/functions for implementation:
- Linear_FixedStateBothEnds.m - this script can be run directly to solve the linear fixed endpoint control problem using the adapted FBSM without acceleration
- Linear_FixedStateBothEnds_Aitken.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear fixed endpoint control problem using the adapted FBSM with the parital Aitken method. Example function call: [Control,Fevals]  = Linear_FixedStateBothEnds_Aitken(1e-10,100,2);
- Linear_FixedStateBothEnds_Steffensen.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear fixed endpoint control problem using the adapted FBSM with the parital Steffensen method. Example function call: [Control,Fevals]  = Linear_FixedStateBothEnds_Steffensen(1e-10,100,7);
- Linear_FixedStateBothEnds_Wegstein.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower' and 'upper' as described above, to solve the linear fixed endpoint control problem using the adapted FBSM with Wegstein's method. Example function call: [Control,Fevals]  = Linear_FixedStateBothEnds_Steffensen(1e-10,100,7);
- Linear_FixedStateBothEnds_Anderson.m - this function can be called with inputs 'tol', 'MaxFevals', 'M' and 'Droptol' as described above, to solve the linear fixed endpoint control problem using the adapted FBSM with Anderson Acceleration. Example function call: [Control,Fevals]  = Linear_FixedStateBothEnds_Anderson(1e-10,100,4,1e10); 

The above scripts are dependent on the following (all contained within the Linear_fixedfinalstate folder):
- Control.m - function for the control     
- State.m - function for the state equations
- Costate.m - function for the costate equations
- FBSM.m - subroutine that performs one iteration of the forward-backward sweep method and updates the control
- Sweeps.m - subroutine that solves the internal FBSM problems without acceleration for a guess, 'theta', of the adapted FBSM 
- Sweeps_Aitken.m - subroutine that solves the internal FBSM problems with the Aitken method for a guess, 'theta', of the adapted FBSM 
- Sweeps_Anderson.m - subroutine that solves the internal FBSM problems with Anderson Acceleration for a guess, 'theta', of the adapted FBSM 
- Sweeps_Steffensen.m - subroutine that solves the internal FBSM problems with the Steffensen method for a guess, 'theta', of the adapted FBSM 
- Sweeps_Wegstein.m - subroutine that solves the internal FBSM problems with the Wegstein method for a guess, 'theta', of the adapted FBSM 

**AML_continuous** - Corresponds to the AML continuous control problem described in Section 3 of Sharp et al. 2021. 

This subfolder contains 5 scripts/functions for implementation:
- AML_base.m - this script can be run directly to solve the AML continuous control problem using the standard FBSM without acceleration
- AML_Aitken.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the AML continuous control problem using the FBSM with the parital Aitken method. Example function call: [Control,Fevals]  = AML_Aitken(1e-10,100,5,0.5);  
- AML_Steffensen.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the AML continuous control problem using the FBSM with the parital Steffensen method. Example function call: [Control,Fevals]  = AML_Steffensen(1e-10,100,5,0.5);  
- AML_Wegstein.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower' and 'upper' as described above, to solve the AML continuous control problem using the FBSM with Wegstein's method. Example function call: [Control,Fevals]  = AML_Wegstein(1e-10,100,6,1,-1,1,0.55);  
- AML_Anderson.m - this function can be called with inputs 'tol', 'MaxFevals', 'M' and 'Droptol' as described above, to solve the AML continuous control problem using the FBSM with Anderson Acceleration. Example function call: [Control,Fevals]  = AML_Anderson(1e-10,100,6,1e10,0.85);  

The above scripts are dependent on the following (all contained within the AML_continuous folder):
- Control.m - function for the control     
- State.m - function for the state equations
- Costate.m - function for the costate equations
- FBSM.m - subroutine that performs one iteration of the forward-backward sweep method and updates the control

**AML_bangbang** - Corresponds to the AML bangbang control problem described in Section 3 of Sharp et al. 2021. 

This subfolder contains 5 scripts/functions for implementation:
- AML_base_BB.m - this script can be run directly to solve the AML bangbang control problem using the standard FBSM without acceleration
- AML_Aitken_BB.m - this function can be called with inputs 'tol', 'MaxFevals', 'm' and 'omega' as described above, to solve the AML bangbang control problem using the FBSM with the parital Aitken method. Example function call: [Control,Fevals]  = AML_Aitken_BB(1e-10,100,1,0.5);  
- AML_Steffensen_BB.m - this function can be called with inputs 'tol', 'MaxFevals', 'm' and 'omega' as described above, to solve the AML bangbang control problem using the FBSM with the parital Steffensen method. Example function call: [Control,Fevals]  = AML_Steffensen_BB(1e-10,100,5,0.5);
- AML_Wegstein_BB.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower', 'upper' and 'omega' as described above, to solve the AML bangbang control problem using the FBSM with Wegstein's method. Example function call: [Control,Fevals]  = AML_Wegstein_BB(1e-10,100,7,1,-1,1,0);
- AML_Anderson_BB.m - this function can be called with inputs 'tol', 'MaxFevals', 'M', 'Droptol' and 'omega' as described above, to solve the AML bangbang control problem using the FBSM with Anderson Acceleration. Example function call: [Control,Fevals]  = AML_Anderson_BB(1e-10,100,7,1e10,0.35);

The above scripts are dependent on the following (all contained within the AML_bangbang folder):
- Control_BB.m - function for the control     
- State.m - function for the state equations
- Costate_BB.m - function for the costate equations
- FBSM_BB.m - subroutine that performs one iteration of the forward-backward sweep method and updates the control

**AML_fixedfinalstate** - Corresponds to the AML continuous control problem with fixed endpoint described in Section 3 of Sharp et al. 2021.

This subfolder contains 5 scripts/functions for implementation:
- AML_FixedStateBothEnds.m - this script can be run directly to solve the AML fixed endpoint control problem using the adapted FBSM without acceleration
- AML_FixedStateBothEnds_Aitken.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the AML fixed endpoint control problem using the adapted FBSM with the parital Aitken method. Example function call: [Sol,Fevals] = AML_FixedStateBothEnds_Aitken(1e-10,1000,5);
- AML_FixedStateBothEnds_Steffensen.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the AML fixed endpoint control problem using the adapted FBSM with the parital Steffensen method. Example function call: [Sol,Fevals] = AML_FixedStateBothEnds_Steffensen(1e-10,1000,5);
- AML_FixedStateBothEnds_Wegstein.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower' and 'upper' as described above, to solve the AML fixed endpoint control problem using the adapted FBSM with Wegstein's method. Example function call: [Sol,Fevals] = AML_FixedStateBothEnds_Wegstein(1e-10,1000,6,1,-1,1); 
- AML_FixedStateBothEnds_Anderson.m - this function can be called with inputs 'tol', 'MaxFevals', 'M' and 'Droptol' as described above, to solve the AML fixed endpoint control problem using the adapted FBSM with Anderson Acceleration. Example function call: [Sol,Fevals] = AML_FixedStateBothEnds_Anderson(1e-10,1000,6,1e10);

The above scripts are dependent on the following (all contained within the AML_fixedfinalstate folder):
- Control.m - function for the control     
- State.m - function for the state equations
- Costate.m - function for the costate equations
- FBSM.m - subroutine that performs one iteration of the forward-backward sweep method and updates the control
- Sweeps.m - subroutine that solves the internal FBSM problems without acceleration for a guess, 'theta', of the adapted FBSM 
- Sweeps_Aitken.m - subroutine that solves the internal FBSM problems with the Aitken method for a guess, 'theta', of the adapted FBSM 
- Sweeps_Anderson.m - subroutine that solves the internal FBSM problems with Anderson Acceleration for a guess, 'theta', of the adapted FBSM 
- Sweeps_Steffensen.m - subroutine that solves the internal FBSM problems with the Steffensen method for a guess, 'theta', of the adapted FBSM 
- Sweeps_Wegstein.m - subroutine that solves the internal FBSM problems with the Wegstein method for a guess, 'theta', of the adapted FBSM 

We note that implementation of Anderson Acceleration involves computing condition numbers of matrices. Matrices that are ill-conditioned (indicated by a large condition number), are close to singular, such that significant numerical error can arise when computing the inverse, or obtaining the solution of a corresponding linear system of equations (Higham DJ, 1995. Condition numbers and their conditiuon numbers). It is known that for ill-conditioned matrices, computation of the condition number can itself be highly sensitive (Highham 1995). As such, computing condition numbers is commonly impacted by underflow and overflow, or rounding errors (Demmel JW, Li X, 1994. Faster numerical algorithms via exception handling, Highham 1995). For this reason, users attempting to reproduce results of the Anderson Acceleration method using different software or hardware may find that in some instances convergence is achieved with a different $\mathcal{N}$ to what is indicated in the tables, depending on whether the estimated condition number at each iteration suggests that the matrix is ill-conditioned. For a thorough discussion of condition numbers and issues arising from floating point arithmetic we direct readers to (Higham NJ, 2002. Accuracy and Stability of Numerical Algorithms).   
