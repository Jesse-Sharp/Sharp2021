# Sharp2021
Code for implementing iterative acceleration techniques for general systems and within the forward-backward sweep method (FBSM). 

This repository contains 5 MATLAB scripts and 6 subfolders, corresponding to each of the linear and acute myeloid leukaemia (AML) control problems considered in the work of Sharp et al. 2021. Where functions can be called by users, we provide example function calls corresponding to results presented in the paper. 

MATLAB scripts in the top level of the repository are provided for implementing the fixed point iteration method and acceleration algorithms considered in the work, for a system of arbitrary size N: 

Example F3 = @(X) [X(1)-(X(1)^2+X(2)^2-5)/4;
 X(2)-(X(1)*X(2)-2)/2;
-(X(1)*X(3)-X(2))/3];

FixedPoint.m
Aitken.m  
Anderson.m 
Steffensen.m
Wegstein.m 

These scripts solve equations of the form _X = F(X)_, for the user specified 'F', initial guess 'X0', convergence tolerance 'tol', and maximum function evaluations 'MaxFevals'. The acceleration methods each require user specified tuning parameters. For the Aitken and Steffensen methods, input 'm' specifies the desired dimension of the _N x m_ difference matrices for the "partial" implementation discussed in Sharp2021; setting _m = N_ corresponds to the standard implementation of the Aitken and Steffensen methods. For Anderson Acceleration, input 'M' determines the maximum number of previous iterations to incorporate in each iteration, while input 'DropTol' specifies the maximum accptable condition number of the residual difference matrix, _dG_. For Wegstein's method, input 'nth' specifies how frequently to update q; every nth iteration, input 'bounding' acts as a switch to turn on or off bounds on _q_, with bounds applied if _bounding=1_. Inputs 'lower' and 'upper' specify the lower and upper bounds to apply when _bounding = 1_.



Each subfolder corresponds to a control problem presented in Sharp et al. 2021. Their contents and description are provided: 

**Linear_continuous** - Corresponds to the linear continuous control problem described in Section 3 of Sharp et al. 2021. 

This subfolder contains 5 scripts/functions for implementation:
- Linear_base.m - this script can be run directly to solve the linear continuous control problem using the standard FBSM without acceleration
- Linear_Aitken.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear continuous control problem using the FBSM with the parital Aitken method. Example: [] = 
- Linear_Steffensen.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear continuous control problem using the FBSM with the parital Steffensen method
- Linear_Wegstein.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower' and 'upper' as described above, to solve the linear continuous control problem using the FBSM with Wegstein's method
- Linear_Anderson.m - this function can be called with inputs 'tol', 'MaxFevals', 'M' and 'DropTol' as described above, to solve the linear continuous control problem using the FBSM with Anderson Acceleration

The above scripts are dependent on the following (all contained within the Linear_continuous folder):
- Control.m - function for the control     
- State.m - function for the state equations
- Costate.m - function for the costate equations
- FBSM.m - subroutine that performs one iteration of the forward-backward sweep method and updates the control

**Linear_bangbang** - Corresponds to the linear bangbang control problem described in Section 3 of Sharp et al. 2021. 

This subfolder contains 5 scripts/functions for implementation:
- Linear_base_BB.m - this script can be run directly to solve the linear bangbang control problem using the standard FBSM without acceleration
- Linear_Aitken_BB.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear bangbang control problem using the FBSM with the parital Aitken method
- Linear_Steffensen_BB.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear bangbang control problem using the FBSM with the parital Steffensen method
- Linear_Wegstein_BB.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower' and 'upper' as described above, to solve the linear bangbang control problem using the FBSM with Wegstein's method
- Linear_Anderson_BB.m - this function can be called with inputs 'tol', 'MaxFevals', 'M' and 'DropTol' as described above, to solve the linear bangbang control problem using the FBSM with Anderson Acceleration

The above scripts are dependent on the following (all contained within the Linear_bangbang folder):
- Control_BB.m - function for the control     
- State.m - function for the state equations
- Costate_BB.m - function for the costate equations
- FBSM_BB.m - subroutine that performs one iteration of the forward-backward sweep method and updates the control

**Linear_fixedfinalstate** - Corresponds to the linear continuous control problem with fixed endpoint described in Section 3 of Sharp et al. 2021.

This subfolder contains 5 scripts/functions for implementation:
- Linear_FixedStateBothEnds.m - this script can be run directly to solve the linear fixed endpoint control problem using the adapted FBSM without acceleration
- Linear_FixedStateBothEnds_Aitken.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear fixed endpoint control problem using the adapted FBSM with the parital Aitken method
- Linear_FixedStateBothEnds_Steffensen.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the linear fixed endpoint control problem using the adapted FBSM with the parital Steffensen method
- Linear_FixedStateBothEnds_Wegstein.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower' and 'upper' as described above, to solve the linear fixed endpoint control problem using the adapted FBSM with Wegstein's method
- Linear_FixedStateBothEnds_Anderson.m - this function can be called with inputs 'tol', 'MaxFevals', 'M' and 'DropTol' as described above, to solve the linear fixed endpoint control problem using the adapted FBSM with Anderson Acceleration

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
- AML_Aitken.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the AML continuous control problem using the FBSM with the parital Aitken method
- AML_Steffensen.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the AML continuous control problem using the FBSM with the parital Steffensen method
- AML_Wegstein.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower' and 'upper' as described above, to solve the AML continuous control problem using the FBSM with Wegstein's method
- AML_Anderson.m - this function can be called with inputs 'tol', 'MaxFevals', 'M' and 'DropTol' as described above, to solve the AML continuous control problem using the FBSM with Anderson Acceleration

The above scripts are dependent on the following (all contained within the AML_continuous folder):
- Control.m - function for the control     
- State.m - function for the state equations
- Costate.m - function for the costate equations
- FBSM.m - subroutine that performs one iteration of the forward-backward sweep method and updates the control

**AML_bangbang** - Corresponds to the AML bangbang control problem described in Section 3 of Sharp et al. 2021. 

This subfolder contains 5 scripts/functions for implementation:
- AML_base_BB.m - this script can be run directly to solve the AML bangbang control problem using the standard FBSM without acceleration
- AML_Aitken_BB.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the AML bangbang control problem using the FBSM with the parital Aitken method
- AML_Steffensen_BB.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the AML bangbang control problem using the FBSM with the parital Steffensen method
- AML_Wegstein_BB.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower' and 'upper' as described above, to solve the AML bangbang control problem using the FBSM with Wegstein's method
- AML_Anderson_BB.m - this function can be called with inputs 'tol', 'MaxFevals', 'M' and 'DropTol' as described above, to solve the AML bangbang control problem using the FBSM with Anderson Acceleration

The above scripts are dependent on the following (all contained within the AML_bangbang folder):
- Control_BB.m - function for the control     
- State.m - function for the state equations
- Costate_BB.m - function for the costate equations
- FBSM_BB.m - subroutine that performs one iteration of the forward-backward sweep method and updates the control

**AML_fixedfinalstate** - Corresponds to the AML continuous control problem with fixed endpoint described in Section 3 of Sharp et al. 2021.

This subfolder contains 5 scripts/functions for implementation:
- AML_FixedStateBothEnds.m - this script can be run directly to solve the AML fixed endpoint control problem using the adapted FBSM without acceleration
- AML_FixedStateBothEnds_Aitken.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the AML fixed endpoint control problem using the adapted FBSM with the parital Aitken method
- AML_FixedStateBothEnds_Steffensen.m - this function can be called with inputs 'tol', 'MaxFevals' and 'm' as described above, to solve the AML fixed endpoint control problem using the adapted FBSM with the parital Steffensen method
- AML_FixedStateBothEnds_Wegstein.m - this function can be called with inputs 'tol', 'MaxFevals', 'nth', 'bounding', 'lower' and 'upper' as described above, to solve the AML fixed endpoint control problem using the adapted FBSM with Wegstein's method
- AML_FixedStateBothEnds_Anderson.m - this function can be called with inputs 'tol', 'MaxFevals', 'M' and 'DropTol' as described above, to solve the AML fixed endpoint control problem using the adapted FBSM with Anderson Acceleration

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
