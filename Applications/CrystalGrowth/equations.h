// The code is based on the model of Wheeler et al.:
// Wheeler, A. A., Murray, B. T. & Schaefer, R. J. Computation of dendrites using a phase field model. Physica D 66, 243-262, (1993).


// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 2

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"phi", "T"}
#define variable_type {"SCALAR", "SCALAR"}
#define variable_eq_type {"PARABOLIC", "PARABOLIC"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define need_val {true, true}
#define need_grad {true, true}
#define need_hess {false, true}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define need_val_residual {true, true}
#define need_grad_residual {true, true}

// =================================================================================
// Define the model parameters and the residual equations
// =================================================================================
// Parameters in the residual equations and expressions for the residual equations
// can be set here. For simple cases, the entire residual equation can be written
// here. For more complex cases with loops or conditional statements, residual
// equations (or parts of residual equations) can be written below in "residualRHS".

// Material parameters (ice-water)
#define T_M 273.15
#define L 3.34e5
#define Cp 4.2e3
#define D 1.3e-7
#define W 1e-7
#define tau 1e-4
#define lambda 10.0
#define epsilon_a 0.02

// Hexagonal anisotropy (6-fold)
#define a_theta (1.0 + epsilon_a * cos(6.0 * theta))

// Residuals
template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<std::vector<modelVariable<dim>>*> & modelVariablesList, 
     std::vector<modelResidual<dim>> & modelResidualsList, 
     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

  // phi
  scalarvalueType phi = (*modelVariablesList[0])[0].scalarValue;
  scalargradType phix = (*modelVariablesList[0])[0].scalarGrad;
  // T
  scalarvalueType T = (*modelVariablesList[0])[1].scalarValue;
  scalargradType Tx = (*modelVariablesList[0])[1].scalarGrad;

  // Compute theta from gradient (for anisotropy)
  double theta = atan2(phix[1], phix[0]); // 2D; for 3D, use spherical coordinates

  // Anisotropy
  double a = a_theta;
  double a2 = a * a;

  // Phase-field equation (Allen-Cahn/Karma)
  double u = (T - T_M) / 1.0; // scale denominator as needed
  double phase_rhs = (phi - phi*phi)*(phi - 0.5 + lambda*u);
  modelResidualsList[0].scalarValueResidual = (phi - old_phi)/tau - W*W*div(a2*grad(phi)) - phase_rhs;

  // Heat equation with latent heat
  double latent = (L/Cp) * (phi - old_phi)/tau;
  modelResidualsList[1].scalarValueResidual = (T - old_T)/tau - D*laplacian(T) - latent;
}

// =================================================================================
// residualLHS (needed only if at least one equation is elliptic)
// =================================================================================
// This function calculates the residual equations for the iterative solver for
// elliptic equations.for each variable. It takes "modelVariablesList" as an input,
// which is a list of the value and derivatives of each of the variables at a
// specific quadrature point. The (x,y,z) location of that quadrature point is given
// by "q_point_loc". The function outputs "modelRes", the value and gradient terms of
// for the left-hand-side of the residual equation for the iterative solver. The
// index for each variable in these lists corresponds to the order it is defined at
// the top of this file (starting at 0), not counting variables that have
// "need_val_LHS", "need_grad_LHS", and "need_hess_LHS" all set to "false". If there
// are multiple elliptic equations, conditional statements should be used to ensure
// that the correct residual is being submitted. The index of the field being solved
// can be accessed by "this->currentFieldIndex".
template <int dim>
void generalizedProblem<dim>::residualLHS(const std::vector<modelVariable<dim>> & modelVarList,
    modelResidual<dim> & modelRes, dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

}

// =================================================================================
// energyDensity (needed only if calcEnergy == true)
// =================================================================================
// This function integrates the free energy density across the computational domain.
// It takes "modelVariablesList" as an input, which is a list of the value and
// derivatives of each of the variables at a specific quadrature point. It also
// takes the mapped quadrature weight, "JxW_value", as an input. The (x,y,z) location
// of the quadrature point is given by "q_point_loc". The weighted value of the
// energy density is added to "energy" variable and the components of the energy
// density are added to the "energy_components" variable (index 0: chemical energy,
// index 1: gradient energy, index 2: elastic energy).
template <int dim>
void generalizedProblem<dim>::energyDensity(const std::vector<modelVariable<dim>> & modelVarList,
    const dealii::VectorizedArray<double> & JxW_value,
    dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {

}
