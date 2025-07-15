// The code is based on the model of Wheeler et al.:
// Wheeler, A. A., Murray, B. T. & Schaefer, R. J. Computation of dendrites using a phase field model. Physica D 66, 243-262, (1993).


// =================================================================================
// Define the variables in the model
// =================================================================================
// The number of variables
#define num_var 2

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define variable_name {"c", "n"}
#define variable_type {"SCALAR","SCALAR"}
#define variable_eq_type {"PARABOLIC","PARABOLIC"}

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

//define Cahn-Hilliard parameters (No Gradient energy term)
#define alpha 400
#define epsilon 0.0025
#define delta 0.5
#define m 0.05

//anisotropy and regularization parameters
#define epsilonM 0.05
#define coeff (epsilon*alpha*delta)

//anisotropy gamma as a function of the components of the normal vector
//current anisotropy has 4-fold or octahedral symmetry
#if problemDIM==1
#define gamma 1.0
#elif problemDIM==2
//writing out powers instead of using std::pow(double,double) for performance reasons
#define gamma (1.0+epsilonM*(4.0*(normal[0]*normal[0]*normal[0]*normal[0]+normal[1]*normal[1]*normal[1]*normal[1])-3.0))
#else
#define gamma (1.0+epsilonM*(4.0*(normal[0]*normal[0]*normal[0]*normal[0]+normal[1]*normal[1]*normal[1]*normal[1]+normal[2]*normal[2]*normal[2]*normal[2])-3.0))
#endif

//derivatives of gamma with respect to the components of the unit normal
#define gammanx (epsilonM*16.0*normal[0]*normal[0]*normal[0])
#define gammany (epsilonM*16.0*normal[1]*normal[1]*normal[1])
#define gammanz (epsilonM*16.0*normal[2]*normal[2]*normal[2])

#if problemDIM==2
#define gammax (gammanx*normalx[0][0]+gammany*normalx[0][1])
#define gammay (gammanx*normalx[1][0]+gammany*normalx[1][1])
#define gammaz (gammanx*normalx[2][0]+gammany*normalx[2][1])

#elif problemDIM==3
#define gammax (gammanx*normalx[0][0]+gammany*normalx[0][1]+gammanz*normalx[0][2])
#define gammay (gammanx*normalx[1][0]+gammany*normalx[1][1]+gammanz*normalx[1][2])
#define gammaz (gammanx*normalx[2][0]+gammany*normalx[2][1]+gammanz*normalx[2][2])
#endif

#define gammanxx (epsilonM*48.0*normal[0]*normal[0]*normalx[0][0])
#define gammanyy (epsilonM*48.0*normal[1]*normal[1]*normalx[1][1])
#define gammanzz (epsilonM*48.0*normal[2]*normal[2]*normalx[2][2])

//Allen-Cahn mobility (isotropic)
#define tau1 (epsilon*epsilon/m)
#define tau2 (epsilon*epsilon)

//Allen-Cahn mobility (anisotropic)
//#define MnV (1.0/(gamma*gamma+1e-10))
#define An (n*(1.0-n)*(n-0.5+30.0*coeff*c*n*(1.0-n)))
#define Bn ((n-oldn)*30.0*n*n*(1.0-n)*(1.0-n))
//define required residuals (aniso defined in model)

#define rcV   (c-constV(1.0/delta)*Bn)
#define rcxV  (constV(-timeStep)*cx)
#define rnV  (n+constV(timeStep/tau1)*(constV(tau2)*aniso+An))
#define rnxV (constV(tau2*timeStep/tau1)*(-gamma*gamma*nx))


// =================================================================================
// residualRHS
// =================================================================================
// This function calculates the residual equations for each variable. It takes
// "modelVariablesList" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of
// that quadrature point is given by "q_point_loc". The function outputs
// "modelResidualsList", a list of the value and gradient terms of the residual for
// each residual equation. The index for each variable in these lists corresponds to
// the order it is defined at the top of this file (starting at 0).template <int dim>
template <int dim>
void generalizedProblem<dim>::residualRHS(const std::vector<std::vector<modelVariable<dim>>*> & modelVariablesList, 
     std::vector<modelResidual<dim>> & modelResidualsList, 
					  dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

//c
 scalarvalueType c = (*modelVariablesList[0])[0].scalarValue;

 scalargradType cx = (*modelVariablesList[0])[0].scalarGrad;
 
//n
 scalarvalueType n = (*modelVariablesList[0])[1].scalarValue;
 scalarvalueType oldn = (*modelVariablesList[1])[1].scalarValue;
 scalargradType nx = (*modelVariablesList[0])[1].scalarGrad;
 scalarhessType nxx = (*modelVariablesList[0])[1].scalarHess;

// anisotropy code
 scalarvalueType normgradn = std::sqrt(nx.norm_square());
 scalargradType normal = nx/(normgradn+constV(1.0e-16));
 scalarhessType normalx = nxx/(normgradn+constV(1.0e-16));
 scalarvalueType gamma_scl = gamma;
 scalarvalueType aniso = constV(0.0);

#if problemDIM>1
 scalargradType dgammadnorm;
 dgammadnorm[0]=gammanx;
 dgammadnorm[1]=gammany;
 
 scalargradType dgammadnormdx;
 dgammadnormdx[0] = gammanxx;
 dgammadnormdx[1] = gammanyy;

 scalargradType gammax_scl;
 gammax_scl[0] = gammax;
 gammax_scl[1] = gammay;

#if problemDIM>2
 dgammadnorm[2]=gammanz;
 dgammadnormdx[2] = gammanzz;
 gammax_scl[2] = gammaz;
#endif
      for (unsigned int i=0; i<problemDIM; ++i){
	aniso += normgradn*(gammax_scl[i]*dgammadnorm[i]
			    +gamma_scl*dgammadnormdx[i]);
	/*	
	if (n[0]>0.1 && n[0]<0.9)
	  std::cout << "aniso: " << aniso[0] << std::endl
		    << "normxx: "<< normalx[0][0][0] << std::endl
		    << "normyy: "<< normalx[1][1][0] << std::endl
		    << "normzz: "<< normalx[2][2][0] << std::endl
		    << "normgrad: "<< normgradn[0] << std::endl
		    << "gammax: " << gammax_scl[i][0] << std::endl
		    << "dgammadnorm: " << dgammadnorm[i][0] << std::endl
		    << "gamma: " << gamma_scl[0] << std::endl
		    << "dgammadnormdx: " << dgammadnormdx[i][0] <<std::endl;
	*/
      }
#endif
// end anisotropy code

modelResidualsList[0].scalarValueResidual = rcV;
modelResidualsList[0].scalarGradResidual = rcxV;

modelResidualsList[1].scalarValueResidual = rnV;
modelResidualsList[1].scalarGradResidual = rnxV;
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
