#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
#include <hermes3d.h>
#include <iostream>

//  This example shows how to solve a time-dependent PDE discretized
//  in time via the implicit Euler method on a fixed mesh. 
//  You will see how to use the solution from previous time step.
//
//  PDE: stationary heat transfer equation
//  dT/dt - Laplace T = f.
//
//  Domain: (0, 0, 1)x(0, 1, 0)x(1, 0, 0), see the file hexahedron.mesh3d. 
//
//  Known exact solution, see unctions fn() and fndd(). 
//
//  IC:  T = 0.
//  BC:  T = 0. ... Dirichlet,
//
//  Time-stepping: implicit Euler.
//
//  The following parameters can be changed:

// Adaptivity threshold.
const int THRESHOLD = 0.0;	
const int INIT_REF_NUM = 3;			                  // Number of initial uniform mesh refinements.
const int P_INIT_X = 2,
          P_INIT_Y = 2,
          P_INIT_Z = 2;                           // Initial polynomial degree of all mesh elements.
const double TAU = 0.05;			                    // Time step in seconds. 
bool solution_output = true;                      // Generate output files (if true).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h)                                                  

// Problem parameters. 
const double FINAL_TIME = 2 * M_PI;		            // Length of time interval in seconds. 

// Global time variable. 
double TIME = TAU;

// Exact solution. 
#include "exact_solution.cpp"

// Boundary condition types.
BCType bc_types(int marker) {
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values. 
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}

#include "forms.cpp"

int criterion(Element* e)
{
  std::cout << e->id << std::endl;
  if(e->id % 3 == 0)
    return 0;
else
  return -1;
}

int main(int argc, char **args) 
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the initial mesh. 
  Mesh mesh1, mesh2;
  H3DReader mesh_loader;
  mesh_loader.load("hexahedron.mesh3d", &mesh1);
  mesh2.copy(mesh1);

  // Perform initial mesh refinement. 
  for (int i=0; i < INIT_REF_NUM; i++)
    mesh1.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
  mesh2.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

  // Some wild refinements to test multimesh.
  mesh2.refine_element(2, H3D_H3D_H3D_REFT_HEX_XYZ);

  mesh2.refine_element(11, H3D_H3D_H3D_REFT_HEX_XYZ);
  mesh2.refine_element(13, H3D_H3D_H3D_REFT_HEX_XYZ);

  mesh2.refine_element(14, H3D_REFT_HEX_X);
  mesh2.refine_element(15, H3D_H3D_H3D_REFT_HEX_XYZ);
  mesh2.refine_element(16, H3D_H3D_REFT_HEX_XZ);
  mesh2.refine_element(31, H3D_H3D_H3D_REFT_HEX_XYZ);

  mesh2.refine_element(43, H3D_H3D_H3D_REFT_HEX_XYZ);

  // Create H1 space with default shapeset.
  H1Space space1(&mesh1, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
  H1Space space2(&mesh2, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));

  // Construct initial solution and set it to zero.
  Solution sln_prev1(&mesh1);
  Solution sln_prev2(&mesh2);
  sln_prev1.set_zero();
  sln_prev2.set_zero();

  // Initialize weak formulation. 
  WeakForm wf(2);
  wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
  wf.add_matrix_form(1, 1, bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
  wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY, &sln_prev1);
  wf.add_vector_form(1, linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY, &sln_prev2);

  // Initialize discrete problem.
  bool is_linear = true;
  Hermes::Tuple<Space *> spaces(&space1, &space2);
  
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the preconditioner in the case of SOLVER_AZTECOO.
  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver)->set_solver(iterative_method);
    ((AztecOOSolver*) solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }

  // Time stepping. 
  int nsteps = (int) (FINAL_TIME/TAU + 0.5);
  for (int ts = 0; ts < nsteps;  ts++)
  {
    info("---- Time step %d, time %3.5f.", ts, TIME);
   
    // Adaptivity loop.
    int as = 1; 
    bool done = false;
    do 
    {
      info("------ Adaptivity step %d.", as);
		  // Assemble the linear problem.
		  Hermes::Tuple<Space *>* ref_spaces = construct_refined_spaces(spaces,0 , H3D_H3D_H3D_REFT_HEX_XYZ);

			DiscreteProblem dp(&wf, *ref_spaces, is_linear);

			Solution ref_sln1((*ref_spaces)[0]->get_mesh());
			Solution ref_sln2((*ref_spaces)[1]->get_mesh());
			Hermes::Tuple<Solution *> ref_slns(&ref_sln1, &ref_sln2);

		  info("Assembling the linear problem (ndof: %d).", Space::get_num_dofs(*ref_spaces));

		  bool rhsonly = (ts > 0);
		  dp.assemble(matrix, rhs, rhsonly);

		  // Solve the linear system. If successful, obtain the solution.
		  info("Solving the linear problem.");
		  if(solver->solve())
		  {
		    Solution::vector_to_solution(solver->get_solution(), (*ref_spaces)[0], &ref_sln1);
		    Solution::vector_to_solution(solver->get_solution(), (*ref_spaces)[1], &ref_sln2);
		  }
		  else
		    error ("Matrix solver failed.\n");

		  // Output solution.
		  if (solution_output)
		  {
		    out_fn_vtk(&ref_sln1, "sln1", ts);
		    out_fn_vtk(&ref_sln2, "sln2", ts);
		  }

		  // Output mesh with polynomial orders.
		  if (solution_output)
		  {
		    out_orders_vtk((*ref_spaces)[0], "order1", ts);
		    out_orders_vtk((*ref_spaces)[1], "order2", ts);
		  }

		  // Project the reference solutions on the coarse mesh.
		  Solution sln1(space1.get_mesh());
		  Solution sln2(space2.get_mesh());
		  Hermes::Tuple<Solution *> slns(&sln1, &sln2);
		  info("Projecting reference solution on coarse mesh.");
		  OGProjection::project_global(&space1, &ref_sln1, &sln1, matrix_solver);
		  OGProjection::project_global(&space2, &ref_sln2, &sln2, matrix_solver);

		  // Calculate error wrt. exact solution. 
		  info("Calculating exact error and error estimate.");
		  ExactSolution esln1((*ref_spaces)[0]->get_mesh(), fndd);
		  ExactSolution esln2((*ref_spaces)[1]->get_mesh(), fndd);
		  double err_exact = (h1_error(&sln1, &esln1) + h1_error(&sln2, &esln2)) * 100; 
		  info("Err. exact: %g%%.", err_exact);

		  Adapt *adaptivity = new Adapt(spaces, Hermes::Tuple<ProjNormType>(HERMES_H1_NORM, HERMES_H1_NORM));
		  bool solutions_for_adapt = true;
		  double err_estimate = adaptivity->calc_err_est(slns, ref_slns, solutions_for_adapt, HERMES_TOTAL_ERROR_ABS);
		  info("Err. estimate: %g%%.", err_estimate);

		  if (err_exact < THRESHOLD)
      {
		    done = true;

				// Next time step.
				sln_prev1 = ref_sln1;
				sln_prev2 = ref_sln2;
      }
		  else 
		  {
		    info("Adapting coarse mesh.");
		    adaptivity->adapt(THRESHOLD);
		  }

		  delete (*ref_spaces)[0]->get_mesh();
		  delete (*ref_spaces)[1]->get_mesh();
		  delete (*ref_spaces)[0];
		  delete (*ref_spaces)[1];
		  delete matrix;
		  delete rhs;
		  delete solver;
		  delete adaptivity;

		  // Increase the counter of performed adaptivity steps.
		  as++;
    } 
    while (!done);
    TIME += TAU;
  }

  // Time measurement.
  cpu_time.tick();

  // Print timing information.
  info("Solutions and mesh with polynomial orders saved. Total running time: %g s", cpu_time.accumulated());

  // Clean up.
  delete matrix;
  delete rhs;
  delete solver;

  return 0;
}
