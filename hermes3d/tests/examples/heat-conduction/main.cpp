#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
#include <hermes3d.h>

// This test makes sure that the example heat-conduction works correctly.

const int INIT_REF_NUM = 1;                         // Number of initial uniform mesh refinements.
const int P_INIT_X = 2, P_INIT_Y = 2, P_INIT_Z = 2; // Initial polynomial degree of all mesh elements.
const double TAU = 0.05;                            // Time step in seconds. 
MatrixSolverType matrix_solver = SOLVER_UMFPACK;    // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.
const double ERR_STOP = 5.0;                        // Stopping criterion for adaptivity (rel. error tolerance between the
// fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                       // Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const double THRESHOLD = 0.5;		                    // error threshold for element refinement

unsigned int DEREFINEMENT_DONE_EVERY_N_STEPS = 1;		// derefinement performed every N steps.
bool SPACE_ADAPTED_SINCE_LAST_DEREF = false;		    // If there was no adaptation since last derefinement, we do not derefine.

// Problem parameters. 
const double FINAL_TIME = 2 * M_PI;		              // Length of time interval in seconds. 

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

int main(int argc, char **args) 
{
  // Test variable.
  int success_test = 1;

  // Load the initial mesh. 
  Mesh mesh;
  H3DReader mesh_loader;
  mesh_loader.load("hexahedron.mesh3d", &mesh);

  // Perform initial mesh refinement. 
  for (int i=0; i < INIT_REF_NUM; i++) 
    mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

  // Create H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));

  // Construct initial solution and set it to zero.
  Solution sln_prev(&mesh);
  sln_prev.set_zero();

  // Initialize weak formulation. 
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
  wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY, &sln_prev);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Error estimate and exact error.
  double err_est, err_exact;

  // Solutions
  Solution sln(space.get_mesh());

  // Time stepping. 
  int nsteps = (int) (FINAL_TIME/TAU + 0.5);
  for (int ts = 0; ts < nsteps;  ts++)
  {
    // Global derefinement.
    if(SPACE_ADAPTED_SINCE_LAST_DEREF && ts % DEREFINEMENT_DONE_EVERY_N_STEPS == 0)
    {
      info("Global space derefinement.");
      SPACE_ADAPTED_SINCE_LAST_DEREF = false;
      space.unrefine_all_mesh_elements();
    }

    info("---- Time step %d, time %3.5f.", ts, TIME);

    // Adaptivity loop.
    int as = 1; 
    bool done = false;
    do 
    {
      info("---- Adaptivity step %d:", as);
      Space* ref_space = construct_refined_space(&space,1 , H3D_H3D_H3D_REFT_HEX_XYZ);

      // Initialize discrete problem.
      bool is_linear = true;
      DiscreteProblem dp(&wf, ref_space, is_linear);

      // Assemble the linear problem.
      info("Assembling the linear problem (ndof: %d).", Space::get_num_dofs(ref_space));

      dp.assemble(matrix, rhs);
      Solution rsln(ref_space->get_mesh());

      // Solve the linear system. If successful, obtain the solution.
      info("Solving the linear problem.");
      if(solver->solve()) 
        Solution::vector_to_solution(solver->get_solution(), ref_space, &rsln);
      else 
        error ("Matrix solver failed.\n");

      info("Project onto the coarse space.");
      OGProjection::project_global(&space, &rsln, &sln);

      // Output solution.
        out_fn_vtk(&rsln, "sln", ts);

      // Exact solution.
      ExactSolution esln(ref_space->get_mesh(), fndd);
  
      info("Calculating error estimate and the exact error.");
      Adapt *adaptivity = new Adapt(&space, HERMES_H1_NORM);
      bool solutions_for_adapt = false;
      err_exact = adaptivity->calc_err_exact(&sln, &esln, solutions_for_adapt, HERMES_TOTAL_ERROR_REL) * 100;
      solutions_for_adapt = true;
      err_est = adaptivity->calc_err_est(&sln, &rsln, solutions_for_adapt, HERMES_TOTAL_ERROR_REL) * 100;
      info("Err. est: %g%%.", err_est);
      info("(Probably wrong) Err. exact: %g%%.", err_exact);

      // Output.
      out_mesh_vtk(ref_space->get_mesh(), "Mesh-fine", ts, as);
      out_orders_vtk(ref_space, "Space-fine", ts, as);
      out_fn_vtk(&rsln, "Solution-fine-vector", ts, as);

      // If err_est_rel is too large, adapt the mesh. 
      if (err_est < ERR_STOP) 
      {
        done = true;
      }	
      else 
      {
        info("Adapting coarse mesh.");
        adaptivity->adapt(THRESHOLD);
        SPACE_ADAPTED_SINCE_LAST_DEREF = true;
      }

      // If we reached the maximum allowed number of degrees of freedom, set the return flag to failure.
      if (Space::get_num_dofs(&space) >= NDOF_STOP)
        done = true;

      // Cleanup.
      delete ref_space;
      delete adaptivity;

      // Increase the counter of performed adaptivity steps.
      as++;
    } 
    while (!done);

    // Next time step.
    sln_prev = sln;
    TIME += TAU;
  }

  // Clean up.
  delete matrix;
  delete rhs;
  delete solver;

  return 0;
}
