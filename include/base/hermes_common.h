// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file hermes_common.h
    \brief File containing includes of all HermesCommon functionality + solvers. Intended to be included.
*/
#include "util/common.h"
#include "algebra/solvers/linear_matrix_solver.h"
#include "algebra/solvers/nonlinear_solver.h"
#include "algebra/solvers/amesos_solver.h"
#include "algebra/solvers/aztecoo_solver.h"
#include "algebra/solvers/epetra.h"
#include "algebra/solvers/mumps_solver.h"
#include "algebra/solvers/newton_solver_nox.h"
#include "algebra/solvers/petsc_solver.h"
#include "algebra/solvers/umfpack_solver.h"
#include "algebra/solvers/superlu_solver.h"
#include "algebra/solvers/precond.h"
#include "algebra/solvers/precond_ifpack.h"
#include "algebra/solvers/precond_ml.h"
#include "algebra/solvers/eigensolver.h"
#include "base/hermes_function.h"
#include "util/compat.h"
#include "util/callstack.h"
#include "util/vector.h"
#include "algebra/tables.h"
#include "util/array.h"
#include "util/qsort.h"
#include "base/ord.h"
#include "util/mixins.h"
#include "util/api.h"