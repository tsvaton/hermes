// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_PRECALC_H
#define __H2D_PRECALC_H

#include "../function/function.h"
#include "../shapeset/shapeset.h"

namespace Hermes
{
  namespace Hermes2D
  {
    enum SpaceType;
    /// @ingroup meshFunctions
    /// \brief Caches precalculated shape function values.
    ///
    /// PrecalcShapeset is a cache of precalculated shape function values.
    ///
    ///
    class HERMES_API PrecalcShapeset : public Function<double>
    {
    public:
      /// Returns type of space
      SpaceType get_space_type() const;

      /// \brief Constructs a standard (master) precalculated shapeset class.
      /// \param shapeset[in] Pointer to the shapeset to be precalculated.
      PrecalcShapeset(Shapeset* shapeset);

      /// Destructor.
      virtual ~PrecalcShapeset();

      /// Activates a shape function given by its index. The values of the shape function
      /// can then be obtained by setting the required integration rule order by calling
      /// set_quad_order() and after that calling get_values(), get_dx_values(), etc.
      /// \param index[in] Shape index.
      void set_active_shape(int index);

    private:
      virtual void set_quad_2d(Quad2D* quad_2d);

      /// \brief Frees all precalculated tables.
      virtual void free();

      /// Returns the index of the active shape (can be negative if the shape is constrained).
      int get_active_shape() const;

      /// Returns a pointer to the shapeset which is being precalculated.
      Shapeset* get_shapeset() const;

      /// Returns the polynomial order of the active shape function on given edge.
      virtual int get_edge_fn_order(int edge);

      Shapeset* shapeset;

      int index;

      int max_index[H2D_NUM_MODES];

      /// Transformed points to the reference domain, used by precalculate.
      double2 ref_points[H2D_MAX_INTEGRATION_POINTS_COUNT];

      virtual void precalculate(int order, int mask);

      void update_max_index();

      friend class RefMap;
      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class Adapt;
      template<typename T> friend class Func;
      template<typename T> friend class Solution;
      template<typename T> friend class DiscontinuousFunc;
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemDGAssembler;
      template<typename T> friend class DiscreteProblemThreadAssembler;
      template<typename T> friend class NeighborSearch;
      friend class CurvMap;
    };
  }
}
#endif
