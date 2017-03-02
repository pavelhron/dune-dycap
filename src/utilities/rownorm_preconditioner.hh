// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_ROWNORM_PRECONDITIONER_HH
#define DUNE_DYCAP_ROWNORM_PRECONDITIONER_HH

#include <dune/common/float_cmp.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/backend/istl.hh>


//! class to change the matrix and the right side
template <class GV, class GFS>
class RowNormPreconditioner
{

public:
  // Grid
  typedef typename GV::Grid Grid;
  typedef typename GV::IndexSet IndexSet;
  typedef GV GridView;
  typedef typename Grid::ctype ctype;
  enum
    { dim = Grid::dimension /* !<dimension of the grid we are using */  };

  // Local Matrix Block
  static const int BlockSize = GFS::Traits::BackendType::BlockSize;

  typedef Dune::FieldVector < ctype, BlockSize > LocalVectorBlock;
  typedef Dune::FieldMatrix < ctype, BlockSize, BlockSize > LocalMatrixBlock;

  typedef Dune::BCRSMatrix < LocalMatrixBlock > Matrix;
  typedef Dune::BlockVector < LocalVectorBlock > SolutionVector;

private:


  typedef typename GridView::template Codim<0>::Iterator Iterator;
  typedef typename GridView::IntersectionIterator IntersectionIterator;

  /*
    This functions are used to scale the both rows in the jacobian matrix. The scale factor is
    the infinity norm in the every row => if the biggest value is on the diagonal,
    the diagonal in jacobian matrix consists of identities
  */
  void getBlockNorm(LocalMatrixBlock &m, LocalVectorBlock & norm)
  {
    typename LocalMatrixBlock::RowIterator fit  = m.begin();
    typename LocalMatrixBlock::RowIterator efit = m.end();

    for(; fit!=efit; ++fit)
      norm[fit.index()] = fit->infinity_norm();
  }

  template <bool diagonal>
  void rowNormalizeBlock(LocalMatrixBlock & m, const LocalVectorBlock & norm)
  {

    typename LocalMatrixBlock::Iterator fit  = m.begin();
    typename LocalMatrixBlock::Iterator efit = m.end();

    for(; fit!=efit; ++fit){
      const ctype n = norm[fit.index()];
      if(n)
        *fit /= n;
      else if(diagonal)
        (*fit)[fit.index()] = 1.0;
    }
  }

  void rowNormalizeRhs(LocalVectorBlock & m, const LocalVectorBlock & norm)
  {

    typename LocalVectorBlock::Iterator fit  = m.begin();
    typename LocalVectorBlock::Iterator efit = m.end();

    for(; fit!=efit; ++fit){
      const ctype n = norm[fit.index()];
      if(n)
        *fit /= n;

    }
  }

public:

  //! Perform equilibration of matrix rows
  RowNormPreconditioner(Matrix &A, SolutionVector &B)
  {

    typedef typename Matrix::RowIterator RowIterator;
    typedef typename Matrix::ColIterator ColIterator;
    typedef typename Matrix::size_type size_type;

    {
      RowIterator rit = A.begin();
      const RowIterator erit = A.end();
      for(; rit != erit; ++rit){
        const size_type r_index = rit.index();

        LocalMatrixBlock & diag = (*rit)[r_index];
        LocalVectorBlock vnorm;

        // Get equilibration coefficients (= row norms)
        getBlockNorm( diag, vnorm);

        // Equilibrate row and rhs
        rowNormalizeBlock<true>(diag, vnorm);
        rowNormalizeRhs(B[r_index], vnorm);

        // Iterate off diagonal elements
        ColIterator cit = rit->begin();
        const ColIterator ecit = rit->end();

        for(; cit != ecit; ++cit){

          const size_type c_index = cit.index();

          // Skipt diagonal
          if(c_index == r_index)
            continue;

          rowNormalizeBlock<false>(*cit, vnorm);

        }// cit
      }// rit
    }

  }


};

#endif
