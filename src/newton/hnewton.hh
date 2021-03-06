// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=4 sts=4:
#ifdef DUNE_PDELAB_NEWTON_HH
#warning ***** WARNING ***** newton.hh was already included ******
#endif

#ifndef DUNE_PDELAB_NEWTON_HH
#define DUNE_PDELAB_NEWTON_HH

#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>

#include <dune/common/exceptions.hh>
#include <dune/common/ios_state.hh>
#include <dune/common/timer.hh>
#include <dune/pdelab/backend/solver.hh>

#include "newton_utilities.hh"

namespace Dune
{
    namespace PDELab
    {
        // Exception classes used in NewtonSolver
        class NewtonError : public Exception {};
        class NewtonDefectError : public NewtonError {};
        class NewtonLinearSolverError : public NewtonError {};
        class NewtonLineSearchError : public NewtonError {};
        class NewtonNotConverged : public NewtonError {};



        // Status information of Newton's method
        template<class RFType>
        struct NewtonResult : LinearSolverResult<RFType>
        {
            //            RFType conv_rate;          // average reduction per Newton iteration
            RFType first_defect;       // the first defect
            RFType defect;             // the final defect
            double assembler_time;     // Cumulative time for matrix assembly
            double linear_solver_time; // Cumulative time for linear sovler
            int linear_solver_iterations; // Total number of linear iterations

            NewtonResult() :
                first_defect(0.0), defect(0.0), assembler_time(0.0), linear_solver_time(0.0), linear_solver_iterations(0) {} // conv_rate(0.0),

        };

        template<class GOS, class TrlV, class TstV, class HELPER>
        class NewtonBase
        {
            typedef GOS GridOperator;
            typedef TrlV TrialVector;
            typedef TstV TestVector;

            typedef typename TestVector::ElementType RFType;
            typedef typename GOS::template MatrixContainer<RFType>::Type Matrix;

            typedef NewtonResult<RFType> Result;

        public:
            void setVerbosityLevel(unsigned int verbosity_level_)
            {
                if (gridoperator.trialGridFunctionSpace().gridView().comm().rank()>0)
                    verbosity_level = 0;
                else
                    verbosity_level = verbosity_level_;
            }

        protected:
            GridOperator& gridoperator;
            TrialVector *u;
            HELPER *helper;
            Result res;
            unsigned int verbosity_level;
            RFType prev_defect;
            RFType linear_reduction;
            bool reassembled;

            NewtonBase(GridOperator& go, TrialVector& u_, HELPER* helper_)
                : gridoperator(go)
                , u(&u_)
                , helper(helper_)
                , verbosity_level(1)
            {
                if (gridoperator.trialGridFunctionSpace().gridView().comm().rank()>0)
                    verbosity_level = 0;
            }

            NewtonBase(GridOperator& go, HELPER* helper_)
                : gridoperator(go)
                , u(0)
                , helper(helper_)
                , verbosity_level(1)
            {
                if (gridoperator.trialGridFunctionSpace().gridView().comm().rank()>0)
                    verbosity_level = 0;
            }

            virtual ~NewtonBase() { }

            virtual bool terminate() = 0;
            virtual void prepare_step(Matrix& A, TestVector& r) = 0;
            virtual void line_search(TrialVector& z, TestVector& r) = 0;
            virtual void defect(TestVector& r) = 0;
        };

        template<class GOS, class S, class TrlV, class TstV, class HELPER>
        class NewtonSolver : public virtual NewtonBase<GOS,TrlV,TstV,HELPER>
        {
            typedef S Solver;
            typedef GOS GridOperator;
            typedef TrlV TrialVector;
            typedef TstV TestVector;

            typedef typename TestVector::ElementType RFType;
            typedef typename GOS::template MatrixContainer<RFType>::Type Matrix;

        public:
            typedef NewtonResult<RFType> Result;

            NewtonSolver(GridOperator& go, TrialVector& u_, Solver& solver_, HELPER* helper_)
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go,u_,helper_)
                , solver(solver_)
                , result_valid(false)
            {}

            NewtonSolver(GridOperator& go, Solver& solver_, HELPER* helper_)
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go, helper_)
                , solver(solver_)
                , result_valid(false)
            {}

            void apply();

            void apply(TrialVector& u_);



            const Result& result() const
            {
                if (!result_valid)
                    DUNE_THROW(NewtonError,
                               "NewtonSolver::result() called before NewtonSolver::solve()");
                return this->res;
            }

        protected:
            virtual void defect(TestVector& r)
            {
                r = 0.0;                                        // TODO: vector interface
                this->gridoperator.residual(*this->u, r);
                TestVector rr = r;
                //this->helper->residual_change(rr);
                this->res.defect = this->solver.norm(rr);                    // TODO: solver interface
                if (!std::isfinite(this->res.defect))
                    DUNE_THROW(NewtonDefectError,
                               "NewtonHelperSolver::defect(): Non-linear defect is NaN or Inf");
            }


        private:
            template<class T>
            bool solutionControl(T& z) const
            {

                double znorm = z.infinity_norm()/static_cast<double>( z.N());
                if (std::isnan(znorm) || std::isinf(znorm) || std::abs(znorm)>1.e15)
                    return false;


                return true;
            }


            void linearSolve(Matrix& A, TrialVector& z, TestVector& r) const
            {
                if (this->verbosity_level >= 4)
                    std::cout << "      Solving linear system..." << std::endl;
                z = 0.0;                                        // TODO: vector interface
                this->solver.apply(A, z, r, this->linear_reduction);        // TODO: solver interface

                ios_base_all_saver restorer(std::cout); // store old ios flags

                if (!this->solver.result().converged)                 // TODO: solver interface
                    DUNE_THROW(NewtonLinearSolverError,
                               "NewtonHelperSolver::linearSolve(): Linear solver did not converge "
                               "in " << this->res.iterations << " iterations");
                if (!solutionControl(z))                 // TODO: solver interface
                    {
                        std::cout << "Bad solution in linear solver" << std::endl;
                    }

                if (this->verbosity_level >= 4)
                    std::cout << "          linear solver iterations:     "
                              << std::setw(12) << solver.result().iterations << std::endl
                              << "          linear defect reduction:      "
                              << std::setw(12) << std::setprecision(4) << std::scientific
                              << solver.result().reduction << std::endl;
            }

            Solver& solver;
            bool result_valid;
        };

        template<class GOS, class S, class TrlV, class TstV, class HELPER>
        void NewtonSolver<GOS,S,TrlV,TstV,HELPER>::apply(TrialVector& u_)
        {
            this->u = &u_;
            apply();
        }

        template<class GOS, class S, class TrlV, class TstV, class HELPER>
        void NewtonSolver<GOS,S,TrlV,TstV,HELPER>::apply()
        {
            this->res.iterations = 0;
            this->res.converged = false;
            this->res.reduction = 1.0;
            this->res.conv_rate = 1.0;
            this->res.elapsed = 0.0;
            this->res.assembler_time = 0.0;
            this->res.linear_solver_time = 0.0;
            this->res.linear_solver_iterations = 0;
            result_valid = true;
            Timer timer;

            try
                {
                    TestVector r(this->gridoperator.testGridFunctionSpace());
                    this->defect(r);
                    this->res.first_defect = this->res.defect;
                    this->prev_defect = this->res.defect;

                    if (this->verbosity_level >= 2)
                        {
                            // store old ios flags
                            ios_base_all_saver restorer(std::cout);
                            std::cout << "  Initial defect: "
                                      << std::setw(12) << std::setprecision(4) << std::scientific
                                      << this->res.defect << std::endl;
                        }

                    Matrix A(this->gridoperator);
                    TrialVector z(this->gridoperator.trialGridFunctionSpace());

                    //  this->helper.set_iteration(0);
                    while (!this->terminate())
                        {
                            if (this->verbosity_level >= 3)
                                std::cout << "  NewtonHelper iteration " << this->res.iterations
                                          << " --------------------------------" << std::endl;

                            Timer assembler_timer;
                            try
                                {
                                    this->prepare_step(A,r);
                                }
                            catch (...)
                                {
                                    this->res.assembler_time += assembler_timer.elapsed();
                                    throw;
                                }
                            double assembler_time = assembler_timer.elapsed();
                            this->res.assembler_time += assembler_time;
                            if (this->verbosity_level >= 3)
                                std::cout << "      matrix assembly time:             "
                                          << std::setw(12) << std::setprecision(4) << std::scientific
                                          << assembler_time << std::endl;

                            Timer linear_solver_timer;
                            try
                                {
                                    //   Dune::printvector(std::cout,r.base(),"residual","row",100,9,1);
                                    //  printmatrix(std::cout,A.base(),"global stiffness matrix","row",9,1);
                                    this->linearSolve(A, z, r);
                                    //  Dune::printvector(std::cout,z.base(),"update","row",100,9,1);
                                }
                            catch (...)
                                {
                                    this->res.linear_solver_time += linear_solver_timer.elapsed();
                                    this->res.linear_solver_iterations += this->solver.result().iterations;
                                    throw;
                                }
                            double linear_solver_time = linear_solver_timer.elapsed();
                            this->res.linear_solver_time += linear_solver_time;
                            this->res.linear_solver_iterations += this->solver.result().iterations;

                            if((this->res.iterations)>0)
                                this->helper->set_newton_iteration(this->res.iterations);
                            try
                                {
                                    this->line_search(z, r);
                                }
                            catch (NewtonLineSearchError)
                                {
                                    if (this->reassembled)
                                        throw;
                                    if (this->verbosity_level >= 3)
                                        std::cout << "      line search failed - trying again with reassembled matrix" << std::endl;
                                    continue;
                                }

                            this->res.reduction = this->res.defect/this->res.first_defect;
                            this->res.iterations++;
                            this->res.conv_rate = std::pow(this->res.reduction, 1.0/this->res.iterations);

                            // store old ios flags
                            ios_base_all_saver restorer(std::cout);
                            if (this->verbosity_level >= 3)
                                std::cout << "      linear solver time:               "
                                          << std::setw(12) << std::setprecision(4) << std::scientific
                                          << linear_solver_time << std::endl
                                          << "      defect reduction (this iteration):"
                                          << std::setw(12) << std::setprecision(4) << std::scientific
                                          << this->res.defect/this->prev_defect << std::endl
                                          << "      defect reduction (total):         "
                                          << std::setw(12) << std::setprecision(4) << std::scientific
                                          << this->res.reduction << std::endl
                                          << "      new defect:                       "
                                          << std::setw(12) << std::setprecision(4) << std::scientific
                                          << this->res.defect << std::endl;
                            if (this->verbosity_level == 2)
                                std::cout << "  Newton iteration " << std::setw(2) << this->res.iterations
                                          << ".  New defect: "
                                          << std::setw(12) << std::setprecision(4) << std::scientific
                                          << this->res.defect
                                          << ".  Reduction (this): "
                                          << std::setw(12) << std::setprecision(4) << std::scientific
                                          << this->res.defect/this->prev_defect
                                          << ".  Reduction (total): "
                                          << std::setw(12) << std::setprecision(4) << std::scientific
                                          << this->res.reduction << std::endl;

                        }
                }
            catch(...)
                {
                    this->res.elapsed = timer.elapsed();
                    throw;
                }
            this->res.elapsed = timer.elapsed();
            //write final number of iterations
            this->helper->set_newton_iteration(-this->res.iterations);
            ios_base_all_saver restorer(std::cout); // store old ios flags

            if (this->verbosity_level == 1)
                std::cout << "  NewtonHelper converged after " << std::setw(2) << this->res.iterations
                          << " iterations.  Reduction: "
                          << std::setw(12) << std::setprecision(4) << std::scientific
                          << this->res.reduction
                          << "   (" << std::setprecision(4) << this->res.elapsed << "s)"
                          << std::endl;
        }

        template<class GOS, class TrlV, class TstV, class HELPER>
        class NewtonTerminate : public virtual NewtonBase<GOS,TrlV,TstV,HELPER>
        {
            typedef GOS GridOperator;
            typedef TrlV TrialVector;

            typedef typename TstV::ElementType RFType;

        public:
            NewtonTerminate(GridOperator& go, TrialVector& u_, HELPER* helper_)
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go,u_, helper_)
                , reduction(1e-8)
                , maxit(40)
                , force_iteration(false)
                , abs_limit(1e-12)
            {}

            NewtonTerminate(GridOperator& go, HELPER* helper_)
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go,helper_)
                , reduction(1e-8)
                , maxit(40)
                , force_iteration(false)
                , abs_limit(1e-12)
            {}

            void setReduction(RFType reduction_)
            {
                reduction = reduction_;
            }

            void setMaxIterations(unsigned int maxit_)
            {
                maxit = maxit_;
            }

            void setForceIteration(bool force_iteration_)
            {
                force_iteration = force_iteration_;
            }

            void setAbsoluteLimit(RFType abs_limit_)
            {
                abs_limit = abs_limit_;
            }

            virtual bool terminate()
            {
                if (force_iteration && this->res.iterations == 0)
                    return false;
                this->res.converged = this->res.defect < abs_limit
                    || this->res.defect < this->res.first_defect * reduction;
                if (this->res.iterations >= maxit && !this->res.converged)
                    DUNE_THROW(NewtonNotConverged,
                               "NewtonHelperTerminate::terminate(): Maximum iteration count reached");
                return this->res.converged;
            }

        private:
            RFType reduction;
            unsigned int maxit;
            bool force_iteration;
            RFType abs_limit;
        };

        template<class GOS, class TrlV, class TstV, class HELPER>
        class NewtonPrepareStep : public virtual NewtonBase<GOS,TrlV,TstV,HELPER>
        {
            typedef GOS GridOperator;
            typedef TrlV TrialVector;

            typedef typename TstV::ElementType RFType;
            typedef typename GOS::template MatrixContainer<RFType>::Type Matrix;

        public:
            NewtonPrepareStep(GridOperator& go, TrialVector& u_, HELPER* helper_)
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go,u_, helper_)
                , min_linear_reduction(1e-3)
                , reassemble_threshold(0.0)
            {}

            NewtonPrepareStep(GridOperator& go, HELPER* helper_)
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go, helper_)
                , min_linear_reduction(1e-3)
                , reassemble_threshold(0.0)
            {}

            void setMinLinearReduction(RFType min_linear_reduction_)
            {
                min_linear_reduction = min_linear_reduction_;
            }

            void setReassembleThreshold(RFType reassemble_threshold_)
            {
                reassemble_threshold = reassemble_threshold_;
            }

            virtual void prepare_step(Matrix& A, TstV& )
            {
                this->reassembled = false;
                if (this->res.defect/this->prev_defect > reassemble_threshold)
                    {
                        if (this->verbosity_level >= 3)
                            std::cout << "      Reassembling matrix..." << std::endl;
                        A = 0.0;                                    // TODO: Matrix interface
                        this->gridoperator.jacobian(*this->u, A);
                        this->reassembled = true;
                    }

                this->linear_reduction = std::min(min_linear_reduction,
                                                  this->res.defect*this->res.defect/
                                                  (this->prev_defect*this->prev_defect));

                this->prev_defect = this->res.defect;

                ios_base_all_saver restorer(std::cout); // store old ios flags

                if (this->verbosity_level >= 4)
                    std::cout << "      requested linear reduction:       "
                              << std::setw(12) << std::setprecision(4) << std::scientific
                              << this->linear_reduction << std::endl;
            }

        private:
            RFType min_linear_reduction;
            RFType reassemble_threshold;
        };

        template<class GOS, class TrlV, class TstV, class HELPER>
        class NewtonLineSearch : public virtual NewtonBase<GOS,TrlV,TstV,HELPER>
        {
            typedef GOS GridOperator;
            typedef TrlV TrialVector;
            typedef TstV TestVector;

            typedef typename TestVector::ElementType RFType;

        public:
            enum Strategy { noLineSearch,
                            hackbuschReusken,
                            hackbuschReuskenAcceptBest };

            NewtonLineSearch(GridOperator& go, TrialVector& u_, HELPER* helper_)
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go,u_,helper_)
                , strategy(hackbuschReusken)
                , maxit(10)
                , damping_factor(0.5)
            {}

            NewtonLineSearch(GridOperator& go, HELPER* helper_)
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go,helper_)
                , strategy(hackbuschReusken)
                , maxit(10)
                , damping_factor(0.5)
            {}

            void setLineSearchStrategy(const std::string strategy_)
            {
                strategy = strategyFromName(strategy_);
            }


            void setLineSearchMaxIterations(unsigned int maxit_)
            {
                maxit = maxit_;
            }

            void setLineSearchDampingFactor(RFType damping_factor_)
            {
                damping_factor = damping_factor_;
            }

            virtual void line_search(TrialVector& z, TestVector& r)
            {
                if (strategy == noLineSearch)
                    {
                        this->u->axpy(-1.0, z);                     // TODO: vector interface
                        this->defect(r);
                        return;
                    }

                if (this->verbosity_level >= 4)
                    std::cout << "      Performing line search..." << std::endl;
                RFType lambda = 1.0;
                RFType best_lambda = 0.0;
                RFType best_defect = this->res.defect;
                TrialVector prev_u(*this->u);  // TODO: vector interface
                unsigned int i = 0;
                ios_base_all_saver restorer(std::cout); // store old ios flags

                while (1)
                    {
                        if (this->verbosity_level >= 4)
                            std::cout << "          trying line search damping factor:   "
                                      << std::setw(12) << std::setprecision(4) << std::scientific
                                      << lambda;
                        // << std::endl;

                        this->u->axpy(-lambda, z);                  // TODO: vector interface
                        try {
                            this->defect(r);
                            if (this->verbosity_level >= 4)
                                std::cout << " new defect " << this->res.defect << std::endl;
                        }
                        catch (NewtonDefectError)
                            {
                                if (this->verbosity_level >= 4)
                                    std::cout << "          Nans detected" << std::endl;
                            }       // ignore NaNs and try again with lower lambda

                        // write linesearch iteration and vector to helper
                        this->helper->ls_iteration(i,(*this->u));
                        if (this->res.defect <= (1.0 - lambda/4) * this->prev_defect)
                            {
                                if (this->verbosity_level >= 4)
                                    std::cout << "          line search converged" << std::endl;
                                break;
                            }

                        if (this->res.defect < best_defect)
                            {
                                best_defect = this->res.defect;
                                best_lambda = lambda;
                            }

                        if (++i >= maxit)
                            {
                                if (this->verbosity_level >= 4)
                                    std::cout << "          max line search iterations exceeded" << std::endl;
                                switch (strategy)
                                    {
                                    case hackbuschReusken:
                                        *this->u = prev_u;
                                        this->defect(r);
                                        DUNE_THROW(NewtonLineSearchError,
                                                   "NewtonHelperLineSearch::line_search(): line search failed, "
                                                   "max iteration count reached, "
                                                   "defect did not improve enough");
                                    case hackbuschReuskenAcceptBest:
                                        if (best_lambda == 0.0)
                                            {
                                                *this->u = prev_u;
                                                this->defect(r);
                                                DUNE_THROW(NewtonLineSearchError,
                                                           "NewtonHelperLineSearch::line_search(): line search failed, "
                                                           "max iteration count reached, "
                                                           "defect did not improve in any of the iterations");
                                            }
                                        if (best_lambda != lambda)
                                            {
                                                *this->u = prev_u;
                                                this->u->axpy(-best_lambda, z);
                                                this->defect(r);
                                            }
                                        break;
                                    case noLineSearch:
                                        break;
                                    }
                                break;
                            }

                        lambda *= damping_factor;
                        *this->u = prev_u;                          // TODO: vector interface
                    }
                if (this->verbosity_level >= 4)
                    std::cout << "          line search damping factor:   "
                              << std::setw(12) << std::setprecision(4) << std::scientific
                              << lambda << std::endl;
            }

        protected:
            /** helper function to get the different strategies from their name */
            Strategy strategyFromName(const std::string & s) {
                if (s == "noLineSearch")
                    return noLineSearch;
                if (s == "hackbuschReusken")
                    return hackbuschReusken;
                if (s == "hackbuschReuskenAcceptBest")
                    return hackbuschReuskenAcceptBest;
                DUNE_THROW(NewtonError,"linear search strategy " << s << " not known");
                return noLineSearch;
            }

        private:
            Strategy strategy;
            unsigned int maxit;
            RFType damping_factor;
        };

        template<typename TrlV>
        NewtonHelperBase<TrlV> getHelperClass() {
            static NewtonHelperBase<TrlV> nhc;
            return  nhc;
        }

        template<class GOS, class S, class TrlV, class TstV = TrlV, class HELPER = NewtonHelperBase<TrlV> >
        class Newton : public NewtonSolver<GOS,S,TrlV,TstV,HELPER>
                     , public NewtonTerminate<GOS,TrlV,TstV,HELPER>
                     , public NewtonLineSearch<GOS,TrlV,TstV,HELPER>
                     , public NewtonPrepareStep<GOS,TrlV,TstV,HELPER>
        {
            typedef GOS GridOperator;
            typedef S Solver;
            typedef TrlV TrialVector;

        public:
            Newton(GridOperator& go, TrialVector& u_, Solver& solver_, HELPER* helper_ = new HELPER())
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go,u_,helper_)
                , NewtonSolver<GOS,S,TrlV,TstV,HELPER>(go,u_,solver_,helper_)
                , NewtonTerminate<GOS,TrlV,TstV,HELPER>(go,u_,helper_)
                , NewtonLineSearch<GOS,TrlV,TstV,HELPER>(go,u_,helper_)
                , NewtonPrepareStep<GOS,TrlV,TstV,HELPER>(go,u_,helper_)
            {}
            Newton(GridOperator& go, Solver& solver_, HELPER* helper_= new HELPER() )
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go,helper_)
                , NewtonSolver<GOS,S,TrlV,TstV,HELPER>(go,solver_,helper_)
                , NewtonTerminate<GOS,TrlV,TstV,HELPER>(go,helper_)
                , NewtonLineSearch<GOS,TrlV,TstV,HELPER>(go,helper_)
                , NewtonPrepareStep<GOS,TrlV,TstV,HELPER>(go,helper_)
            {}


            Newton(GridOperator& go, TrialVector& u_, Solver& solver_, HELPER& helper_)
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go,u_,&helper_)
                , NewtonSolver<GOS,S,TrlV,TstV,HELPER>(go,u_,solver_,&helper_)
                , NewtonTerminate<GOS,TrlV,TstV,HELPER>(go,u_,&helper_)
                , NewtonLineSearch<GOS,TrlV,TstV,HELPER>(go,u_,&helper_)
                , NewtonPrepareStep<GOS,TrlV,TstV,HELPER>(go,u_,&helper_)
            {}
            Newton(GridOperator& go, Solver& solver_, HELPER& helper_)
                : NewtonBase<GOS,TrlV,TstV,HELPER>(go,&helper_)
                , NewtonSolver<GOS,S,TrlV,TstV,HELPER>(go,solver_,&helper_)
                , NewtonTerminate<GOS,TrlV,TstV,HELPER>(go,&helper_)
                , NewtonLineSearch<GOS,TrlV,TstV,HELPER>(go,&helper_)
                , NewtonPrepareStep<GOS,TrlV,TstV,HELPER>(go,&helper_)
            {}

        };


    }
}

#endif // DUNE_PDELAB_NEWTON_HELPER_HH
