#ifndef __NEWTON_IMP_H__
#define __NEWTON_IMP_H__


/*!
 \param x0 the initial guess for the optimization method
 \param tolerance the tolerance used as stopping criterion for the iterative optimization method
 \param max_iter the maximum number of iterations
 \param ch a reference to a Checker object, used to set the reason of termination of the iterations.
 \param GCV_v a reference to the vector of GCV values evaluated during the iterative procedure
 \param lambda_v a reference to the vector of lambda values explored during the iterative procedure
 \return std::pair<Tuple, UInt>, a pair which containns the optimal lambda found and the number of iterations to reach the tolerance
*/
template <typename Tuple, typename Hessian, typename ...Extensions>
std::pair<Tuple, UInt> Newton_ex<Tuple, Hessian, Extensions...>::compute (const Tuple & x0, const Real tolerance, const UInt max_iter, Checker & ch, std::vector<Real> & GCV_v, std::vector<Real> & lambda_v)
{
       // Initialize the algorithm
       Tuple x_old;
       Tuple x      = x0;
       UInt  n_iter = 0;
       Real  error  = std::numeric_limits<Real>::infinity();

       // Debugging purpose
       // Rprintf("\n Starting Initializing lambda phase\n");

       /* Start from 6 lambda and find the minimum value of GCV to start from it the newton's method*/
       Real valmin, valcur, lambda_min;
       UInt Nm = 6;
       std::vector<Real> vals={5.000000e-05, 1.442700e-03, 4.162766e-02, 1.201124e+00, 3.465724e+01, 1.000000e+03};
       valcur = this->F.evaluate_f(vals[0]);
       lambda_min = 5e-5;
       valmin = valcur;

       for(UInt i=1; i<Nm; i++)
       {
               valcur = this->F.evaluate_f(vals[i]);
               if(valcur<valmin)
               {
                       valmin = valcur;
                       lambda_min = vals[i];
               }
       }

       if(x>lambda_min/4 || x<=0)
       {
               x = lambda_min/8;
       }

       Rprintf("\n Starting Newton's iterations: starting point lambda=%f\n",x);

       //only the first time applied here
       Real    fx  = this->F.evaluate_f(x);
       Tuple   fpx = this->F.evaluate_first_derivative (x);
       Hessian fsx = this->F.evaluate_second_derivative(x);

       while(n_iter < max_iter)
       {
               // Debugging purpose f(x)
               GCV_v.push_back(fx);
               lambda_v.push_back(x);

               if(Auxiliary<Tuple>::isNull(fsx))
               {
                       // Debug message
                       // std::cout << "Division by zero" << std::endl;
                       return {x, n_iter};
               }

               ++n_iter;


               x_old = x;
               Auxiliary<Tuple>::divide(fsx, fpx, x);
               x = x_old - x;

               if (x<=0)
               {
                       Rprintf("\nProbably monotone increasing GCV function\n");

                       fx = this->F.evaluate_f(x);
                       return {x_old, n_iter};
               }

               // Put here the updates in order to compute error on the correct derivative and to have z_hat updated for the solution

               fpx = this->F.evaluate_first_derivative (x);

               error = Auxiliary<Tuple>::residual(fpx);

               Rprintf("\nStep number %d  of EXACT-NEWTON: residual = %f\n", n_iter, error);

               if(error<tolerance)
               {
                       /* fpx=this->F.evaluate_f(x-0.5);
                       fsx=this->F.evaluate_f(x+0.5);
                       if (std::abs(fpx-fsx)/fsx<=0.01)
                               Rprintf("GCV has a non standard shape");*/

                       ch.set_tolerance();
                       fx = this->F.evaluate_f(x);
                       return {x, n_iter};
               }

               fx  = this->F.evaluate_f(x);
               fsx = this->F.evaluate_second_derivative(x);
       }

       fx = this->F.evaluate_f(x);
       ch.set_max_iter();
       return {x, n_iter};
}


/*!
 \param x0 the initial guess for the optimization method
 \param tolerance the tolerance used as stopping criterion for the iterative optimization method
 \param max_iter the maximum number of iterations
 \param ch a reference to a Checker object, used to set the reason of termination of the iterations.
 \param GCV_v a reference to the vector of GCV values evaluated during the iterative procedure
 \param lambda_v a reference to the vector of lambda values explored during the iterative procedure
 \return std::pair<Real, UInt>, a pair which containns the optimal lambda found and the number of iterations to reach the tolerance
*/
template <typename ...Extensions>
std::pair<Real, UInt> Newton_fd<Real, Real, Extensions...>::compute (const Real & x0, const Real tolerance, const UInt max_iter, Checker & ch, std::vector<Real> & GCV_v, std::vector<Real> & lambda_v)
{
        // Initialize the algorithm
        Real x_old;
        Real x      = x0;
        UInt n_iter = 0;
        Real error  = std::numeric_limits<Real>::infinity();
        Real h      = 4e-6;

        // Debugging purpose
        // Rprintf("\n Starting Initializing lambda phase"); /*! Start from 6 lambda and find the minimum value of GCV to start from it the newton's method*/

        Real valmin, valcur, lambda_min;
        UInt Nm = 6;
        std::vector<Real> vals = {5.000000e-05, 1.442700e-03, 4.162766e-02, 1.201124e+00, 3.465724e+01, 1.000000e+03};
        valcur=this->F.evaluate_f(vals[0]);
        lambda_min = 5e-5;
        valmin = valcur;

        for(UInt i=1; i<Nm; i++)
        {
                valcur = this->F.evaluate_f(vals[i]);

                if(valcur<valmin)
                {
                        valmin = valcur;
                        lambda_min = vals[i];
                }
        }


        if (x>lambda_min/4 || x<=0)
        {
                x = lambda_min/8;
        }
        Rprintf("\n Starting Newton's iterations: starting point lambda=%f\n",x);

        // Only the first time applied here
        // Rprintf("Forward: \n");
        Real fxph = this->F.evaluate_f(x+h);
        //Rprintf("Backward: \n");
        Real fxmh = this->F.evaluate_f(x-h);
        // Rprintf("Center: \n");
        Real fx  = this->F.evaluate_f(x);

        Real fpx = (fxph-fxmh)/(2*h);
        // Rprintf("fp(x): %f\n", fpx);

        Real fsx = (fxph+fxmh-(2*fx))/(h*h);
        // Rprintf("fs(x): %f\n", fsx);

        while(n_iter < max_iter)
        {
                GCV_v.push_back(fx);
                lambda_v.push_back(x);

                //Debugging purpose f(x)
                if (Auxiliary<Real>::isNull(fsx))
                {
                        // Debug message
                        // std::cout << "Division by zero" << std::endl;
                        return {x, n_iter};
                }

                ++n_iter;


                x_old = x;
                Auxiliary<Real>::divide(fsx, fpx, x);
                x = x_old - x;

                if (x<=0)
                { //too small value
                        Rprintf("\nProbably monotone increasing GCV function\n");

                        fx = this->F.evaluate_f(x_old);
                        return {x_old, n_iter};
                }

                // Put here the updates in order to compute error on the correct derivative and to have z_hat updated for the solution
                // Rprintf("Forward:\n");
                fxph = this->F.evaluate_f(x+h);
                // Rprintf("Backward: \n");
                fxmh = this->F.evaluate_f(x-h);

                fpx = (fxph-fxmh)/(2*h);

                error = Auxiliary<Real>::residual(fpx);

                Rprintf("\nStep number %d  of FD-NEWTON: residual = %f\n", n_iter, error);

                if (error < tolerance)
                {
                        ch.set_tolerance();
                        fx  = this->F.evaluate_f(x); //eventuale miglioramento: va fatto altirmenti prende gli z:hat di quellos sbagliato.
                        return {x, n_iter};
                }

                // Rprintf("Center: \n");
                fx  = this->F.evaluate_f(x);

                fsx = (fxph+fxmh-(2*fx))/(h*h);

                // Rprintf("fp(x): %f\n", fpx);
                // Rprintf("fs(x): %f\n", fsx);
        }
        fx  = this->F.evaluate_f(x);
        ch.set_max_iter();
        return {x, n_iter};
}

#endif
