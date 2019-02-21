# temperature-distributions
Modeling the heat distribution of a 2d metal square at equilibrium with known boundary temperatures

### Background, described in 'n' dimensions:

The flow of heat through an n-dimensional object can be determined by a partial differential equation called the [heat equation](https://en.wikipedia.org/wiki/Heat_equation "Wikipedia - Heat Equation"):

![du/dt - α∇²u = 0](https://wikimedia.org/api/rest_v1/media/math/render/svg/3edc07e9067b68e6057723653f7c3e7403889598 "the heat equation")

If there are no external sources or sinks of heat, then eventually the temperature at any point in the object will converge on a constant temperature (this will be the actual temperature of the point after an infinite amount of time, but after finite amounts of time the value will become asymptotically similar to the eventual value).

Then, &nbsp;<img src="figures/frac_del-u_del-t.png" height="36" alt="du/dt = 0" title="">

so, &nbsp;<img src="figures/nabla2_u_0.png" height="18" alt="∇²u = 0" title="">

This is [Laplace's equation](https://en.wikipedia.org/wiki/Laplace%27s_equation "Laplace's equation"), and therefore the final (or equilibrium) temperature distribution is determined by Laplace's equation.
The solutions to Laplace's equation are called the [harmonic functions](https://en.wikipedia.org/wiki/Harmonic_function "Harmonic Function"), and so the function u that describes the equilibrium heat distribution will be a harmonic function.

More specifically, if the temperatures at the boundary of the n-dimensional object are known, the problem of finding the equilibrium temperature at any point in the n-dimensional object whose boundary values are known is an example of the [Dirichlet problem](https://en.wikipedia.org/wiki/Dirichlet_problem "Dirichlet problem").

The solution to a Dirichlet problem is not only a harmonic function (as all solutions to Laplace's equation are), but is also unique, meaning that there is only one possible equilibrium heat distribution over a given n-dimensional object whose boundary temperatures are known.

One useful property of the harmonic functions is that [they satisfy the mean value property](https://en.wikipedia.org/wiki/Harmonic_function#The_mean_value_property "Wikipedia - Harmonic function : The mean value property"): for any point, the temperature must be the average of the temperatures at nearby points.
The [continuous mean value property](https://sites.math.washington.edu/~morrow/336_13/mvp1.pdf "University of Washington : Math - Mean Value Property") describes this situation perfectly: The value at any point P in the n-dimensional object will be the average of all equidistant points a very small distance r from P. Considering the infinitely many nearby points as r becomes infinitely small, for every point P, we can obtain a continuous function over the n-dimensional object.

Finding a perfect continuous function that satisfies our Dirichlet problem would take some very difficult calculus and be very computationally expensive. Even in lower dimensions, we'd like to avoid having to calculate it.
Consider, instead, a *discrete* mean value property, letting the temperature at any point be the algebraic mean of the temperatures of a finite number of nearby points.

This can be done with an n-dimensional grid, letting the temperature at any point be the average of all points point that are exactly one unit away on the grid.

This proves to be a shockingly good approximation even with a very small number of points, and as the distance between points is decreased (tending towards 0), the values found with the discrete mean value property tend towards the values of the continuous one.

## The Problem at Hand:
For simplicity (and the sake of your CPU) we'll consider a 2-dimensional conductive object, ie. a very thin metal plate.
We'll assume the plate to be a unit square, and take the temperature on the left half of the plate's boundary to be 0°, while the temperature on the right half of the plate's boundary will be 1°.
The discrete mean value property can be applied to an imaginary grid laid across our plate such that the temperature at each point on the grid will be the average of the four nearest points, such that

<img src="figures/f_pxy_equals.png" height="36" alt="f(P(x,y)) = 1/4 ( f(P(x+1,y)) + f(P(x-1,y)) + f(P(x,y+1)) + f(P(x-1,y-1)) )">

### A 2x2 Matrix
The previous function leads to the set of linear equations for a four-point grid, numbered left-right top-bottom, to be

<img src="figures/t1_lin_sys.png" height="24" alt="t<sub>1</sub>= (1/4) (t<sub>2</sub> + t<sub>3</sub> + 0 + 0)" title="">

<img src="figures/t2_lin_sys.png" height="24" alt="t<sub>2</sub>= (1/4) (t<sub>1</sub> + t<sub>4</sub> + 1 + 1)" title="">

<img src="figures/t3_lin_sys.png" height="24" alt="t<sub>3</sub>= (1/4) (t<sub>1</sub> + t<sub>4</sub> + 0 + 0)" title="">

<img src="figures/t4_lin_sys.png" height="24" alt="t<sub>4</sub>= (1/4) (t<sub>2</sub> + t<sub>3</sub> + 1 + 1)" title="">

this is equivalent to

<img src="figures/t_exp_lin_sys.png" height="72" alt="t=Mt+b expanded" title="">

or simply

<img src="figures/t_smol_lin_sys.png" height="30" alt="t=Mt+b expanded" title="">

solving for t, this leaves

<img src="figures/t_I-M_-1.png" height="30" alt="t=(I-M)<sup>-1</sup> b" title="">


If you choose to calculate this out by hand (I used Guass-Jordan elimination), you'll find

<img src="figures/I-M_inv_eq.png" height="72" alt="what (I-M)<sup>-1</sup> equals" title="">

which determines the temperature approximations at our four points to be t<sub>1</sub>=t<sub>3</sub>=0.25 and t<sub>2</sub>=t<sub>4</sub>=0.75.

Ask your friendly neighborhood physics PhD and they could tell you that the exact temperatures in reality are t<sub>1</sub>=t<sub>3</sub>=0.2371 and t<sub>2</sub>=t<sub>4</sub>=0.7629. That means that with only 4 grid points, our margin of error from the real temperatures was only about 5%.
And if we increase the number of grid points, our approximations will only get better.

### The need for automation
Perhaps you, like me, just absolutely refuse to do the Gauss-Jordan elimination process by hand on anything more complex than the previous 2x2 matrix.
Thankfully, [Jacobi's method](https://www.maa.org/press/periodicals/loci/joma/iterative-methods-for-solving-iaxi-ibi-jacobis-method "Mathematical Association of America - Iterative Methods for Solving Ax = b - Jacobi's Method") exists.
Jacobi's method is an iterative method for finding x when Ax=b, and A is a diagonally dominant matrix, where the approximation as the iteration number increases converges on the true solution.
Jacobi's method can be expressed as

<img src="https://www.maa.org/sites/default/files/images/cms_upload/JacobisMethod435220.gif" height="30" alt="x^k+1 = D^-1 ((-L-U) x^k + b)" title="Jacobi's method">

Where D is the matrix with only the diagonal component of A, L and U are the lower and upper triangular portions of A, and the superscript of x is the iteration number, not a power operator.

Our original matrix equation is actually very similar: where D is the identity matrix, M can step in for (-L-U), and I+M would be a diagonally dominant matrix. Therefore, to avoid giving the computer any further work, we can simply use our original equation modified only so that the lefthand side simply being the next iteration:

<img src="figures/t_k+1.png" height="30" alt="t^k+1 = Mt^k + b" title="">

This is now in a form that can easily be looped over with numpy:

```python
t = numpy.matmul(M, t) + b
```

Try this on our 2x2 matrix with an initial t<sup>0</sup> set to the 0 vector
and after just <a href="figures/5_iter.png" >5 iterations</a> we have the values t<sub>1</sub>=t<sub>3</sub>=0.2344 and t<sub>2</sub>=t<sub>4</sub>=0.7344, which are already great approximations of the t vector we found earlier for the 2x2 grid.
