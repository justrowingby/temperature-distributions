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
