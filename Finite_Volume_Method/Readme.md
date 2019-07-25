# Solve the PDE problem

- Accurate results: refer to the photos
- Finite Volume method: refer to the Matlab Code

======================================================================

# Way of thinking:
- initialization
1. t* = 2 (see photos).
2. Finite volume method: integrate about x for both sides of equation.
3. Calculate `u_ave` for each x step.
- evolution
4. using 3-order Runge-Kutta calculate for the next time step (periodic condition).
- reconsitution 
5. use u_ave to calculate u at endpoints.

======================================================================

There are some discussions about how to choose the analitical value.
You can refer to the photos.
