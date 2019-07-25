# Problem: 
### Using different numerical methods to calcualte PDE.
=============================================================

# Symbol: 

- x: space position
- u1, u2: same functions in pdf
- t: time
=============================================================

#Thinking in code:

Choose `xstep(N)`, `cfl` and then make loops about time.
1. initialize u10, u20;
2. calculate accurate results based on u(x,t)=u0(x-at);
3. initialize all these numerical method using periodic condition;
4. calculate them:
- Crank Nicolson and Global finite difference are calculated in matrix;
- Upwind, leap frog and lax wendroff are not difficult;
- Lax Friedrich requires some temporary variable (first step is Upwind).
- make plots.

