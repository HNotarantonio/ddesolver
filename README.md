# **INTRODUCTION**

We introduce **ddesolver**, a Maple package designed for solving *Discrete Differential
Equations (DDEs)*. These equations are functional equations relating algebraically a formal
power series $F(t, u)$ with polynomial coefficients in a “catalytic” variable $u$, with
specializations of it with respect to the catalytic variable. Such equations appear in
enumerative combinatorics, for instance in the enumeration of maps. Bousquet-Mélou and
Jehanne showed in $2006$ that when these equations are of a fixed point type in $F$, then $F$ is
an *algebraic series* (that is, there exists $R_u\in\mathbb{Q}[X, t, u]\setminus\{0\}$ such that $R_u(F(t, u), t, u)=0$).
In the same paper, they proposed a systematic method for computing
annihilating polynomials of these series: note that the usual interest lies in the algebraicity of 
one of the specialized series (say of $F(t, 1)$). Bostan, Safey El Din and the author of this package
designed in $2023$ new efficient algorithms for computing these witnesses of algebraicity. 

The package **ddesolver** contains the implementations of these new algorithms.

Want to know more about ddesolver? Have a look at the dedicated [tutorial paper](https://mathexp.eu/notarantonio/papers/ddesolver.pdf)!   


# **INSTALLATION**

ddesolver is available in the “.mla” format.

The Maple variable libname shall be set so that “ddesolver.mla” is
located in a visible place.

```
libname := ”/home/notarantonio/ddesolver/lib”, libname:
```

Once libname has been correctly set up, one executes in Maple

```
with(ddesolver);
```

in order to load and use the package.


# **MOTIVATING EXAMPLE**

  Consider the DDE associated to the enumeration of $3$-constellations:
    
  ```math
     \begin{align*} F(t, u) = 1 + tuF(t, u)^3 &+ tu(2F(t, u) + F(t, 1))\frac{F(t, u)-F(t, 1)}{u-1}\\
           &+ tu\frac{F(t, u) - F(t, 1) - (u-1)\partial_uF(t, 1)}{(u-1)^2}.
     \end{align*}
  ```
            
  Multiplying the above functional equation by $(u-1)^2$ yields
  
  ```math
     \begin{align*}  0=(u-1)^2(1-F(t, u) + tuF(t, u)^3) &+ tu(u-1)(2F(t, u) + F(t, 1))(F(t, u)-F(t, 1))\\
             &+  tu(F(t, u) - F(t, 1) - (u-1)\partial_uF(t, 1)).
     \end{align*}
  ```
      
  The above rewrites in the form $P(F(t, u), F(t, 1), \partial_uF(t, 1), t, u)=0$,
  where $P\in\mathbb{Q}[x, z_0, z_1, t, u]$ is given by
  
  ```math
    P := (u-1)^2(1-x+tux^3) +tu(u-1)(2x+z_0)(x-z_0)+tu(x-z_0-(u-1)z_1).
  ```
    
  We continue the analysis with Maple
    
    P := (u-1)^2(1-x+tux^3) +tu(u-1)(2x+z0)(x-z0)+tu(x-z0-(u-1)z1):
    
    with(ddesolver):

    annihilating_polynomial(P, 2);
 
    (16tz0^2-8tz0+t-16)(81t^2z0^3-81t^2z0^2+27t^2z0+18tz0^2-3t^2-66tz0+47t+z0-1)

  Thus $R := (16tz_0^2-8tz_0+t-16)(81t^2z_0^3-81t^2z_0^2+27t^2z_0+18tz_0^2-3t^2-66
    tz_0+47t+z_0-1)$ annihilates the series $F(t, 1)$.
 

