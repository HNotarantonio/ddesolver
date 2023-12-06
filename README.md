# **INTRODUCTION**

We introduce **ddesolver**, a Maple package designed for solving *Discrete Differential
Equations (DDEs)*. These equations are functional equations relating algebraically a formal
power series $F(t, u)$ with polynomial coefficients in a “catalytic” variable $u$, with
specializations of it with respect to the catalytic variable. Such equations appear in
enumerative combinatorics, for instance in the enumeration of maps. [Mireille Bousquet-Mélou](https://www.labri.fr/perso/bousquet/) and
[Arnaud Jehanne](https://www.math.u-bordeaux.fr/~ajehanne/) showed in $2006$ that when these equations are of a fixed point type in $F$, then $F$ is
an *algebraic series* (that is, there exists a nonzero $R_u\in\mathbb{Q}[X, t, u]$ such that $R_u(F(t, u), t, u)=0$).
In the same paper, they proposed a systematic method for computing
annihilating polynomials of these series: note that the usual interest lies in the algebraicity of 
one of the specialized series, say of $F(t, 1)$. [Alin Bostan](https://mathexp.eu/bostan/), [Mohab Safey El Din](https://www-polsys.lip6.fr/~safey/) and the author of this package
designed in $2023$ new efficient algorithms for computing these witnesses of algebraicity. 

The package **ddesolver** contains the implementations of these new algorithms.

Want to know more about ddesolver? Have a look at the dedicated [tutorial paper](https://mathexp.eu/notarantonio/papers/ddesolver.pdf)!   
Want to more about the algorithms? Have a look [here](https://mathexp.eu/notarantonio/papers/reviewed_BoNoSa2023--HAL-ARXIV.pdf)!

# **INSTALLATION**

ddesolver is a [Maple](https://fr.maplesoft.com/) package available in the “.mla” format.

The Maple variable [libname](https://fr.maplesoft.com/support/help/maple/view.aspx?path=libname&L=F) shall be set so that “ddesolver.mla” is
located in a visible place.

```
libname := ”PATH_TO/ddesolver/lib”, libname:
```

Once libname has been correctly set up, one executes in Maple

```
with(ddesolver);
```

in order to load and use the package.

# **INPUT/OUTPUT FORMAT**
The function **annihilating_polynomial** takes as input $(P, k)$, where:

- $P\in\mathbb{Q}[x, z_0, z_1, \ldots, z_{k-1}, t, u]$
is a nonzero polynomial such that
```math
P(F(t, u), F(t, 1), \partial_uF(t, 1), \ldots, \partial_u^{k-1}F(t, 1), t, u)=0,
```
- $k$ is the number of specialized series.

Note that the specialization of the series can be at any point of $a\in\mathbb{Q}$. Also, it is important that the series of interest $F(t, 1)$ corresponds to $z_0$ and that the full series $F(t, u)$ corresponds to $x$.

The function **annihilating_polynomial** outputs a nonzero polynomial $R\in\mathbb{Q}[t, z_0]$
such that $R(t, F(t, 1))=0$. 

Note that in order to have a nonzero polynomial $R\in\mathbb{Q}[t, z_0]$ for $F(t, 1)$, it is mandatory that the assumptions related to Section $4.2$ in the [tutorial paper](https://mathexp.eu/notarantonio/papers/ddesolver.pdf)
are full filled. 

# **EXAMPLE**

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
 
# **OPTIONS**

It is possible to benefit in practice from two options: the choice of the `algorithm`, and the
choice of the `variable` (either $t$ or $z_0$) on which we perform evaluation–interpolation

**Choice of the algorithm:** Four algorithms (described in Section $4$ of the [tutorial paper](https://mathexp.eu/notarantonio/papers/ddesolver.pdf)) are implemented: "duplication", "elimination",
"geometry" (only when $k = 2$) and "hybrid". By default, the algorithm used is "elimination".

**Choice of the variable:** The two choices are $t$ and $z_0$. The default choice is $t$.
The choice of the algorithm and of the variable on which we perform evaluation–
interpolation can be made by executing `annihilating_polynomial(P , k, algorithm, variable)`.
Note that these two options must be either not specified at all, or specified in the same call.

**Example:** We continue with the enumeration of $3$-constellations, this time with the options `algorithm = "geometry"` and `variable=z0`.

    annihilating_polynomial(P , 2, "geometry" , z0);
    (16tz0^2-8tz0+t-16)(81t^2z0^3-81t^2z0^2+27t^2z0+18tz0^2-3t^2-66tz0+47t+z0-1)
    
