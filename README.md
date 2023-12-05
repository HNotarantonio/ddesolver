# **INTRODUCTION**

We introduce DDE-Solver, a Maple package designed for solving Discrete Differential
Equations (DDEs). These equations are functional equations relating algebraically a formal
power series F(t, u) with polynomial coefficients in a “catalytic” variable u, with special-
izations of it with respect to the catalytic variable. Such equations appear in
enumerative combinatorics, for instance in the enumeration of maps. Bousquet-Mélou and
Jehanne showed in 2006 that when these equations are of a fixed point type in F, then F is
an algebraic series. In the same paper, they proposed a systematic method for computing
annihilating polynomials of these series. Bostan, Safey El Din and the author of this package
recently designed new efficient algorithms for computing these witnesses of algebraicity. 
The package ddesolver contains the implementations of these new algorithms.

Have a look at the dedicated [tutorial paper](https://mathexp.eu/notarantonio/papers/ddesolver.pdf)!   


# **INSTALLATION**

DDE-Solver is available in the “.mla” format.

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
