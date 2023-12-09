ddesolver := module()

description "DDE-Solver, a Maple package for Discrete Differential Equations";

option package;

export annihilating_polynomial:

local Euclide, rational_interpolation, rat_interp, stickelberger2modp, elim_without_dup_modp, 
Reconstruction_Elimination_Pol_QQ1, elim_pol_1dim_main, elim_pol_1dim,
FormatOutputMSolve,GetRootsFromMSolve,ToMSolve,GetOptions,MSolveGroebner,MSolveelim,
MSolveRealRoots, pol_eval, remainder, linearize, GenerateDAC, block_elim_lex,
brut_force_elim, hgp, elim_without_dup, stickelberger2, generate:




#					README:
#
#	The present maple worksheet aims at providing an implementation of Sections 3, 4, 5, 6
#	of the 2023 ISSAC paper "Fast Algorithms for Discrete Differential Equations"
# 	by Bostan, Notarantonio and Safey El Din.
#
#	The implemented algorithms compute a witness of algebraicity for the specialized
#	solution F(t, a), associated with the variable z0.
#
#	The ddesolver package is associated with the article
#	"DDE-Solver: A Maple package for Discrete Differential Equations".
#
#	Regarding Grobner bases computations, we use the command by default Groebner[Basis].
#	It is however possible to use as in the tutorial paper the C library msolve by 
#	Berthomieu, Eder and Safey El Din, available at https://msolve.lip6.fr/ 
#	For people whishing to use msolve, it is possible to modify
#	(and use for private purposes only) the present worksheet by commenting
#	out the corresponding lines, then executing the file build.mpl





#############################################################################################
#############################################################################################
Euclide := proc(A, B, degstop, X, p)
############################################
#	-Input: A       : polynomial in Fp[X],
#   	        B       : polynomial in Fp[X],
#	  	degstop : a bound at which we stop iterating the division,
#	        X       : a variable,
#	        p       : the field characteristic.
# 	-Output: The remainder R and the cofactor of the usual extended euclidean algorithm.
############################################
	local R, Q, U, V, i, buff0, buff1;

	## Initialization + code with modp1 functions for faster computations.
	buff0 := modp1(ConvertIn(0, X), p);
	buff1 := modp1(ConvertIn(p-1, X), p);
	R := [A, B];
	U := [modp1(ConvertIn(1, X), p), buff0];
	V := [buff0, modp1(ConvertIn(1, X), p)];
	i := 2;
	Q := [buff0];
	
	do
		Q := [op(Q), modp1(Quo(R[i-1], R[i]), p)];
		R := [op(R), modp1(Rem(R[i-1], R[i]), p)];
		U := [op(U),
			modp1(Add(U[i-1],
			modp1(Multiply(buff1, modp1(Multiply(Q[i], U[i]), p)), p)), p)];
		V := [op(V),
			modp1(Add(V[i-1],
			modp1(Multiply(buff1, modp1(Multiply(Q[i], V[i]), p)), p)), p)];
		i := i + 1;
	until (modp1(Degree(R[i]), p) < degstop);
	

	return [modp1(ConvertOut(R[i], X), p),
	        modp1(ConvertOut(U[i], X), p),
	   	modp1(ConvertOut(V[i], X), p)]
	
end proc: # Euclide
#############################################################################################
#############################################################################################











#############################################################################################
#############################################################################################
rational_interpolation := proc(List_evaluation_points, List_taken_values, X, p)
############################################
#      -Input: List_evaluation_points : list of evaluation points,
#   	       List_taken_values      : list of values taken,
#	       X                      : variable to interpolate,
#	       p                      : the field characteristic.
#	   
#      -Output: A rational fraction F in the variable X interpolated
#	        from the input points and taken values in the finite field Fp with p elements.
############################################
	local m, g, R, i, n, sd;
	n := nops(List_evaluation_points);
	m := mul(X - List_evaluation_points[i], i = 1..n) mod p;
	sd := time();
	g := Interp([seq(List_evaluation_points[i], i = 1..n)],
	     			[seq(List_taken_values[i], i = 1..n)], X) mod p;
	if degree(g, X) = 0 then
	   	  return g
	end if;
	if (degree(g, X) + 1 <> nops(List_evaluation_points)) then
	
	   	 return g
	else
		  R := Euclide(modp1(ConvertIn(m, X), p), modp1(ConvertIn(g, X), p),
	     				       	       			   iquo(n, 2), X, p);
	end if;
	return normal(R[1]/R[3]) mod p
end proc: # rational_interpolation
#############################################################################################
#############################################################################################




















#############################################################################################
#############################################################################################
rat_interp := proc(List_Spe_Pols, points, principal_var, second_var, p)
############################################
#	-Input: List_Spe_Pols : A list of the specialized polynomials mod p,
#               points        : a list of points at which List_Spe_Pols are specialized,
#	        principal_var : the variable in List_Spe_Pols,
#	        second_var    : the variable at which we want to interpolate,
#	        p             : the characteristic of the field,
#	      
#	-Output: A rational interpolation of this polynomial,
#                  and we take the squarefree part of the numerator of this interpolation.
#	         We denote by pol this output.
############################################
	   local pol, deg, i, L, j, Pol_Buffer, F;
	   
	   pol := 0;
	       #lprint(List_Spe_Pols, points, p);
	   deg := degree(List_Spe_Pols[1], principal_var);

	   for i from 0 to deg do
	       L := [seq(coeff(List_Spe_Pols[j], principal_var, i),
	       	    			      j = 1..nops(List_Spe_Pols))];
	       pol := pol +
	       (rational_interpolation(points, L, second_var, p))*(principal_var)^i;
	   od;
	   ## rational interpolation done


	   pol := numer(pol) mod p;
	   ## we now select only the terms in which principal_var appears.
	   F := Factors(pol) mod p;
	   L := [seq(F[2][j][1], j = 1..nops(F[2]))];
	   pol := expand(mul(select(has, L, principal_var))) mod p;

	   return pol;
end proc: # rat_interp
#############################################################################################
#############################################################################################























#############################################################################################
#############################################################################################
block_elim_lex := proc(S, Var1, Var2, p, save_basis)
############################################
#	-Input:  S          : a list of polynomials in the variables Var1, Var2,
#	   	 Var1       : The variables which we potentially want to eliminate,
#		 Var2       : The variables which we don't want to eliminate,
#		 p          : characteristic of the field,
#		 save_basis : 0 if we want to keep Var1 in the basis, 1 else.
#	-Output:   Denote Var1 = x1, ..., xr; it returns 
#		   the Groebner basis for x1 > x2 > ... > xr > {Var2} if save_basis = 0,
#	           and with Var1 eliminated if save_basis = 1. (for {Var2} a degrevlex order
#                  over the variables Var2).
############################################
	local GB1, GB2, Buffer1, Buffer2, List1, Listelim, i;
	
	Listelim := [];
	Buffer1 := [];
	Buffer2 := [];
	
	## The aim is to eliminate variables one by one from Var1, this naturally yields a
	## disjunction of cases.
	if nops(Var1) = 1 then
	
 ## Choose one between the next two lines depending on if you are using msolve or not.
	     # -->   GB1 := MSolveelim(S, p, [op(Var1), op(Var2)], 1, 1);
	        GB1 := Groebner[Basis](S, lexdeg(Var1, Var2), characteristic = p);
	   
		if save_basis = 0 then
		   return GB1
		else
		   return remove(has, GB1, {op(Var1)})
		end if;
	else
		Listelim := [Var1[1]];
		List1 := [op({op(Var1), op(Var2)} minus {op(Listelim)})];
		
 ## Choose one between the next two lines depending on if you are using msolve or not.
	    # -->    GB1 := MSolveelim(S, p, [op(Listelim), op(List1)], 1, 1);
	       GB1 := Groebner[Basis](S, lexdeg(Listelim, List1), characteristic = p);
		
	       ## Depending of the save_basis option, call recursively the function by putting
  	       ## the proper input.
	       if save_basis = 0 then 
	              return [op(block_elim_lex(GB1, [op({op(Var1)} minus {op(Listelim)})],
	      	     			        	Var2, p, 0)), op(Buffer1)]
	       elif save_basis = 1 then
	              return block_elim_lex(remove(has, GB1, {op(Listelim)}),
	      	                [op({op(Var1)} minus {op(Listelim)})], Var2, p, 1)
               end if;
	
	end if;
end proc: # block_elim_lex
#############################################################################################
#############################################################################################














#############################################################################################
#############################################################################################
brut_force_elim := proc(S, x, y)
############################################
# 	 -Input: S a list of polynomials containing the variable x, defined over Q(y),
#	 	  such that the ideal generated by S has dimension 0 over Q(y).
#	 -Output: a generator of <S> \cap Q[x, y].
############################################

	 local elim_pol_1dim, sd;
	 
	 sd := time[real]();
	 elim_pol_1dim := collect(
	 	     op(block_elim_lex(S, [op(indets(S) minus {x, y})], [x, y], 0, 1)),
		     x, factor);
	 sd := time[real]()-sd;
	 
	 return Record(
	 'R' = op(block_elim_lex(S, [op(indets(S) minus {x, y})], [x, y], 0, 1)),
	 'time' = sd)
end proc: # brut_force_elim
#############################################################################################
#############################################################################################




















#############################################################################################
#############################################################################################
Reconstruction_Elimination_Pol_QQ1 := proc(P, M, principal_var, specialized_var)
############################################
#	-Input: P               : a bivariate polynomial in principal_var and specialized_var,
#	        M               : an integer such that P is computed mod M,
#	        Principal_Var   : the principal variable of P,
#	        Specialized_Var : the second variable of P.
#
#	-Output: A lift P1 of P over the rationals (possibly FAIL if M is not big enough),
############################################
	local i, TermesP, P1, sd;
	
	P1 := expand(P):
	TermesP := [op(P1)] mod M:
	       # TermesP contains all monomials of expand(P).
	    
	for i from 1 to nops(TermesP) do
	       P1 := P1 - TermesP[i]:
	       
	       # iratrecon is a function to recover a fraction a/b from its modular expression
	       # c mod d such that a/b = c mod d. 
	       P1 := P1 + iratrecon(coeffs(subs(principal_var = 1, specialized_var = 1,
	       	          TermesP[i]) mod M), M)*(TermesP[i]/subs(principal_var = 1,
		              specialized_var = 1, TermesP[i])):
	end do:
	
	return P1
end proc: # Reconstruction_Elimination_Pol_QQ1
#############################################################################################
#############################################################################################
















#############################################################################################
#############################################################################################
elim_pol_1dim := proc(S, principal_var, second_var)
######################################
#	Input: S: a polynomial system, expected of dimension 1,
#	       second_var: the variable on which we evaluate/interpolate,
#	       		   note that at second_var specialized, the ideal
#			   induced by S shall be of dimension 0,
#	       principal_var: the variable of which we want to compute the
#	                   specialized elimination polynomial, that is:
#			   <S> \cap \mathbb{Q}[principal_var, second_var]
#	Output: The elimination polynomial over \mathbb{Q} which generates
#		the elimination ideal of interest.
#		We denote by PolQQ1 this output.
#		Way of computing it: evaluation interpolation on second_var,
#		       		     use of multi-modular computation,
#				     computation of a degrevlex GB basis
#				     and conversion via FGLM onto a lex
#				     basis.
#####################################

	local GB, GB1, Hilb, borne_t, BufferS, sd, deg_ideal, point, ListeVar, pol, R;


 	   # Initialization of the timing (in terms of real timing: made for
	   # using 1 thread. Else it is not relevant since maple can't take into account
	   # the use of msolve.
	sd := time[real]();



	
	   # Number of points for the evaluation/interpolation
	   # on second_var. The number of point is determine by computing the precise
	   # value that we wish to compute in the end:
	   # computing a grevlex basis, and converting it into a lex basis in order
	   # to obtain the univariate elimination polynomial in second_var
	   # (only the elimination polynomial is computed, not the all GB).
	BufferS := subs(principal_var = rand(), S);
	ListeVar := [op(indets(BufferS) minus {second_var}), second_var];
	
### Choose one between the next two lines depending on if you are using msolve or not.
        #GB := MSolveGroebner(BufferS, 65521, ListeVar, 1);
	GB := Groebner[Basis](BufferS, tdeg(op(ListeVar)), characteristic = 65521);
	GB1 := Groebner[FGLM](GB, tdeg(op(ListeVar)), plex(op(ListeVar)),
			          characteristic = 65521,
			          `=`(stoppingcondition, s->has(s, indets(ListeVar)
				  minus {second_var})));
	pol := remove(has, GB1, indets(subs(principal_var = rand(), S)) minus
					{second_var});
	   # Because we do a rational interpolation:
        borne_t := 2*degree(op(pol)) + 2;


        	# Having access to the ideal's degree, denoted by 'deg_ideal'.
	BufferS := subs(second_var = rand(), S);
### Choose one between the next two lines depending on if you are using msolve or not.
        #GB := MSolveGroebner(BufferS, 65521, [op(indets(BufferS))], 1);
        GB := Groebner[Basis](BufferS, tdeg(op(indets(BufferS))), characteristic = 65521);

	Hilb := Groebner[HilbertSeries](GB, tdeg(op(indets(BufferS))), s,
	     				    characteristic = 65521):
	deg_ideal := subs(s = 1, Hilb);


	local L, p, points, Pol_list_eval, Pol_list, N, PolQQ1, PolQQ2, Pol, i;
	L := [];
	N := 1;
	PolQQ1 := 0;
	PolQQ2 := 0;
	Pol := 1;
	do
		points := []; # List of evaluation points
		Pol_list_eval := []; # List of the evaluated polynomials
			      	     # at the ptn in points.
		do
			p := nextprime(randomize() mod 2^29);
		until(evalb(p in L) = false): # We draw different primes.
		L := [op(L), p];


		   # We draw enough points for the rational
		   # evaluation/interpolation step.
		do
			point := randomize() mod p;
			if evalb(point in points) = false then
			   	   points := [op(points), point];
			end if;
		until(nops(points) = borne_t): 


	   	    # Number of points for the evaluation/interpolation
		    # on second_var. The number of point is determine by computing the precise
		    # value that we wish to compute in the end:
		    # computing a grevlex basis, and converting it into a lex basis in order
	     	    # to obtain the univariate elimination polynomial in principal_var
	     	    # (only the elimination polynomial is computed, not the all GB).
		for point in points do
### Choose one between the next two lines depending on if you are using msolve or not.
    	   ListeVar := [op(indets(S) minus {second_var, principal_var}), principal_var];
		    #	GB := MSolveGroebner(subs(second_var = point, S), p, ListeVar, 1);
			GB := Groebner[Basis](subs(second_var = point, S),
			      		tdeg(op(ListeVar)),
					characteristic = p);
			#lprint(GB);
		       	GB := Groebner[FGLM](GB, tdeg(op(ListeVar)), plex(op(ListeVar)),
			      			  characteristic = p,
			  	 `=`(stoppingcondition, m->has(m, indets(ListeVar) minus
			   			  {principal_var})));
			pol := remove(has, GB, indets(subs(second_var = point, S)) minus
					{principal_var});
			Pol_list_eval := [op(Pol_list_eval), op(pol)];
		od; # mod p, we computed enough specialization of the desired polynomial.



		   # We perform rational interpolation with all the computed polynomials mod p
		pol := rat_interp(Pol_list_eval, points, principal_var, second_var, p);

		if nops(L) > 1 then # If we draw more than one prime: 
		   	Pol := chrem([pol, Pol], [p, N]);
		else
			Pol := pol;
		end if;
		N := N*p; # We give N the value of the product of primes.


		PolQQ2 := Reconstruction_Elimination_Pol_QQ1(Pol, N, principal_var,
		       	  				    second_var);
		#	print(PolQQ2);
	#	if nops(L) > 1 then
	#	      # If two lifts corresponds successively:
	#	   if PolQQ2 = PolQQ1 and PolQQ2 <> FAIL then
	#	      	    sd := time[real]() - sd;
	#                    R := Record(
	#		    'R' = factor(numer(PolQQ2)),
	#		    'nbre_primes_used' = nops(L),
	#		    'deg_ideal' = deg_ideal,
	#		    'nbre_eval_pts_second_var' = borne_t,
	#		    'time' = sd
	#		    );
	#		    return R
	#	   end if;
	#	end if;
		if nops(L) > 1 then
		      # If two lifts correspond successively:
		   if PolQQ2 = PolQQ1 and PolQQ2 <> FAIL then
		      	    sd := time[real]() - sd;
	                    return numer(factor(numer(PolQQ2)))
		   end if;
		end if;

		PolQQ1 := PolQQ2;
	until(PolQQ1 != PolQQ2);
end proc; #elim_pol_1dim
#############################################################################################
#############################################################################################





















#############################################################################################
#############################################################################################
elim_pol_1dim_main := proc(S, principal_var, second_var)
##############################################
#	Input: S: a polynomial system, expected of dimension 1,
#	       second_var: the variable on which we evaluate/interpolate,
#	       		   note that at second_var specialized, the ideal
#			   induced by S shall be of dimension 0,
#	       principal_var: the variable of which we want to compute the
#	                   specialized elimination polynomial, that is:
#			   <S> \cap \mathbb{Q}[principal_var, second_var]
#	       k: the order of the DDE.
#	Output: Print of the function
#		elim_pol_1dim(S, principal_var, second_var, k)
##############################################
	local P;
	P := elim_pol_1dim(S, principal_var, second_var);
	print("Computed polynomial", numer(P:-R));
	print("Number of prime used", P:-nbre_primes_used);
	print("Degree of the ideal at second_var specialized", P:-deg_ideal);
	print("Number of evaluation points needed in the second variable",
		      P:-nbre_eval_pts_second_var);
	print("Timing", P:-time);
	return P
end proc; # elim_pol_1dim_main
#############################################################################################
#############################################################################################

 





















#############################################################################################
#############################################################################################
elim_without_dup_modp := proc(S, principal_var, k, p)
############################################
#	Input: S: a polynomial system over QQ in the variables x, u, z_1, ..., z_{k-1} and
#	       	    principal_var,
#	       principal_var: either z_0 or t,
#	       k: order of the DDE,
#	       p: a prime number < 2^31
#	Output:
#		a (possibly equal to 0) elimination polynomial of F(t, a)
#		   	     	   in Fp[t, z_0]
#	Assumptions: The fiber of interest shall be finite to output a non-zero
#		     bivariate polynomial.
#
#	Example: (for solving 3-constellations)
# S := [
#				# P = numer(DDE eqn)
#(u*(u-1)^2*x^3+2*u*(u-1)*x^2-u*(u*z0-z0-1)*x-u*(u*z0^2+u*z1-z0^2+z0-z1))*t-(u-1)^2*x+(u-1)^2,
#				    # diff(P, x)
#(3*u*(u-1)^2*x^2+4*u*(u-1)*x-u*(u*z0-z0-1))*t-(u-1)^2,
#			            # diff(P, u)
#((u-1)^2*x^3+2*u*(u-1)*x^3+2*(u-1)*x^2+2*x^2*u-(u*z0-z0-1)*x-x*u*z0-u*z0^2-u*z1+z0^2-z0+z1
#-u*(z0^2+z1))*t-2*(u-1)*x+2*u-2,
#			      # saturation condition 
#m*t*u*(u-1)-1];
#Output := elim_without_dup_modp(subs(t = rand(), S), z0, 2, 65521):
############################################
	local i, Gm, Gx, Gu, LeadingCoeffsx, LeadingCoeffsu, LTlowdegrees:
	local H, pol, F, ct, g_list, R, j, Ru;
	LeadingCoeffsx, LeadingCoeffsu := [], [];
	
	### The following lines allow us to determine inequations in order to
	### focus on the points in the projection (and not on the Zariski closure of it).
	Gm := block_elim_lex(S, [m], [x, u, op(indets(S) minus {m, x, u})],
	      				   	p, 1):
	Gx := block_elim_lex(Gm, [x], [u, op(indets(Gm) minus {x, u})], p, 0):
	LeadingCoeffsx := [seq(coeff(select(has, Gx, x)[i], x,
		       degree(select(has, Gx, x)[i], x)),
		       			  i=1..nops(select(has, Gx, x)))]:
	Gu := block_elim_lex(remove(has, Gx, x), [u],
	      					 [op(indets(Gx) minus {x, u})], p, 0):
	LeadingCoeffsu := [seq(coeff(select(has, Gu, u)[i], u,
		       degree(select(has, Gu, u)[i], u)),
		       			  i=1..nops(select(has, Gu, u)))]:

	Ru := op(block_elim_lex(remove(has, Gx, x),
	      [cat(z, k-1)], [u, principal_var, seq(cat(z, i), i=1..(k-2))],
	      				    	    p, 1)):
	Ru := discrim(Ru, u) mod p:
	F := Factors(Ru) mod p;
	Ru := mul(F[2][i][1], i=1..nops(F[2])) mod p:

	### All polynomials of degree in u lower than k shall be forced to zero.
	LTlowdegrees := []:
	for i to nops(Gu) do
	    if degree(Gu[i], u) < k then
	       	     LTlowdegrees := [op(LTlowdegrees), coeffs(Gu[i], u)]:
	    end if:
	od:
	

 ### Choose one between the next two lines depending on if you are using msolve or not.
 ### Possible comments can be commented out to take into account the saturation
 ### and thus the points that belong only to the projection (not to its Zariski closure).
 

		(*	   H := MSolveGroebner(
			 [op(Gu), op(LTlowdegrees),
		      	 m*Ru
			 *add(LeadingCoeffsx[i]*rand(), i=1..nops(LeadingCoeffsx))
			 *add(LeadingCoeffsu[i]*rand(), i=1..nops(LeadingCoeffsu))
			  -1],
			  p, [m, op(indets(Gu) minus {principal_var}), principal_var], 1):
		*)
			  
	         H := Groebner[Basis]([op(Gu), op(LTlowdegrees),
		      	 m*Ru*add(LeadingCoeffsx[i]*rand(), i=1..nops(LeadingCoeffsx))
			     *add(LeadingCoeffsu[i]*rand(), i=1..nops(LeadingCoeffsu))
			     -1],
			  tdeg(m, op(indets(Gu) minus {principal_var}), principal_var),
			  characteristic = p);
			  
		 pol := op(Groebner[FGLM](H, tdeg(m, op(indets(H) minus {m, principal_var}),
		     			      principal_var),
					     plex(m, op(indets(H) minus {m, principal_var}),
					      principal_var),
		     			     characteristic = p,
		     			     stoppingcondition = (s->has(s, indets(H) minus
					     				{principal_var}))));
		 F := Sqrfree(pol) mod p:
		 pol := mul(F[2][i][1], i=1..nops(F[2]));
		 if nops(indets(pol)) = 1 then
		    return [pol]
		 end if:

	Gu := block_elim_lex([op(Gu), op(LTlowdegrees)], [u],
	      					 [op(indets(Gu) minus {u})], p, 0):
						 
	## All polynomials of degree in u lower than k shall be forced to zero.
	LTlowdegrees := []:
	for i to nops(Gu) do
	    if degree(Gu[i], u) < k then
	       	     LTlowdegrees := [op(LTlowdegrees), coeffs(Gu[i], u)]:
	    end if:
	od:

	Gu := block_elim_lex([op(Gu), op(LTlowdegrees)], [u],
	      					 [op(indets(Gu) minus {u})], p, 0):
						 
	## All polynomials of degree in u lower than k shall be forced to zero.
	LTlowdegrees := []:
	for i to nops(Gu) do
	    if degree(Gu[i], u) < k then
	       	     LTlowdegrees := [op(LTlowdegrees), coeffs(Gu[i], u)]:
	    end if:
	od:
	
	## Identify the position of the polynomial in Gu with least degree >= k.
	ct := 0;
	for i to nops(Gu) do
	    if degree(Gu[i], u) >= k then ct := i; break; end if:
	od:

	
	### Investigate the conjunction of polynomial equations/inequations.
	if degree(Gu[ct], u) = k then
	   	g_list := [[0, [coeff(Gu[ct], u, degree(Gu[ct], u)),
		       [discrim(Gu[ct], u) mod p]]]];
	else
		g_list := [[0, [coeff(Gu[ct], u, degree(Gu[ct], u)),
		       		     map(s->numer(s) mod p,
				     Minors_determinants(Hermite_Quadratic_Form([Gu[ct]],
				     [u], [op(indets(Gu[ct]) minus {u})], p), k))]]]:
	end if:
### Choose one between the next two lines depending on if you are using msolve or not.
	     (* H := MSolveGroebner([op(Gu), op(LTlowdegrees), op(g_list[1][1]),
	      	 			 m*Ru*(g_list[1][2][1]*rand() +
					 add(g_list[1][2][2][j]*rand(),
					 j=1..nops(g_list[1][2][2]))
					 + add(LeadingCoeffsx[j]*rand(),
					 j=1..nops(LeadingCoeffsx))
					 + add(LeadingCoeffsu[j]*rand(),
					 j=1..nops(LeadingCoeffsu)))-1], p,
	      [m, op(indets(Gu) minus {principal_var}), principal_var], 1);
	      *)
	      
	      H := Groebner[Basis]([op(Gu), op(LTlowdegrees), op(g_list[1][1]),
	      	 			 m*Ru*(g_list[1][2][1]*rand() +
					 add(g_list[1][2][2][j]*rand(),
					 j=1..nops(g_list[1][2][2]))
					 + add(LeadingCoeffsx[j]*rand(),
					 j=1..nops(LeadingCoeffsx))
					 + add(LeadingCoeffsu[j]*rand(),
					 j=1..nops(LeadingCoeffsu)))-1],
			 tdeg(m, op(indets(Gu) minus {principal_var}), principal_var),
			 characteristic = p);

	      ## Assuming that H generates an ideal of dimension at most -1 or 0, compute
	      ## the elimination polynomial in the variable principal_var.
	    pol := op(Groebner[FGLM](H, tdeg(m, op(indets(Gu) minus {principal_var}),
	    	principal_var), plex(m, op(indets(Gu) minus {principal_var}), principal_var),
		characteristic = p,
	    	stoppingcondition = (s ->has(s, indets(H) minus	{principal_var}))));
	      ## select only the squarefree part as we are only interested in an annihilating
	      ## polynomial.
	    F := Sqrfree(pol) mod p;
	    pol := mul(F[2][j][1], j=1..nops(F[2]));
	    R := pol;#print(R);
	for i from ct+1 to nops(Gu) do
	    if degree(Gu[i], u) = k then
	        g_list := [op(g_list),
		       	  # equations:
		       [[seq(coeff(Gu[j], u, degree(Gu[j], u)), j=ct..i-1)],
		          # inequations given by the Groebner variant of the extension theorem
		       [coeff(Gu[i], u, degree(Gu[i], u)),
		       		     [discrim(Gu[i], u) mod p]]]]:
	    else
		g_list := [op(g_list),
	    	     # equations:
		       [[seq(coeff(Gu[j], u, degree(Gu[j], u)), j=ct..i-1)],
	    	     # inequations given by the Groebner variant of the extension theorem
		       [coeff(Gu[i], u, degree(Gu[i], u)),
		       		     map(s->numer(s) mod p,
				     Minors_determinants(Hermite_Quadratic_Form([Gu[i]],
				     [u], [op(indets(Gu[i]) minus {u})], p), k))]]]:
	    end if:

		     # computation of a Groebner basis respecting all the equations
		     # and saturation conditions.
 ### Choose one between the next two commands depending on whether you are using msolve or not.
	    (* H := MSolveGroebner([op(Gu), op(LTlowdegrees), op(g_list[i-ct+1][1]),
	      	 			 m*Ru*(g_list[i-ct+1][2][1]*rand() +
					 add(g_list[i-ct+1][2][2][j]*rand(),
					 j=1..nops(g_list[i-ct+1][2][2]))
					 + add(LeadingCoeffsx[j]*rand(),
					 j=1..nops(LeadingCoeffsx))
					 + add(LeadingCoeffsu[j]*rand(),
					 j=1..nops(LeadingCoeffsu)))-1], p,
	      [m, op(indets(Gu) minus {principal_var}), principal_var], 1):
	      *)
	      
	      H := Groebner[Basis]([op(Gu), op(LTlowdegrees), op(g_list[i-ct+1][1]),
	      	 			 m*Ru*(g_list[i-ct+1][2][1]*rand() +
					 add(g_list[i-ct+1][2][2][j]*rand(),
					 j=1..nops(g_list[i-ct+1][2][2]))
					 + add(LeadingCoeffsx[j]*rand(),
					 j=1..nops(LeadingCoeffsx))
					 + add(LeadingCoeffsu[j]*rand(),
					 j=1..nops(LeadingCoeffsu)))-1],
			tdeg(m, op(indets(Gu) minus {principal_var}), principal_var),
			characteristic = p):
	      
	      ## Assuming that H generates an ideal of dimension at most -1 or 0, compute
	      ## the elimination polynomial in the variable principal_var.
	    pol := op(Groebner[FGLM](H, tdeg(m, op(indets(Gu) minus {principal_var}),
	    	principal_var), plex(m, op(indets(Gu) minus {principal_var}), principal_var),
		characteristic = p,
	    	stoppingcondition = (s ->has(s, indets(H) minus	{principal_var}))));
	      ## Select only the squarefree part as we are only interested in an annihilating
	      ## polynomial.
	      
	    F := Sqrfree(pol) mod p;
	    pol := mul(F[2][j][1], j=1..nops(F[2]));
	    R := R*pol;
	    if pol = 1 then
	       break;
	    end if:
	od:
	F := Sqrfree(R) mod p;
	R := mul(F[2][i][1], i=1..nops(F[2])) mod p;
	return [R]:

end proc: # elim_without_dup_modp
#############################################################################################
#############################################################################################






















#############################################################################################
#############################################################################################
elim_without_dup := proc(S, principal_var, second_var, k)
############################################
#	Input: S: a polynomial system over QQ in the variables x, u, z_1, ..., z_{k-1} and
#	       	    principal_var,
#	       principal_var: either z_0 or t,
#	       k: order of the DDE,
#	Output:
#		
#	Assumptions: The fiber of interest shall be finite to output a non-zero
#		     bivariate polynomial.
#
#	Example: (for solving 3-constellations)
# S := [
#				# P = numer(DDE eqn)
#(u*(u-1)^2*x^3+2*u*(u-1)*x^2-u*(u*z0-z0-1)*x-u*(u*z0^2+u*z1-z0^2+z0-z1))*t-(u-1)^2*x+(u-1)^2,
#				    # diff(P, x)
#(3*u*(u-1)^2*x^2+4*u*(u-1)*x-u*(u*z0-z0-1))*t-(u-1)^2,
#				    # diff(P, u)
#((u-1)^2*x^3+2*u*(u-1)*x^3+2*(u-1)*x^2+2*x^2*u-(u*z0-z0-1)*x-x*u*z0-u*z0^2-u*z1+z0^2-z0+z1
#-u*(z0^2+z1))*t-2*(u-1)*x+2*u-2,
#				# saturation condition 
#m*t*u*(u-1)-1];
#Output := elim_without_dup(S, z0, t, 2);
############################################
	 local R, pol_spe, pol_interp, Pol_list_eval, PolQQ1, PolQQ2:
	 local i, N, p, pts, pt, sd, bound_second_var, L, Pol:
	 sd := time[real]();
	 
	 ### Compute a bound for the number of points needed for the rational interpolation
	 ### in second_var which must be done after all evaluations mod p.
	 bound_second_var := 2*degree(elim_without_dup_modp(subs(principal_var = rand(), S),
	 		     				second_var, k, 65521)[1])+2:
	#	print("upper-bound on second_var", bound_second_var);
	 L := []:
	 N := 1;
	 do
		pts := []; ### List of evaluation points mod p
		Pol_list_eval := []; ### List of specialized polynomials mod p
		do
			p := nextprime(randomize() mod 2^29);
		until(evalb(p in L) = false): ### We draw different primes.
		L := [op(L), p];
	
		### Determine bound_second_var distinct points.
		do
			pt := randomize() mod p;
			if evalb(pt in pts) = false then
			   	   pts := [op(pts), pt];
			end if;
		until(nops(pts) = bound_second_var):### We draw enough points for the rational
				                    ### evaluation/interpolation step.
		
		for pt in pts do ### collecting enough specializations mod p
		       	 pol_spe := elim_without_dup_modp(subs(second_var = pt, S),
			     					principal_var, k, p);
			 Pol_list_eval := [op(Pol_list_eval), pol_spe[1]];
			 #print(pol_spe[1]);

		od:

		### computing a rational interpolation coefficients by coefficients.
	        pol_interp := rat_interp(Pol_list_eval, pts, principal_var, second_var, p);
		#print(pol_interp, nops(L));
		### Using the Chinese Remainder Theorem
		if nops(L) > 1 then 
		   	Pol := chrem([pol_interp, Pol], [p, N]);
		else
			Pol := pol_interp;
		end if;
		N := N*p; ### We give N the value of the product of primes.

		### Trying to lift the modular image of R to the rational numbers QQ.
		PolQQ2 := Reconstruction_Elimination_Pol_QQ1(Pol, N, principal_var,
		       	  				    second_var);

					
		if nops(L) > 1 then

		      if PolQQ2 = PolQQ1 and PolQQ2 <> FAIL then
		      	    sd := time[real]() - sd;
	                    return numer(factor(numer(PolQQ2)))
		   end if;

		end if;
		PolQQ1 := PolQQ2;
	until(PolQQ1 != PolQQ2);
	return R
end proc: # elim_without_dup
#############################################################################################
#############################################################################################























#############################################################################################
#############################################################################################
pol_eval := proc(V, k, F, P, ord, const)
   local Pstock, i, n, G;
   n := nops(V);
   G := normal(F);
   if n <> k + 1 then
       print("Problem, sizes of the arguments do not match");
   end if;
   Pstock := series(P, t, ord);
   Pstock := series(subs(V[1] = F, Pstock), t, ord);
   Pstock := series(subs(V[2] = subs(u = const, normal(G)), Pstock), t, ord);
   if 2 <= k then
       for i from 3 to k + 1 do
           G := diff(G, u);
	   Pstock := series(subs(V[i] = subs(u = const, G), Pstock), t, ord);
       end do;
   end if;
   return convert(series(Pstock, t, ord), polynom);
end proc: # pol_eval
#############################################################################################
#############################################################################################














#############################################################################################
#############################################################################################
remainder := proc(C, P, F, k, n, y0, V, const)
   local Reste, dery0, i;
   Reste := C - y0*pol_eval(V, k, F, diff(P, V[1]), n, const);
   dery0 := y0;
   for i from 2 to k + 1 do
      Reste := Reste - pol_eval(V, k, F, diff(P, V[i]), n, const)*subs(u = const, dery0);
      dery0 := diff(dery0, u);
   end do;
   return Reste;
end proc: # remainder
#############################################################################################
#############################################################################################




















#############################################################################################
#############################################################################################
linearize := proc(C, P, F, k, n, V, const)
   local d, Reste, h1, h2;
   if n = 1 then
      return quo(C, subs(t = 0, diff(P, V[1])), u);
   else
      d := 1/2*n;
      h1 := linearize(C, P, F, k, d, V, const);
      Reste := remainder(C, P, F, k, n, h1, V, const);
      Reste := quo(Reste, t^d, t);
      h2 := linearize(Reste, P, F, k, d, V, const);
      return convert(series(h1 + t^d*h2, t, n), polynom);
   end if;
end proc: # linearize
#############################################################################################
#############################################################################################
















#############################################################################################
#############################################################################################
GenerateDAC := proc(P, V, k, n, c)
   local d, y0, Reste, h;
   if n = 1 then
      return solve(subs(t = 0, P), x);
   else
      d := 1/2*n;
      y0 := GenerateDAC(P, V, k, d, c);
      Reste := convert(-pol_eval(V, k, y0, P, n, c), polynom);
      Reste := quo(Reste, t^d, t);
      h := linearize(Reste, P, y0, k, d, V, c);
      return convert(series(y0 + t^d*h, t, n), polynom);
   end if;
end proc: # GenerateDAC
#############################################################################################
#############################################################################################


















#############################################################################################
#############################################################################################
#	       P: the numerator of a DDE OF A FIXED-POINT TYPE,
#	       k: the order of the DDE,
#	       c: the specialization point,
# 	       N: a power of 2.
#	Output: The first N terms of F(t, u) unique solution of P.
generate := proc(P, k, c, N)
   local i, y0, V, Fspe;
   V := [x, seq(cat(z, i), i=0..(k-1))];
   y0 := GenerateDAC(P, V, k, N, c);
   return convert(y0, polynom)
end proc: # generate
#############################################################################################
#############################################################################################





















#############################################################################################
#############################################################################################
hgp := proc(S, P, a, k)
############################################
#	Input: S: the initial polynomial system
#	       	      (e.g. S := [P, diff(P, x), diff(P, u), m*t*(u-a)*u-1] for
#		      	    P the numerator equation of a DDE),
#	       P: the numerator of a DDE (which would typically be S[1]),
#	       a: the specialization point for the series,
#	       k: the order of the DDE,
#	Output: . A polynomial R in Q[t, z0] such that R(t, F(t, a)) = 0,
#		. Intermediate data (e.g. bounds, prime associated to the bounds,
#		  	       intermediate timings)
#	Assumption: . It is possible to apply the function elim_without_dup(S, z0, t, k);
#		    . We have the equality subs(t = 0, diff(P, x)) = \pm (u-a)^k
#	Example:
#S := [(u*(u-1)^2*x^3+2*u*(u-1)*x^2-u*(u*z0-z0-1)*x-u*(u*z0^2+u*z1-z0^2+z0-z1))*t
#-(u-1)^2*x+(u-1)^2, (3*u*(u-1)^2*x^2+4*u*(u-1)*x-u*(u*z0-z0-1))*t-(u-1)^2,
#((u-1)^2*x^3+2*u*(u-1)*x^3+2*(u-1)*x^2+2*u*x^2-(u*z0-z0-1)*x
#-u*z0*x-u*z0^2-u*z1+z0^2-z0+z1-u*(z0^2+z1))*t-2*(u-1)*x+2*u-2, m*t*u*(u-1)-1];
#	hgp(S, S[1], 1, 2, 0);
############################################
	      local bt, bz0, M, Fa, p, sd, R;
	      local time_bt, time_bz0, time_expand, time_guess, time_prove;
	      sd := time[real]();
	      
	      p := nextprime(randomize() mod 2^29);

	      bt := degree(elim_without_dup_modp(subs(z0 = rand(), S), t, k, p)[1]);
	      time_bt := time[real]() - sd; sd := time[real]();
	      
              bz0 :=degree(elim_without_dup_modp(subs(t = rand(), S), z0, k, p)[1]);
	      time_bz0 := time[real]() - sd; sd := time[real]();

	      
	      Fa := subs(u=a, generate(P, k, a, 2^(floor(log(2*bt*bz0+1)/log(2))+1)));
	      Fa := series(Fa, t, degree(Fa, t));
	      time_expand := time[real]() - sd;  sd := time[real]();

	  
	      M := subs(T(t) = z0, gfun[seriestoalgeq](Fa, T(t))[1]);
	      time_guess := time[real]() - sd; sd := time[real]();

	      if (convert(series(subs(z0 = Fa, M), t, 2*bt*bz0+1), polynom) = 0) then
	      	  time_prove := time[real]() - sd;
	      	  return factor(M)
	      else
	         print("Error while executing the guess and prove strategy");
	      	 return 0;
	      end if;
end proc: # hgp
#############################################################################################
#############################################################################################
















#############################################################################################
#############################################################################################
stickelberger2modp := proc(S, principal_var, p)
############################################
#	-Input: S             : a polynomial system over QQ in the variables
#	       		          x, u, z_1, ..., z_{k-1} and principal_var,
#	       principal_var : either z_0 or t,
#	       p             : a prime number < 2^31.
# 
#	-Output: a non-zero univariate polynomial R in the variable principal_var,
#		  modular specialization (mod p) of a nonzero bivariate polynomial
#		  annihilating F(t, a).
# 
#	-Assumptions: The assumptions from Proposition 6.4 from \cite{BoNoSa23}
# 
#	-Example: (for solving 3-constellations)
# S := [
#                            # P = numer(DDE eqn):
# (u*(u-1)^2*x^3+2*u*(u-1)*x^2-u*(u*z0-z0-1)*x-u*(u*z0^2+u*z1-z0^2+z0-z1))*t-(u-1)^2*x+(u-1)^2,
#                               # diff(P, x):
# (3*u*(u-1)^2*x^2+4*u*(u-1)*x-u*(u*z0-z0-1))*t-(u-1)^2,
#                               # diff(P, u):
# ((u-1)^2*x^3+2*u*(u-1)*x^3+2*(u-1)*x^2+2*x^2*u-(u*z0-z0-1)*x-x*u*z0-u*z0^2-u*z1+z0^2-z0+z1
# -u*(z0^2+z1))*t-2*(u-1)*x+2*u-2,
#                          # saturation condition:
# m*t*u*(u-1)-1];
# Procedure to execute: stickelberger2modp(subs(t = rand(), S), z0, 65521);
############################################
	local GB, Mu, sat, D1, F2, F1, L1, L2, M, chi, Pol, i, num, den;

	#GB := MSolveelim(S, p, [m, x, u, z1, principal_var], 4, 1);
	GB := Groebner[Basis](S, lexdeg([m, x, u, z1], [principal_var]), characteristic = p);
	
	L1, L2 := Groebner[NormalSet](subs(principal_var = rand(), GB), tdeg(m, x, u, z1)):
	M := Groebner[MultiplicationMatrix](z1, L1, L2,
	     					GB, tdeg(m, x, u, z1), characteristic = p):
	chi := LinearAlgebra[CharacteristicPolynomial](M, T) mod p:

	Mu := Groebner[MultiplicationMatrix](u, L1, L2,
	     					GB, tdeg(m, x, u, z1), characteristic = p):
	sat := discrim(numer(LinearAlgebra[CharacteristicPolynomial](Mu, T)), T) mod p:
	F2 := Sqrfree(sat) mod p;
	sat := mul(F2[2][i][1], i=1..nops(F2[2])) mod p;

	num := numer(chi) mod p;
	den := denom(chi) mod p;
	
	Pol := discrim(num, T) mod p;
	F1 := Sqrfree(Pol) mod p;
	Pol := mul(F1[2][i][1], i=1..nops(F1[2])) mod p;
	Pol := op(block_elim_lex([Pol, m*sat*den-1], [m], [principal_var], p, 1));

	F1 := Sqrfree(Pol) mod p;
	Pol := mul(F1[2][i][1], i=1..nops(F1[2])) mod p;

	return Pol
end proc: # stickelberger2modp
#############################################################################################
#############################################################################################



















#############################################################################################
#############################################################################################
stickelberger2 := proc(S, principal_var, second_var)
############################################
#	-Input: S: a polynomial system over QQ in the variables x, u, z_1, ..., z_{k-1} and
#	       	    principal_var,
#	       principal_var: either z_0 or t,
#            second_var: the other variable amoung z0 and t (not principal_var)
# 
#	-Output: A non-zero polynomial R annihilating F(t, a) + some computational data.
# 
#	-Assumptions: The assumptions from Proposition 6.4 from \cite{BoNoSa23}
# 
#	-Example: (for solving 3-constellations)
# S := [
#				# P = numer(DDE eqn)
#(u*(u-1)^2*x^3+2*u*(u-1)*x^2-u*(u*z0-z0-1)*x-u*(u*z0^2+u*z1-z0^2+z0-z1))*t-(u-1)^2*x+(u-1)^2,
#				    # diff(P, x)
# (3*u*(u-1)^2*x^2+4*u*(u-1)*x-u*(u*z0-z0-1))*t-(u-1)^2,
#			            # diff(P, u)
# ((u-1)^2*x^3+2*u*(u-1)*x^3+2*(u-1)*x^2+2*x^2*u-(u*z0-z0-1)*x-x*u*z0-u*z0^2-u*z1+z0^2-z0+z1
# -u*(z0^2+z1))*t-2*(u-1)*x+2*u-2,
#			      # saturation condition 
# m*t*u*(u-1)-1];
# Procedure to execute: stickelberger2(S, z0, t);
############################################
	 local R, pol_spe, pol_interp, Pol_list_eval, PolQQ1, PolQQ2:
	 local i, product_primes, p, pts, pt, sd, bound_second_var, Lprimes, Pol:
	 
	 sd := time[real]();
	 ## Compute a bound for the number of points needed for the rational interpolation
	 ## in second_var which must be done after all evaluations mod p.
	 bound_second_var := 2*degree(stickelberger2modp(subs(principal_var = rand(), S),
	 		     				second_var, 65521))+2;
			      
	 Lprimes := []:
	 product_primes := 1;
	 
	 do
		pts := []; ### List of evaluation points mod p
		Pol_list_eval := []; ### List of specialized polynomials mod p
		do
			p := nextprime(randomize() mod 2^29);
		until(evalb(p in Lprimes) = false): ### We draw different primes.
		Lprimes := [op(Lprimes), p];
	
		## Determine bound_second_var distinct points.
		do
			pt := randomize() mod p;
			if evalb(pt in pts) = false then
			   	   pts := [op(pts), pt];
			end if;
		until(nops(pts) = bound_second_var):
				### We draw enough points for the rational
				### evaluation/interpolation step over second_var.
		
	       for pt in pts do ### collecting enough specializations mod p
		       	 pol_spe := stickelberger2modp(subs(second_var = pt, S),
			     				 principal_var, p);
			 Pol_list_eval := [op(Pol_list_eval), pol_spe];
			 #print(pol_spe);

		od:

		## Computing a rational interpolation coefficients by coefficients.
	        pol_interp := rat_interp(Pol_list_eval, pts, principal_var, second_var, p);
		
		## Using the Chinese Remainder Theorem with chrem
		if nops(Lprimes) > 1 then 
		   	Pol := chrem([pol_interp, Pol], [p, product_primes]);
		else
			Pol := pol_interp;
		end if;
		product_primes := product_primes*p;
			       	  ### We give N the value of the product of primes.

		## Trying to lift the modular image of R to the rational numbers QQ.
		PolQQ2 := Reconstruction_Elimination_Pol_QQ1(Pol, product_primes,
						principal_var,  second_var);
					
		if nops(Lprimes) > 1 then

		      ## If the lift over the rational specializes for two consecutive lap
		      ## Of the do/until loop
		      if PolQQ2 = PolQQ1 and PolQQ2 <> FAIL then
		      	    sd := time[real]() - sd;
	                    return numer(factor(numer(PolQQ2)))
		      end if;

		end if;
		PolQQ1 := PolQQ2;
	until(PolQQ1 != PolQQ2);
	return R
end proc: # stickelberger2
#############################################################################################
#############################################################################################





















#############################################################################################
#############################################################################################
#	Input: P: The polynomial associated to a numerator DDE,
#	       k: the order k of the DDE
#	       Options:
#		 algo: {identical("elimination"), identical("geometry"),
#			identical("hgp"), identical("duplication")},
#		 var : z0 or t
#	NB: Do not hesitate to modify the imput polynomial system, especially the saturation equation.
annihilating_polynomial := proc(P, k, algo:={}, var:={})
   local S, sat, i, j, algorithm, algorithm2, principal_var, second_var, Q;

   if k = 1 then 
	Q := factors(discrim(P, x));
	return factor(discrim(mul(Q[2][i][1], i=1..nops(Q[2])), u));
   end if:
   # Default implementation
   algorithm := "elimination";
   second_var := t;

   # Modifications of the default implementation
   if nops({algo, var}) = 2 then
      algorithm2, second_var := algo, var;
   end if;
   if evalb(algorithm2 in {"elimination", "hgp", "geometry", "duplication"}) = true then
      algorithm := algorithm2;
   end if;
   
   principal_var := op({z0, t} minus {second_var});

   sat := factors(subs(t = 0, diff(P, x)));
   sat := mul(sat[2][i][1], i=1..nops(sat[2]));
   
   if algorithm = "duplication" then
      S := [seq(op(subs(x = cat(x, i), u = cat(u, i), [P, diff(P, x), diff(P, u)])), i = 1..k),
      	   	 m*mul(mul((cat(u, i) - cat(u, j)), i = 1..(j - 1)), j = 1..k)
		  *mul(subs(u = cat(u, j), sat*u), j = 1..k)*t-1];
      return elim_pol_1dim(S, principal_var, second_var, k)
   end if;

   if algorithm = "hgp" then
      S := [P, diff(P, x), diff(P, u), m*sat*u*t-1];
      return hgp(S, P, op({solve(sat, u)}), k)
   end if;

   if algorithm = "geometry" then
      S := [P, diff(P, x), diff(P, u), m*sat*u*t-1];
      return stickelberger2(S, principal_var, second_var)
   end if;

   S := [P, diff(P, x), diff(P, u), m*sat*u*t-1];
   return elim_without_dup(S, principal_var, second_var, k)

end proc: # annihilating_polynomial
#############################################################################################
#############################################################################################







end module:
