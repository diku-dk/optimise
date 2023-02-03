-- | Library for solving linear programming problems using the
-- (revised) Simplex method.
--
-- The result of applying `mk_simplex`@term is a module implementing
-- both Dantzig's original Simplex method and the revised method,
-- which is more efficient than the original method.
--
-- The implementations are based on the matrix formulation described
-- in [1, Chapter 5].
--
-- [1] Frederick S. Hillier, Gerald J. Lieberman. Introduction to
-- Operations Research. 9th edition. Mc Graw Hill. 2010.
local module type simplex = {
  -- | The scalar type.
  type t

  -- | The result type.
  type result[n] = #Ok (t,[n]t) | #Unbounded

  -- | Returns the optimal solution to the objective function (a
  -- vector) and the optimal objective value, given `A`, a
  -- representation of the constraint coefficients (m rows, n
  -- columns), `b`, a representation of the constraint values (length
  -- m), and `c`, a representation of the objective coefficients
  -- (length n). The `simplex` function implements the revised simplex
  -- method, which iteratively and effciently modifies the inverse
  -- basis, without directly computing the matrix inverse.
  val simplex      [n][m] : (A:[m][n]t) -> (b:[m]t) -> (c:[n]t) -> result[n]

  -- | Returns the optimal solution to the objective function (a
  -- vector) and the optimal objective value, given `A`, a
  -- representation of the constraint coefficients (m rows, n
  -- columns), `b`, a representation of the constraint values (length
  -- m), and `c`, a representation of the objective coefficients
  -- (length n). The `simplex_orig` function implements the original
  -- simplex method, which involves computing the inverse of the basis
  -- matrix at every iteration.
  val simplex_orig [n][m] : (A:[m][n]t) -> (b:[m]t) -> (c:[n]t) -> result[n]
}

import "../linalg/linalg"

-- | Given some numeric module, produce a simplex module.
module mk_simplex (T:field) : simplex with t = T.t = {

module L = mk_linalg T

type t = T.t

type var = { name: i64 }

---------------------
-- Some utilities
---------------------

let novar: var = { name = -1 }

-- [eye n] returns the identity matrix of size n x n.

let zero : t = T.i64 0
let inf : t = T.(i64 1 / zero)

let bool b : t = T.i64(i64.bool b)

let eye (n:i64) : *[n][n]t =
  tabulate_2d n n (\i j -> bool(i == j))

let vecmatmul [n][m] (v:[m]t) (M:[m][n]t) : *[n]t =
  map (\(c:[m]t) -> L.dotprod v c) (transpose M)

let matvecmul [n][m] (M:[n][m]t) (v:[m]t) : *[n]t =
  L.matvecmul_row M v

let indexof [n] (xs: [n]var) (x: var): i64 =
    xs |> map2 (\i x' -> if x == x' then i else n) (iota n) |> reduce i64.min n

-- [init A b c] returns (xB,cB,B,bs,nbs), where xB is the initial
-- vector of basic variable values, cB is the initial vector of basic
-- variable function coefficients, B is the initial basis matrix, bs
-- is the initial vector of basic variables (ids) and nbs is the
-- initial vector of nonbasic variables (ids). Variable ids are
-- enumerated from 0, starting with the original variables, then the
-- slack variables.

let init [m] (n:i64) (b:[m]t) : (*[m]t,*[m]t,*[m][m]t,*[m]var,*[n]var) =
  (copy b,replicate m zero,eye m,map (\name -> { name = name + n }) (iota m),map (\name -> { name }) (iota n))

let entering [n][m] (A:[m][n]t) (Binv:[m][m]t) (cB:[m]t) (c:[n]t) (nbs:[n]var) : var =
  let cB_Binv : *[m]t = vecmatmul cB Binv
  let orig : [n]t = map2 (T.-) (vecmatmul cB_Binv A) c
  let slack : [m]t = cB_Binv
  let res: [n](var, t) = map (\v ->
               let c = if v.name < n then orig[v.name] else slack[v.name - n]
               in if T.(c < zero) then (v, c) else (novar, zero)
		) nbs
  in reduce (\x y -> if T.(x.1 < y.1) then x else y)
	    (novar,zero) res
     |> (.0)

let leaving [n][m] (A:[m][n]t) (xB:[m]t) (e:i64) : i64 =
  map3 (\i a b -> T.(if i64 0 < a[e] then (i, b/a[e]) else (-1, inf)))
       (iota m) A xB
  |> reduce (\x y -> if T.(x.1 < y.1) then x else y)
	    (-1,inf)
  |> (.0)

let continue (_count:i64) (_bound:i64) : bool = true
  -- count < bound

type status = #FlagOk | #FlagUnbounded

-- [simplexO A b c] returns intermediate data necessary for computing
-- the optimal solution to the objective function (a vector) and the
-- optimal objective value, given A, a representation of the
-- constraint coefficients (m rows, n columns), b, a representation of
-- the constraint values (length m), and c, a representation of the
-- objective coefficients (length n). The simplexO function implements
-- the original simplex method, which involves computing the inverse
-- of the basis matrix at every iteration.

let simplexO [n][m] (A:[m][n]t) (b:[m]t) (c:[n]t) =
  let bound = n+m
  let (xB:[m]t,cB:[m]t,B:[m][m]t,bs:[m]var,nbs:[n]var) = init n b
  let I = eye m
  let k = entering A B cB c nbs
  in loop (B,xB,cB,k,bs,nbs,count,status:status) = (B,xB,cB,k,bs,nbs,0,#FlagOk)
     while status == #FlagOk && k != novar && continue count bound
     do let ki = indexof nbs k
        let ri = leaving A xB ki
        let r = if ri == -1 then novar else bs[ri]
	in if r == novar
	   then (B,xB,cB,k,bs,nbs,count+1,#FlagUnbounded)
	   else let B' =
		  let col = if k.name < n then A[:,k.name] else I[:,k.name-n]
		  in B with [:,ri] = col
  		let Binv' = L.inv B'
        let bs' = bs with [ri] = k
        let nbs' = nbs with [ki] = r
  		let xB' = matvecmul Binv' b
		let cB' = cB with [ri] = c[ki]
        let k' = entering A Binv' cB' c nbs'
		in (B',xB',cB',k',bs',nbs',count+1,#FlagOk)

-- [simplexR A b c] returns intermediate data necessary for computing
-- the optimal solution to the objective function (a vector) and the
-- optimal objective value, given A, a representation of the
-- constraint coefficients (m rows, n columns), b, a representation of
-- the constraint values (length m), and c, a representation of the
-- objective coefficients (length n). The simplexR function implements
-- the revised simplex method, which iteratively and effciently
-- modifies the inverse basis, without directly computing the matrix
-- inverse.

let simplexR [n][m] (A:[m][n]t) (b:[m]t) (c:[n]t) =
  let bound = n+m
  let (xB:[m]t,cB:[m]t,Binv:[m][m]t,bs:[m]var,nbs:[n]var) = init n b
  let k = entering A Binv cB c nbs
  in loop (Binv:[m][m]t,xB,cB,k,bs,nbs,count,status:status) = (Binv,xB,cB,k,bs,nbs,0,#FlagOk)
     while status == #FlagOk && k != novar && continue count bound
     do let ki = indexof nbs k
        let ri = leaving A xB ki
        let r = if ri == -1 then novar else bs[ri]
	in if r == novar
	   then (Binv,xB,cB,k,bs,nbs,count+1,#FlagUnbounded)
	   else let a_k' = if k.name < n
			   then matvecmul Binv A[:,k.name]
			   else Binv[:,k.name-n]
		let a_rk' = a_k'[ri]
		let Binv_r = Binv[ri]
		let Binv' = tabulate_2d m m
					T.(\i j -> if i == ri
						   then Binv_r[j] / a_rk'
						   else Binv[i,j] -
							a_k'[i] * Binv_r[j] / a_rk'
					  )
        let bs' = bs with [ri] = k
        let nbs' = nbs with [ki] = r
  		let xB' = matvecmul Binv' b
		let cB' = cB with [ri] = c[ki]
		let k' = entering A Binv' cB' c nbs'
		in (Binv',xB',cB',k',bs',nbs',count+1,#FlagOk)

type result [n] = #Ok (t,[n]t) | #Unbounded

let runner [m][n]
	   (looper: [m][n]t -> [m]t -> [n]t ->
	      ([m][m]t,[m]t,[m]t,var,[m]var,[n]var,i64,status))
           (A:[m][n]t) (b:[m]t) (c:[n]t) : result[n] =
  let (_B,xB,cB,_,bs,_nbs,_count,status) = looper A b c
  in if status == #FlagUnbounded then #Unbounded
     else let h = n+m
	  let sol = scatter (replicate h zero) (map (.name) bs) xB
	  let value = L.dotprod cB xB
	  in #Ok (value,sol[:n])

let simplex [n][m] (A:[m][n]t) (b:[m]t) (c:[n]t) : result[n] =
  runner simplexR A b c

let simplex_orig [n][m] (A:[m][n]t) (b:[m]t) (c:[n]t) : result[n] =
  runner simplexO A b c
}
