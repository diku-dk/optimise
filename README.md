# Optimisation package for Futhark [![CI](https://github.com/diku-dk/optimise/workflows/CI/badge.svg)](https://github.com/diku-dk/optimise/actions) [![Documentation](https://futhark-lang.org/pkgs/github.com/diku-dk/optimise/status.svg)](https://futhark-lang.org/pkgs/github.com/diku-dk/optimise/latest/)

A collection of optimisation libraries for Futhark:

- `simplex`: A library for solving linear programming problems.

More optimisation libraries will appear in the future.

## Installation

```
$ futhark pkg add github.com/diku-dk/optimise
$ futhark pkg sync
```

## Usage (using `futhark repl`)

A linear programming problem is in *standard form* if it is phrased as follows:

$$ {\rm maximise}~\mathbf{c}^T x~
   {\rm where}~\mathbf{A} \mathbf{x} \leq \mathbf{b}, ~\mathbf{x} \geq 0
$$

In the example below, we have

- $\mathbf{A} = [[1,0],[0,2],[3,2]]$
- $\mathbf{b} = [4,12,18]$
- $\mathbf{c} = [3,5]$

```
$ futhark repl
[0]> import "lib/github.com/diku-dk/optimise/simplex"
[1]> module simplex_f64 = mk_simplex f64
[2]> simplex_f64.simplex [[1.0,0.0],[0.0,2.0],[3.0,2.0]] [4.0,12.0,18.0] [3.0,5.0]
#Ok (36.0, [2.0, 6.0])
```

## See also

* https://github.com/diku-dk/linalg
