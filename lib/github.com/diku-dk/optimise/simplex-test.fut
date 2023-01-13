
import "simplex"

module S = mk_simplex f64

-- test 0
-- ==
-- input { [[1.0,0.0],[0.0,2.0],[3.0,2.0]] [4.0,12.0,18.0] [3.0,5.0] }
-- output { 36.0 }

-- test 1 (Introduction to Algorithms, Second Edition, p. 791)
-- ==
-- input { [[1f64,1f64,3f64],[2f64,2f64,5f64],[4f64,1f64,2f64]] [30f64,24f64,36f64] [3f64,1f64,2f64] }
-- output { 28f64 }

-- test 2
-- ==
-- input { [[93.44461,71.88667,74.54679],[79.33898,75.88172,79.05727],[92.16154,57.39885,65.48771]]
--         [31.72006,23.99789,47.660904]
--         [76.21876,69.50258,46.55636] }
-- output { 23.054108 }

-- test 3
-- ==
-- input { [[2.0, 3.0, 1.0], [4.0, 1.0, 2.0], [3.0,4.0,2.0]]
--         [5.0, 11.0, 8.0]
--         [5.0,4.0,3.0] }
-- output { 13.0 }

-- test 4
-- ==
-- input { [[91.0,70.0],[60.0,74.0]]
--         [43.0,72.0]
--         [31.0,11.0] }
-- output { 14.648352 }

local open S

let main [m] [n] (A:[m][n]t) (b:[m]t) (c:[n]t) =
  let x = match simplex A b c case #Ok p -> p.0 case #Unbounded -> -1
  let y = match simplex_orig A b c case #Ok p -> p.0 case #Unbounded -> -1
  in (x+y) / 2.0
