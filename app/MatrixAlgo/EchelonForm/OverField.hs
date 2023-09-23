module MatrixAlgo.EchelonForm.OverField where

import Data.Matrix
import Data.List (findIndex)
import Data.Maybe
import Data.Ratio

type MQ = Matrix Rational
type Q = Rational
type MI = Matrix Integer
type I  = Integer

isAllZero :: (Num a, Eq a) => Matrix a -> Bool
isAllZero = all (== 0)

nonSquareIdentity :: (Num a) => Int -> Int -> Matrix a
nonSquareIdentity n m
    | n == m    = identity n
    | n <  m    = identity n <|> zero n (m - n)
    | otherwise = identity m <-> zero (n - m) m

minimalIndex :: (Num a, Eq a) => Matrix a -> Int
minimalIndex mA = fromJust (findIndex (/= 0) (toList mA))



-- <<Algorithms for Matrix Canonical Forms>>, Chapter 2, Algorithm 2.8.
-- Fixed index error and added missing d1.
--
-- A fraction free index Gauss-Jordan transform of a triple (mA, k, d0)
-- is a 5-tuple (mU, mP, r, h, d) satisfy:
-- mU is nonsingular, mP is a permutation matrix
--         (1 / d) * mU                  * (1 / d0) * mP * mA   = mR
--             ==                                 ==              ==
--     (id k <|> m? <|> zero)                        m?           m?
-- <-> (zero <|> m? <|> zero)                   <->  mM       <-> mE
-- <-> (zero <|> m? <|> id (n - k - r))         <->  m?       <-> zero
-- *block decompositions above are conformal*
--
--             mP
--             ==
--     (id k <|> zero <|> zero)
-- <-> (zero <|> m?   <|> zero)
-- <-> (zero <|> zero <|> id h)
--
-- mE = submatrix (k + 1) (n - k - r) 1 m mR
-- mE is in reduced row echelon form
-- r  == rank of rows [(k + 1) .. n] of mA
-- h  == maximal integer such that rows [(k + 1) .. (n - h)] have rank r
-- mM = submatrix (k + 1) (n - k - r) 1 m 
-- d  == d0 * det (mM)
--
-- if input is (mA, 0, d0), then
-- mU * mP * mA == d * d0 * (reducedRowEchelonForm mA)

gaussJordanFractionFreeTrans :: (Fractional a, Num a, Eq a) => Matrix a -> Int -> a -> (Matrix a, Matrix a, Int, Int, a)
gaussJordanFractionFreeTrans mA k d0
    | n <= k    = (d0mI, mI, 0, n - k, d0)
    | isAllZero (submatrix (k + 1) n 1 m mA)  = (d0mI, mI, 0, n - k, d0)
    | m == 1    = (mU, mP, r, h, d)
    | otherwise = gaussJordanFractionFreeTransHelper mA k d0
    where
        n      = nrows mA
        m      = ncols mA
        mI     = identity n
        d0mI   = scaleMatrix d0 mI
        i      = (k + 1) + minimalIndex (submatrix (k + 1) n 1 1 mA)
        mP     = permMatrix n (k + 1) i
        mPmA   = mP * mA
        r      = 1
        h      = n - i
        d      = mPmA ! (k + 1, 1)
        mU     = setElem d0 (k + 1, k + 1) (mapCol (\r x -> - (mPmA ! (r, 1))) (k + 1) (scaleMatrix d mI))

gaussJordanFractionFreeTransHelper :: (Fractional a, Num a, Eq a) => Matrix a -> Int -> a -> (Matrix a, Matrix a, Int, Int, a)
gaussJordanFractionFreeTransHelper mA k d0
    | m <= 1 = undefined
    | otherwise = (mU, mP, r, h, d)
    where
        n = nrows mA
        m = ncols mA
        mI  = identity n
        m1 = m `div` 2
        m2 = m - m1
        mA1 = submatrix 1 n 1 m1 mA
        mB  = submatrix 1 n (m1 + 1) m mA
        (mU1, mP1, r1, h1, d1) = gaussJordanFractionFreeTrans mA1 k d0
        mA2 = scaleMatrix (1 / d0) (mU1 * mP1 * mB)
        (mU2, mP2, r2, h2, d2) = gaussJordanFractionFreeTrans mA2 (k + r1) d1
        mU = scaleMatrix (1 / d1) (mU2 * (mP2 * (mU1 - scaleMatrix d1 mI) + scaleMatrix d1 mI))
        mP = mP2 * mP1
        r  = r1 + r2
        h  = min h1 h2
        d  = d2

gaussJordanTrans :: MI -> (MI , MI, Int, I)
gaussJordanTrans mA = (toI mU, toI mP, r, numerator d)
    where
        toI = fmap numerator
        (mU, mP, r, h, d) = gaussJordanFractionFreeTrans (fmap (% 1) mA) 0 1

reducedRowEchelonForm :: MI -> (MI , Int)
reducedRowEchelonForm mA = (fmap (`div` d) (mU * mP * mA), r)
    where
        (mU, mP, r, d) = gaussJordanTrans mA


-- input  : A matrix in echelon form,
--          rank of this matrix
-- output : A list, with i-th element be the column of pivot in i-th row.
rankProfile :: MI -> Int -> [Int]
rankProfile mE r =  map (\i -> 1 + minimalIndex (submatrix i i 1 (ncols mE) mE)) [1..r]


-- input : m (ncols), rank, and rankProfile of an n * m matrix
-- output : a matrix mQ, which satisfy
--          mE * mQ
--             ==
--         (id r <|> m?)
--     <-> (m?   <|> m?)
-- where mE is reducedRowEchelonForm of the n * m matrix
principalSubMatrixTrans :: Int -> Int -> [Int] -> MI
principalSubMatrixTrans m r js = foldr (uncurry switchCols) mI qs
    where
        mI = identity m
        qs = zip [1..r] js

-- li = [15, 10, 10, 93, 10, 16, 6, 10, 10, 10, 7, 10, 12, 11, 11, 137, 11, 16, 9, 11, 11, 11, 14, 11, 11, 7, 7, 94, 7, 6, 8, 7, 7, 7, 9, 7, 4, 0, 0, 14, 0, 0, -1, 0, 0, 0, -2, 0, 8, 9, 9, 56, 9, 16, 5, 9, 9, 9, 10, 9, 15, 10, 10, 196, 10, 12, 16, 10, 10, 10, 17, 10, 3, 5, 5, 107, 5, 3, 9, 5, 5, 5, 8, 5]
-- mAA = fromList 7 12 li
-- mA  = transpose mAA
-- (mE, r) = reducedRowEchelonForm mA
-- js = rankProfile mE r


