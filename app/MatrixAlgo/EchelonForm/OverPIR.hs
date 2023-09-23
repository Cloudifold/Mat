module MatrixAlgo.EchelonForm.OverPIR where

import MatrixAlgo.BasicOperation
import Data.Matrix hiding (trace)

import Debug.Trace
import Data.Maybe (fromJust)
import Data.List (find, findIndex)

debug = flip trace

-- Corollary 3.2.
-- input : (t + k) * m matrix with last k rows in row echelon form
--          t <= m, k <= m, m is power of 2
-- output : (t * k) * (t * k) matrix U such that U * A is in row echelon form
echelonFormDivideByColsRecur :: (BasicOperation r, Show r) => Matrix r -> Int -> Int -> (Matrix r, Matrix r)
echelonFormDivideByColsRecur mA t k
    | t == 0                             = (identity k, mA)
    | t == 1 && k == 0                   = (identity t, mA)
    | ncols mA == 1 && t == 1 && k == 1  = let (g, s, t, u, v) = gcdExOp (mA ! (1, 1)) (mA ! (2, 1))
                                            in (fromLists [[s, t], [u, v]], fromLists [[g], [0]])
    | ncols mA == 1                      = (identity (nrows mA), mA)
    | otherwise = echelonFormDivideByColsFullRecurCase mA t k


-- 1 <= t, 0 <= k
echelonFormDivideByColsFullRecurCase :: (BasicOperation r, Show r)
                                        => Matrix r -> Int -> Int -> (Matrix r, Matrix r)
echelonFormDivideByColsFullRecurCase mA t k = (mU, mT)
    where
        scanrow rowstart rowend colstart colend mat = case
            find
                (\ x -> all (== 0) (submatrix x x colstart colend mat))
                [rowstart .. rowend]
            of
            Nothing -> rowend + 1
            Just z -> z

        augmentTo matrix1 rowstart rowend colstart colend (i,j) a =
            if (rowstart <= i && i <= rowend) && (colstart <= j && j <= colend)
                then matrix1 ! (i - rowstart + 1, j - colstart + 1)
                else a

        mIn = identity n
        n = t + k

        m = ncols mA

        -- partition of cols : m1
        m1 = m `div` 2
        m2 = m - m1

        -- partition of rows : t1, t, k1

        t1 = t `div` 2
        k1 = scanrow (t + 1) n 1 m1 mA - 1



        (mSubU1, _) = echelonFormDivideByColsRecur (submatrix (t1 + 1) k1  1 m1 mA) (t - t1) (k1 - t)
        mU1    = mapPos (augmentTo mSubU1 (t1 + 1) k1 (t1 + 1) k1) mIn
        mA1 = mU1 * mA

        -- use a new partition of rows : t1, y2, k2
        y2 = scanrow (t1 + 1) n 1 m1 mA1 - 1
        k2 = k1

        (mSubU2, _) = echelonFormDivideByColsRecur (submatrix (y2 + 1) n (m1 + 1) m mA1) (k2 - y2) (n - k2)
        mU2 = mapPos (augmentTo mSubU2 (y2 + 1) n (y2 + 1) n) mIn
        mA2 = mU2 * mA1

        -- use a new partition of rows : t1, y2
        (mSubU3, _) = echelonFormDivideByColsRecur (submatrix 1 y2 1 m1 mA2) t1 (y2 - t1)
        mU3 = mapPos (augmentTo mSubU3 1 y2 1 y2) mIn
        mA3 = mU3 * mA2

        -- use a new partition of rows : y4, y2, k4
        y4 = scanrow 1 n 1 m1 mA3 - 1
        k4 = scanrow (y2 + 1) n (m1 + 1) m mA3 - 1

        (mSubU4, _) = echelonFormDivideByColsRecur (submatrix (y4 + 1) k4  (m1 + 1) m mA3) (y2 - y4) (k4 - y2)
        mU4    = mapPos (augmentTo mSubU4 (y4 + 1) k4 (y4 + 1) k4) mIn
        mU  = mU4 * mU3 * mU2 * mU1
        mT  = mU4 * mA3

echelonFormDivideByCols :: (BasicOperation r, Show r) => Matrix r -> Int -> Int -> (Matrix r, Matrix r)
echelonFormDivideByCols mB t k = (mU, mT)
    where
        n     = nrows mB
        m     = ncols mB
        sizeM = 2 ^ ceilingLog m 2
        z     = sizeM - m
        mA    = mB <|> zero n z
        (mU, mQ) = echelonFormDivideByColsRecur mA t k
        mT = submatrix 1 n 1 m mQ


-- Corollary 3.4.
echelonFormLeftToRightRecur :: (BasicOperation r, Show r) => Matrix r -> Int -> (Matrix r, Int)
echelonFormLeftToRightRecur mUAi z = (mUi, z')
    where
        -- assume m == n
        n = nrows mUAi
        m = ncols mUAi
        mB = submatrix (1 + n - z) n 1 m mUAi
        (mV, mVB) = echelonFormDivideByCols mB z 0

        mUi = mapPos (\(a, b) c -> if a > (n - z) && b > (n - z)
                                    then mV ! (a - n + z, b - n + z)
                                    else c) (identity n)

        nTailingZeroRows mat nrow ncol = case
            find (\ x -> any (/= 0) (submatrix x x 1 ncol mat))
                (reverse [1 .. nrow])
            of
                Just k -> nrow - k
                Nothing -> nrow

        z' = nTailingZeroRows mVB z m

rowEchelonFormLeftToRight :: (BasicOperation r, Show r) => Matrix r -> Matrix r
rowEchelonFormLeftToRight mA' = mU
    where
        n = nrows mA'
        m' = ncols mA'
        pad = n - (m' `rem` n)
        mA = if pad == 0
                then mA'
                else mA' <|> zero n pad

        m = ncols mA
        l = m `div` n

        listAi = map (\i -> submatrix 1 n (1 + n * i) (n * i + n) mA) [0..(l - 1)]


        func (mU,z) mAi = let (mUi, z') = echelonFormLeftToRightRecur (mU * mAi) z in (mUi * mU, z')

        (mU,z) = foldl func (identity n, n) listAi

-- Corollary 3.6.
echelonFormBottomToTopRecur :: (BasicOperation r, Show r) => Matrix r -> Int -> Int -> (Matrix r, Int, Int)
echelonFormBottomToTopRecur mA z r
    | z == n    = (identity n, z, r)
    | otherwise = (mU * mUi, z', r')
    where
        n = nrows mA
        m = ncols mA

        d = min (max r 1) (n - z)
        z' = z + d
        mB = submatrix (n - z' + 1) n 1 m mA
        mV = rowEchelonFormLeftToRight mB
        mVB = mV * mB

        mUi = mapPos (\(a, b) c -> if a > (n - z') && b > (n - z')
                                    then mV ! (a - n + z', b - n + z')
                                    else c) (identity n)

        nonzeroRows mat nrow ncol = length (filter id
                        (map (\ x -> any (/= 0) (submatrix x x 1 ncol mat))
                            [1 .. nrow]))

        r' = nonzeroRows mVB z' m

        (mU, zz, rr) = echelonFormBottomToTopRecur (mUi * mA) z' r'

rowEchelonFormBottomToTop :: (BasicOperation r, Show r) => Matrix r -> Matrix r
rowEchelonFormBottomToTop mA = mU
    where
        (mU, _, _) = echelonFormBottomToTopRecur mA 1 1



