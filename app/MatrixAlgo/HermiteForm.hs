module MatrixAlgo.HermiteForm where

import MatrixAlgo.BasicOperation
import MatrixAlgo.EchelonForm.OverPIR
import Data.Matrix hiding (trace)

import Debug.Trace
import Data.Maybe (fromJust)
import Data.List (find, findIndex)

-- Corollary 3.11.
-- Assuming mA is upper triangular, phi i j (a + (phi i j a) * mA ! (j, j)) == 0
indexReductionTransform :: (BasicOperation r, Show r) => (Int -> Int -> r -> r) -> Matrix r -> Int  -> Matrix r
indexReductionTransform phi mA 0 = identity (nrows mA)
indexReductionTransform phi mA 1 = mapPos (\(a, b) c -> if b == nrows mA && a < nrows mA
                                                        then phi a b (mA ! (a, b))
                                                        else c) (identity (nrows mA))
indexReductionTransform phi mA k = mU
    where
        n = nrows mA
        k2 = k `div` 2
        k1 = k - k2
        mA1 = submatrix 1 (n - k2) 1 (n - k2) mA
        mV1 = indexReductionTransform phi mA1 k1
        mU1 = mapPos (\(a, b) c -> if a <= (n - k2) && b <= (n - k2)
                                    then mV1 ! (a, b)
                                    else c) (identity n)
        mU2 = indexReductionTransform phi (mU1 * mA) k2
        mU = mU2 * mU1

hermiteReductionTransform :: (BasicOperation r, Show r) =>  Matrix r -> Matrix r
hermiteReductionTransform mA = mU
    where
        n = nrows mA
        phi i j a = - (quoOp a (mA ! (j, j)))
        mU = indexReductionTransform phi mA n

-- Corollary 3.16
hermiteForm :: (BasicOperation r, Show r) => Matrix r -> Matrix r
hermiteForm mA = mU
    where
        n = nrows mA
        m = ncols mA
        mU1 = rowEchelonFormBottomToTop mA
        mT1 = mU1 * mA

        scanrow rowstart rowend colstart colend mat = case
            find
                (\ x -> all (== 0) (submatrix x x colstart colend mat))
                [rowstart .. rowend]
            of
            Nothing -> rowend + 1
            Just z -> z

        r = scanrow 1 n 1 m mT1 - 1
        minimalIndex mX = fromJust (findIndex (/= 0) (toList mX))
        rankProfile mE r =  map (\i -> 1 + minimalIndex (submatrix i i 1 (ncols mE) mE)) [1..r]
        js = rankProfile mT1 r
        mP = foldr (uncurry switchCols) (identity m) (zip [1..r] js)
        mE1 = mT1 * mP
        mE =  submatrix 1 r 1 r mE1
        mU2 = diagonalList r 0 (map (\i -> unitOp (mE ! (i, i))) [1..r])
        mU3 = hermiteReductionTransform (mU2 * mE)
        
        mU32 = mapPos (\(a, b) c -> if a <= r && b <= r
                                    then (mU3 * mU2) ! (a, b)
                                    else c) (identity n)

        mU = mU32 * mU1