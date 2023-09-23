module Main where
import MatrixAlgo.EchelonForm.OverField
import Data.Matrix


li = [15, 10, 10, 93, 10, 16, 6, 10, 10, 10, 7, 10, 12, 11, 11, 137, 11, 16, 9, 11, 11, 11, 14, 11, 11, 7, 7, 94, 7, 6, 8, 7, 7, 7, 9, 7, 4, 0, 0, 14, 0, 0, -1, 0, 0, 0, -2, 0, 8, 9, 9, 56, 9, 16, 5, 9, 9, 9, 10, 9, 15, 10, 10, 196, 10, 12, 16, 10, 10, 10, 17, 10, 3, 5, 5, 107, 5, 3, 9, 5, 5, 5, 8, 5]
mAA = fromList 7 12 li
mA  = transpose mAA
(mE, r) = reducedRowEchelonForm mA
js = rankProfile mE r

main :: IO ()
main = print (mE,r)
