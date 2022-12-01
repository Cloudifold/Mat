{-# LANGUAGE TypeFamilies #-}
module SSet.StandardNSimplex where


import SSet

newtype NSimplex = NSimplex {simplexDimension :: Int}

instance Show NSimplex where
    show (NSimplex n) = "Δ^" ++ show n

newtype NSimplexSimplices = NSimplexSimplices [Int]
    deriving (Eq, Ord, Show)

isOrdered :: (Ord a) => [a] -> Bool
isOrdered [] = True
isOrdered [x] = True
isOrdered (x : y : xs) = (x <= y) && isOrdered (y : xs)

deleteAt :: Int -> [a] -> [a]
deleteAt _ [] = []
deleteAt 0 (x : xs) = xs
deleteAt i (x : xs) = x : deleteAt (i - 1) xs

instance SSet NSimplex where
    type GeomSimplex NSimplex = [Int]

    isGeomSimplex (NSimplex n) vs = (length vs - 1) <= n && isOrdered vs

    dimGeomSimplex (NSimplex n) vs = length vs - 1

    faceGeomSimplex (NSimplex n) vs i = Degen [] (deleteAt i vs)

choose :: Int -> [a] -> [[a]]
choose 0 _ = [[]]
choose i [] = []
choose i (x : xs) = fmap (x :) (choose (i -1) xs) ++ choose i xs

instance FiniteSSet NSimplex where
    geomSimplices (NSimplex n) i = choose (i + 1) [0 .. n]