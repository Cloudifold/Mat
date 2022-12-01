{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}


module SSet where

import Data.Kind (Type)
import Data.List
import GHC.Base (ap)
import qualified Control.Category as Cat


data FormalDegen a = Degen [Int] a
    deriving Eq
    deriving Functor

instance Show a => Show (FormalDegen a) where
    show (Degen [] a) = show a
    show (Degen (j : xs) a) = "s_" ++ show j ++ show (Degen xs a)

instance Applicative FormalDegen where
    pure = Degen []
    (<*>) = ap

instance Monad FormalDegen where
    (Degen [] a) >>= f = f a
    (Degen xs a) >>= f = case f a of
        Degen ns b -> Degen (degenListNormalise (xs ++ ns)) b

degenListAppendOne :: [Int] -> Int -> [Int]
degenListAppendOne (j : xs) i | i <= j = j + 1 : degenListAppendOne xs i
degenListAppendOne xs i = i : xs

degenListNormalise :: [Int] -> [Int]
degenListNormalise [] = []
degenListNormalise [j] = [j]
degenListNormalise (j : xs) = degenListAppendOne (degenListNormalise xs) j

degenAppendOne :: FormalDegen a -> Int -> FormalDegen a
degenAppendOne (Degen xs a) i = Degen (degenListAppendOne xs i) a

degenNormalise :: FormalDegen a -> FormalDegen a
degenNormalise (Degen xs a) = Degen (degenListNormalise xs) a

degenCount :: FormalDegen a -> Int
degenCount (Degen xs a) = length xs



type Simplex name = FormalDegen (GeomSimplex name)


class Eq (GeomSimplex name) => SSet name where
    type GeomSimplex name

    isGeomSimplex :: name -> GeomSimplex name -> Bool
    isGeomSimplex _ _ = True

    dimGeomSimplex :: name -> GeomSimplex name -> Int
    -- dim|_GeomSimplex

    faceGeomSimplex :: name -> GeomSimplex name -> Int -> Simplex name

isSimplex :: SSet name => name -> Simplex name -> Bool
isSimplex n (Degen xs s) = isGeomSimplex n s

facesOfGeomSimplex :: SSet name => name -> GeomSimplex name -> [Simplex name]
facesOfGeomSimplex a s =
    let d = dimGeomSimplex a s
        in if d == 0 then [] else fmap (faceGeomSimplex a s) [0 .. d]

dim :: SSet name => name -> Simplex name -> Int
dim n (Degen xs s) = dimGeomSimplex n s + length xs


face :: SSet name => name -> Simplex name -> Int -> Simplex name
face n (Degen [] s) i = faceGeomSimplex n s i
face n (Degen (j : xs) s) i
    | i < j     =  degenAppendOne (face n (Degen xs s) i) (j - 1)
    | i > j + 1 = degenAppendOne (face n (Degen xs s) (i - 1)) j
    | otherwise = Degen xs s

hasFace :: SSet a => a -> GeomSimplex a -> GeomSimplex a -> Bool
hasFace a t s = Degen [] s `elem` facesOfGeomSimplex a t
-- t hasFace s


class SSet name => FiniteSSet name where
    geomSimplices :: name -> Int -> [GeomSimplex name]
    -- geomSimplices name dimension

class SSet name => Pointed name where
    basepoint :: name -> GeomSimplex name

basepointSimplex :: (Pointed name) => name -> Simplex name
basepointSimplex n = Degen [] (basepoint n)



{- (Degen i s, Degen j t)
    | i == j = degen (prodNormalise (s, t)) i
    | i > j =
        let p = prodNormalise (s, Degen j t)
            in fmap (\(s', t') -> (Degen (i - degenCount p) s', t')) p
    | i < j =
            let p = prodNormalise (Degen i s, t)
                in fmap (\(s', t') -> (s', Degen (j - degenCount p) t')) p
normaliseProd s = NonDegen s
-}



newtype SMorphism a b = SMorphism {f :: GeomSimplex a -> FormalDegen (GeomSimplex b)}

morphismOnSimplex :: SMorphism a b -> FormalDegen (GeomSimplex a) -> FormalDegen (GeomSimplex b)
morphismOnSimplex (SMorphism f) s = s >>= f





