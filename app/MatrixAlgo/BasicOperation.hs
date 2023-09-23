{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE DataKinds #-}
module MatrixAlgo.BasicOperation where

import Data.Mod
import GHC.Natural
import GHC.TypeNats hiding (Mod)

class (Num r, Eq r) => BasicOperation r where
    -- input  : a, b in ring r
    -- output : g,s,t,u,v such that  s * v - t * u is a unit in ring r,
    --          and s * a + t * b = g, u * a + v * b = 0 
    gcdExOp :: r -> r -> (r, r, r, r, r)

    gcdOp :: r -> r -> r

    -- input  : a in ring r
    -- output : u in r*, such that u * a in nonAssociateSet(r)
    unitOp  :: r -> r

    -- input  : a, b, n
    -- output : c such that gcd(a, b, n) = gcd(a + c * b, n)
    stabOp  :: r -> r -> r -> r

    -- input : a, b
    -- ouput : q such that a - qb in residue(a,b)
    quoOp   :: r -> r -> r


powMod :: (Integral a, Integral b) => a -> b -> a -> a
powMod x y m
    | m <= 0    = error "powModInt: non-positive modulo"
    | y <  0    = error "powModInt: negative exponent"
    | otherwise = f (x `rem` m) y 1 `mod` m
    where
        f _ 0 acc = acc
        f b e acc = f (b * b `rem` m) (e `quot` 2)
            (if odd e then b * acc `rem` m else acc)

exGcd :: Integral i => i -> i -> (i,i,i)
exGcd a 0 = (a, 1, 0)
exGcd 0 b = (b, 0, 1)
exGcd a b = (g, t - (b `div` a) * s, s)
    where
        (g, s, t) = exGcd (b `mod` a) a

ceilingLog :: Integral i => i -> i -> i
ceilingLog 0 b = error "ceilingLog: logarithm of zero"
ceilingLog n b
    | n <= 1    = 0
    | n <= b    = 1
    | otherwise = 1 + ceilingLog (n `div` b + n `rem` b) b

instance BasicOperation Integer where
    gcdExOp a b
        | a == 0 && b == 0 = (0, 1, 0, 0, 1)
        | otherwise        = (g, s, t, u, v)
        where
            (g, s, t) = exGcd a b
            u = - (b `div` g)
            v = a `div` g
    gcdOp = Prelude.gcd
    unitOp a
        | a < 0     = - 1
        | otherwise = 1

    stabOp a b n
        | f == 1    = n `div` gcdOp (powMod a (ceilingLog n 2) n) n
        | otherwise = stabOp (a `div` f) (b `div` f) (n `div` g)
            where
                f = gcdOp a b
                g = gcdOp f n

    quoOp a b = a `div` b

toI = toInteger . unMod

instance (KnownNat n) => BasicOperation (Mod n) where
    gcdExOp a b = (fromInteger g, fromInteger s, fromInteger t, fromInteger u, fromInteger v)
        where
            (g, s, t, u, v) = gcdExOp (toI a) (toI b)
    gcdOp a b = fromInteger (gcdOp (toI a) (toI b))
    unitOp a1 = fromInteger (t * (s + stabOp s u n * u))
        where
            a = toI a1
            n = toInteger $ natVal a1
            (g, s, _, u, _) = gcdExOp a n
            t = unitOp g
    stabOp a1 b1 d1 = fromInteger (stabOp a b (gcdOp d n))
        where
            (a, b, d) = (toI a1, toI b1, toI d1)
            n = toInteger $ natVal a1

    quoOp a1 b1 = fromInteger (a `div` b)
        where
            (a, b) = (toI a1, toI b1)
