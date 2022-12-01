{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module ChainComplex where

import Data.Functor.Product

newtype FormalSum a = FormalSum {sum :: [(Int,a)]}
    deriving (Eq, Functor)

class Eq (Basis name) => ChainComplex name where
    type Basis name

    isBasis :: name -> Basis name -> Bool

    degree  :: name -> Basis name -> Int
    diff    :: name -> Basis name -> FormalSum (Basis name)

data CCMorphism a b = CCMorphism {morDegree :: Int, f :: Basis a -> FormalSum (Basis b)}

