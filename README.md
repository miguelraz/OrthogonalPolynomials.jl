# OrthogonalPolynomials

[![Build Status](https://travis-ci.org/miguelraz/OrthogonalPolynomials.jl.svg?branch=master)](https://travis-ci.org/miguelraz/OrthogonalPolynomials.jl) [![Coverage Status](https://coveralls.io/repos/miguelraz/OrthogonalPolynomials.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/miguelraz/OrthogonalPolynomials.jl?branch=master) [![codecov.io](http://codecov.io/github/miguelraz/OrthogonalPolynomials.jl/coverage.svg?branch=master)](http://codecov.io/github/miguelraz/OrthogonalPolynomials.jl?branch=master) [![Build status](https://ci.appveyor.com/api/projects/status/tgs1489v06dlk3u1/branch/master?svg=true)](https://ci.appveyor.com/project/miguelraz/orthogonalpolynomials-jl/branch/master)

## Introduction

Hello and welcome to OrthogonalPolynomials v0.0.1!

This is an open source project with a few simple goals:

1. Provide an easy interface for generating arbitrary order generalized Orthogonal Polynomials and efficient evaluation and using them in Julia while being robustly tested.
2. Show the process of producing a working software that is shared for others for academic and applied purposes, whilst building a compendium of handy software practices and principles.

In essence, you rarely get to see anything but the highlights of other's work. This project is an ode to all the pulled hairs and silly mistakes and drowning frustration to getting. this. darn. code. to run properly.

This goal will be met by [video tutorials](https://www.youtube.com/channel/UC840v4b_71e78fmPHiCPQVg) and [livestreaming](https://www.twitch.tv/brainrpg) where you can join in on the action if you so desire, and observe the entirety of the software design process.

## Code and Roadmap

This is how we well try to meet goal 1 above, in about ~100 lines of Julia code.
We will use some:

1. Metaproramming in Julia á la [@horner's macro](https://youtu.be/dK3zRXhrFZY) and some generated functions

2. A very convenient formula that we found in [Abramowitz and Stegun](http://people.math.sfu.ca/~cbm/aands/page_789.htm) that lends itself to multiple dispatch and generating the coefficients we want for our `@horner` formula at compile time.

3. Multiple dispatch design.

In short, the user dream is this:

```
julia> Pkg.add("OrthogonalPolynomials");

julia> using OrthogonalPolynomials, BenchmarkTools;

julia> L₁₂(x) = a(Laguerre(12));

julia> @btime L₁₂(.5)
  1.449 ns (0 allocations: 0 bytes)
-0.23164963886602852
```

## Motivation and related work

Orthogonal Polynomials are polynomials that appear almost everywhere in applied mathematics and sciences due to their centuries-old connections established to differential equations, group theory, numerical analysis, etc.

People find them very interesting for many different reasons, and you can find them in some well established software like `MATLAB` or `Mathematica` ( if you want to pay for it) or in Julialand through packages like [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl/tree/master/src/Spaces). If you're interested in more thorough manipulation of Orthogonal Polynomial spaces, they are probably your best bet right now.

So why another package?

First, I wanted to build my own package, and do not mind upstreaming it later on into SpecialFunctions.jl or ApproxFun.jl, if it so happens to be found a happy home.

Second of all, I wanted to learn by myself and have a lightweight standalone no-dependency repository where I could simply point to others if they wanted a package that does one thing, and does it very well. Most importantly, no other packages, as known to the author at time of writing, offer the full optimizations possible for user chosen parameters, e.g., offering Legendre 17th order with α = 3, β = 9.

Third, I wanted to compare my approach to established literature and other working code and see how I fared, all of this in parallel to producing a living document of videos and notes so that others may learn from my mistakes. As a comparison, here are some benchmarks for the 12th order Laguerre polynomial with different implementations:

```
# Hand-Crafted implementation of L₁₂ - writing out the entire polynomial
julia> @btime L12(.5)
  478.579 ns (0 allocations: 0 bytes)
-0.23164963886551918

# Approx Fun implementation
julia> @btime h12(.5)
  416.910 ns (4 allocations: 64 bytes)
-0.2316496388655192

# Naive recurrence formula
julia> @btime laguerre(12,.5)
  430.729 ns (52 allocations: 832 bytes)
-0.23164963886551906

# Recurrence formula guarding the 0th case
julia> @btime laguerre_corrected(12,.5)
  442.377 ns (52 allocations: 832 bytes)
-0.23164963886551906

# Written Horner method case for degree 12
julia> @btime goal12(.5)
  1.449 ns (0 allocations: 0 bytes)
-0.23164963886602852
```

The key takeaway is that we can generate efficient code (since all our polynomials' coefficients can be computed at compile time) and thus cut our computation time by almost 2 orders of magnitude. Note as well that because we are reducing dramatically the number of floating point operations, we are losing less precision than other methods.
It is also worth noting that the Julia JIT compiler performs quite well for simply writing the most simplistic code possible, as in the `laguerre` code. (Credits to `@yingboma`, who contributed that function in slack.) It doesn't scale into the `@horner` land of ridiculous speed due to not being spoonfed the appropriate structure to SIMD the generated code.

## Approach

This is how we will solve the problem as was hinted above in the roadmap.
There is a couple of things you need to learn about first.

1. How the `@horner` macro works, explained in the video link above.

2. How the formula found matches this approach, and how we our code will build up incrementally towards it.

3. Why using `@generated` functions is necessary at all.

The `@horner` macro works by taking a given list of coefficients to compute a polynomials and rearranging it algebraically so that it a) faster b) more precise, a classical result found in Knuth, TAOCP II. We will have to modify the macro just a tinybit so that it fits with our master formula.

The master formula was found by ~procrastination~ leafing through the book. After reading a discourse post by [Steven G Johnson](https://discourse.julialang.org/t/simple-metaprogramming-exercises-challenges/731/6) and looking at the formula, it is a clear match if you just rearrange it a bit. To wit:

$$ a\_{m-1}(x) = 1 - \frac{b\_m}{c\_m}f(x) a\_m(x)$$

Fits just right in the `muladd` patter if we rearrange the `1` to be the addition. This is why I found the idea of this package exciting: we can take classical orthogonal polynomials and give them a reach and speed not yet taken care of in Julia.

From here, we want to be able to produce functions for the coefficients that we want for a given polynomial. This means that, based on the polynomial we want, we will have to have different behavior for each coefficient, depending on the type of polynomial (this screams multiple dispatch! Do you see how?)

Eventually, we will want to produce these coefficients at compile time so that a new function can optimized by the JIT compiler when we actually call it. This will require the use of `@generated` functions, which can only see types as arguments. This means that any and all computation we want to do in our generated functions like cooking up the desired polynomial coefficients will have to be crammed into the types and the functions have to work on the types alone.

Lastly, The codebase in this package is very much alpha (development stage). It is meant to showcase the basic building blocks and should not be used in production.

Our starting point is to have a Minimal Working Example (MWE) of the vanilla Laguerre polynomial - even if it is a bit wonky.

First Lesson: build a basic example where your code runs, then build on top of it.

The code is all stashed in `src/OrthogonalPolynomials.jl` so go and have a look at the basic machinery before moving on to `v0.0.2`.

Be sure to look at the tests in `tests/runtests.jl` where we verified our workflow along the way.
