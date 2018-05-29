1. Testing for all functions in an array:
@test all(f(theta)==1 for f in fs)
Credit to @chrisrackauckas

2. You an also `0 .|> fs` for "reverse broadcasting" a value to an array of functions.
Credit to @MichaelBorregaard

3. Generated functions ONLY SEE TYPE INFO. You gotta use dispatch for all the goodies.

4. SIMD needs a for loop - that will make it preferable to use the for loop approach over the recursive one!

5. Ternary operator requires a `()` around a return statemnt for it to work.
`x == zero(x) ? (return 1) : nothing

6. Use `@test_broken to know that your tests are broken, and you won't forget later!`

7. 
