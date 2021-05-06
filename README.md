# Code for Weights of Local Isomorphism Classes
This code is used to compute the weights in Tables 1-4 in an upcoming paper about counting D_4 quartic fields by conductor.

For each isomorphism class of quadratic extensions over Q_2, the code computes global representatives `a` of `K_p^x/K_p^{x2}` for each `p|2`. And, for each `a` computes the 2-part of the norm of the relative discriminant of `K(\sqrt{a})/K` and compares that to the 2-part of the discriminant of `Q(\sqrt{Nm(a)})`.

## Running the Code
The code was written in Sage 9.1 but likely would work with earlier versions. To run from the command line:
```
sage local_weights.sage
```
