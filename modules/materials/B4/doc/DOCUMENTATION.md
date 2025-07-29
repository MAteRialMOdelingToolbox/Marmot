Linear Elastic Model based on B4 model by Bazant et al. (2015)
---
\page b4 B4 Model

## Implementation

The implementation can be found in \ref Marmot::Materials::B4.

## Theory 

The present model is based on the B4 model proposed by Bazant et al. \cite Bazant2015 generalized to 3D. 
The constitutive law is given in total form as

\f[
  \displaystyle \sig = \Cel : \epsE  = 
    \Cel : \left( \eps - \epsVE - \epsF - \epsDC -\epsSHR \right), 
\f]  

relating the nominal stress tensor \f$ \sig \f$ 
to the elastic strain tensor\f$ \epsE \f$ 
in terms of the fourth order stiffness tensor \f$ \Cel \f$.
The rate equations decribing the evolution of the respective 
strain components can be found a recent publication by Dummer et al. \cite Dummer2022.
However, nonlinear creep is not accounted for in the present formulation.
