Linear Elastic Model based on B4 model by Bazant et al. (2015)
---
\page b4 B4 Model

## Implementation

The implementation can be found in \ref Marmot::Materials::B4

## Theory 

This model is based on the B4 model proposed by Bazant et al. \cite Bazant2015. The constitutive law is given in total form as

\f[
  \displaystyle \sig = \Cel : \epsE  = 
    \Cel : \left( \eps - \epsVE - \epsF - \epsDC -\epsSHR \right), \f]  

relating the nominal stress tensor \f$ \sig \f$ 
to the elastic strain tensor\f$ \epsE \f$ 
in terms of the fourth order stiffness tensor \f$ \Cel \f$.

The evolution of the strain components are governed by the following rate equations:

\f[ 
 \displaystyle \epsERate(t) =  q_1 \, \CelUnitInv \, \sigRate(t) \\
 \displaystyle \epsVERate (t)  =  \dfrac{1}{v(t)}\, \int_0^t \dot{ \Phi }(t-t')\, \CelUnitInv : \mbox{d} \sig (t')\\ 
 \displaystyle \epsFRate(t)  =  \dfrac{q_4}{t}\, \CelUnitInv : \boldsymbol{ \sigma } (t) \\
\f]
