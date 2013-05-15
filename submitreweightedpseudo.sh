function submit()
{
  qsub -cwd -S /bin/sh -l s_cpu=1:50:0 -o logreweightedpseudoexp.o -e logreweightedpseudoexp.e -V reweightedpseudo
}

rm logreweightedpseudoexp*

export APPROX_MEAS=FALSE  

export REWEIGHT_VAR=M


  export PSEUDOWEIGHT_ZPRIME=FALSE
  export PSEUDOWEIGHT_SQUARE=FALSE 
  export DONOTREWEIGHT=FALSE 

  export PSEUDOWEIGHT_INVERT=TRUE
  export PSEUDOWEIGHT_CENTER=TRUE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=TRUE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=FALSE
  export PSEUDOWEIGHT_CENTER=TRUE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=FALSE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=FALSE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=TRUE
  submit

  export PSEUDOWEIGHT_INVERT=TRUE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=TRUE
  submit
  
  export PSEUDOWEIGHT_SQUARE=TRUE # at this point the other variables dont matter anymore
  submit
  export PSEUDOWEIGHT_SQUARE=FALSE 
  
  export DONOTREWEIGHT=TRUE # at this point the other variables dont matter anymore
  submit
  export DONOTREWEIGHT=FALSE

export REWEIGHT_VAR=PT


  export PSEUDOWEIGHT_ZPRIME=FALSE
  export PSEUDOWEIGHT_SQUARE=FALSE 

  export PSEUDOWEIGHT_INVERT=TRUE
  export PSEUDOWEIGHT_CENTER=TRUE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=TRUE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=FALSE
  export PSEUDOWEIGHT_CENTER=TRUE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=FALSE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=FALSE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=TRUE
  submit

  export PSEUDOWEIGHT_INVERT=TRUE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=TRUE
  submit

  export PSEUDOWEIGHT_SQUARE=TRUE # at this point the other variables dont matter anymore
  submit
  export PSEUDOWEIGHT_SQUARE=FALSE 
  
  export DONOTREWEIGHT=TRUE # at this point the other variables dont matter anymore
  submit
  export DONOTREWEIGHT=FALSE

export REWEIGHT_VAR=Y

  export PSEUDOWEIGHT_ZPRIME=FALSE
  export PSEUDOWEIGHT_SQUARE=FALSE 

  export PSEUDOWEIGHT_INVERT=TRUE
  export PSEUDOWEIGHT_CENTER=TRUE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=TRUE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=FALSE
  export PSEUDOWEIGHT_CENTER=TRUE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=FALSE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=FALSE
  submit

  export PSEUDOWEIGHT_INVERT=FALSE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=TRUE
  submit

  export PSEUDOWEIGHT_INVERT=TRUE
  export PSEUDOWEIGHT_CENTER=FALSE
  export PSEUDOWEIGHT_NEGATIVE=TRUE
  submit
  
  export PSEUDOWEIGHT_SQUARE=TRUE # at this point the other variables dont matter anymore
  submit
  export PSEUDOWEIGHT_SQUARE=FALSE 
  
  export DONOTREWEIGHT=TRUE # at this point the other variables dont matter anymore
  submit
  export DONOTREWEIGHT=FALSE
