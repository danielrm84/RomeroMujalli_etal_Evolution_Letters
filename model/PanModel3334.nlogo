;**************************************************************************************
;          Model Name:                       PanModel33
;
;                                  UNIVERSITY OF POTSDAM    (V1_0         )
;                                  UNIVERSITY OF GREIFSWALD (V2_0 and V3_0)
;
;
; Author(s):       Romero-Mujalli (2019, PhD thesis at University of Potsdam)
;                  Romero-Mujalli D. et al Ecology and Evolution 11 (2021)
;                  Romero-Mujalli D. et al Evolution Letters (2024)
;
; Written by:      Daniel Romero Mujalli
; Maintained by:   Daniel Romero Mujalli
;
; Written:         27.06.2016    (dd.mm.yy)
; Last update:     02.03.2024
;
; Type of model:   Agent-/Individual-based model (ABM, IBM)
;
; Summary:         The purpose of the model is:
;                  to study in-situ responses (genetic changes and
;                  phenotypic plasticity) (version 1_0)
;                  and the factors limiting or facilitating dispersal into new habitats
;                  (versions > 1_0)
;                  of different kinds of organisms (life history strategies) under
;                  scenarios of environmental change
;
;
;              NOTES / COMMENTS / QUESTIONS:
;
;
;
;
;**************************************************************************************
;**************************************************************************************
;                            MIT License
;
;                   Copyright (c) 2019 Daniel Romero Mujalli
;
;   Permission is hereby granted, free of charge, to any person obtaining a copy of
;   this software and associated documentation files (the "Software"), to deal in the
;   Software without restriction, including without limitation the rights to use, copy,
;   modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
;   and to permit persons to whom the Software is furnished to do so, subject to the
;   following conditions:
;
;   The above copyright notice and this permission notice shall be included in all
;   copies or substantial portions of the Software.
;
;   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
;   INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
;   PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
;   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
;   CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
;   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
;
;**************************************************************************************
; load extensions
extensions
[
;  table     ; to count unique genetic polymorphisms as given by x(tau = 0)
  ;r
  vid bitmap
]

;global parameters of the model:
globals
[
  DEBUG?                ; debug the code

  starting-seed         ; random seed being used

  genetic-mean          ; initial genetic mean of the population
  env-effect-mean       ; mean random environmental effect on phenotype (standard model and
                        ; implicit genetics only)
  initial-env-optimum   ; the initial phenotypic optimum as given by the environment
  size-of-environment   ; roaming range
  female-male-ratio     ; controls the sex ratio in the population. Example: a value of 0.7
                        ; means 70% females and 30% males in the population
  env-effect-variance   ; variance of environmental effect (standard model and implicit
                        ; genetics only)
  extinction?           ; boolean: true if the population is extinct
  pop-gv                ; population level additive genetic variance
  density-compensation  ; governs the population dynamics defining the life strategy
  strength-selection    ; strength of selection used in fitness function
  mu-dist-mut-effects   ; mean distribution of effect size from % of beneficial-mutations
  ;std-of-k              ; carrying capacity standard devitaton (for stochastic carrying capacity)
  ;time-limit            ; desired time limit (in generations)
  ; evolution of dispersal
  ;dispersal-mv          ; mutational variance of the dispersal trait (implicit model)
  max-dispersal-distance ; maximum distance a turtle can disperse

  ; STARTING LOCATION
  START-LOCATION   ; starting location of the population on the map

  ; NEUTRAL LANDSCAPE MODEL
  HEIGHT-MAP-SIZE  ; Map size, depends on map-scaling-factor
  HABITAT-GOOD     ; suitable habitat color code
  HABITAT-BAD      ; unsuitable habitat color code
  HABITAT-NEUTRAL  ; patch initial color -> white
  RN-MIN           ; random generator, lower bound
  RN-MAX           ; random generator, upper bound
  PROP-HABITAT-GOOD; for random uniform generated map only
  GRADIENT-SD      ; expected latitudinal deviation from mean env-optimum (at the bottom)

  CUSTOM-PATCH-SIZE; custom patch size
  CUSTOM-TURTLE-SIZE; custom turtle size
  CUSTOM-TURTLE-COLOR; custom color for turtles
  CUSTOM-TURTLE-SHAPE; custom shape for turtles
]

; properties of the individuals
turtles-own
[
  phenotype            ; phenotype 'z'
  genetic-component    ; genetic component 'a' of the phenotype
  environmental-effect ; environmental effect 'e'

  ; evolution of dispersal
  dispersal-propensity ; propensity to disperse (probability)

  ; evolution of plasticity
  sensibility-factor

  ; neutral marker
  neutral-marker

  fitness              ; fitness wi of individual i
  stress               ; individual stress = 1 - fitness (rho)
  fecundity            ; fecundity lambda of individual i
  stage                ; whether "adult" or "juvenile", as string
  sex                  ; sex: "female" or "male"
  reproduced?          ; boolean

  ; explicit genetics: haplotypes
  dna-strain1          ; first strain of the chromosome
  dna-strain2          ; second strain of the chromosome
  ; evolution of plasticity (haplotypes)
  reg-strain1
  reg-strain2
  ; neutral marker
  neutral-strain1
  neutral-strain2
]

; properties of the environment (patches)
patches-own
[
  optimum              ; phenotypic optimum tita as given by the environment
  optimum-0               ; initial optimum
  noise                ; variance of the environmental optimum
  degree-maladaptation ; degree of maladaptation
  time-to-extinction   ; store the time when the population went extinct
  ;mean-env-optimum     ; the mean environmental optimum
  carrying-capacity    ; proxy of "environmental quality"

  ;NEUTRAL LANDSCAPE MODEL
  rank                 ; value affecting the habitat type of the patch
]






;**************************************************************************************
;        TO SETUP
;**************************************************************************************
; This procedure
; - clears the interface, plots, and reset tickets
; - sets the global parameters values
; - initializes the environmental conditions
; - creates a locally adapted population of turtles
; For the model, the environment consists of a group of patches with colors. Typically,
; those on black are considered suitable habitats.
; Each individual is created with own values for the genetic and environmental
; components. These two components are used to calculate the phenotype.
to setup

  ; clear interface, reset ticks
  ;clear-all
  ;r:clear
  clear-turtles
  clear-drawing ; clear all lines and stamps drawn by turtles
  clear-output
  clear-all-plots
  reset-ticks

  ; set values for global parameters and setup the landscape
  set-global-parameters
  ;setup-landscape

  ; set the environment according to the type-of-model?:
  ; non-spatial -> in-situ responses only
  ; spatial-explicit -> accounts for dispersal
  if-else(type-of-model? = "non-spatial")
  [
  ; patch at the center (0,0) and those in radius 10 from the center are use as
  ; the environment to test for local adaptation
  ; more than one patch were selected for better visualization of the world
    ask patch 0 0
    [
      ;set pcolor green
      set optimum initial-env-optimum
      set optimum-0 optimum
      set noise 0
      ;set degree-maladaptation compute-degree-maladaptation
      set time-to-extinction 0
      ;set mean-env-optimum initial-env-optimum

      ; set carrying capacity (proxy of environmental quality)
      set carrying-capacity k-capacity ; initial conditions without stochasticity

      ask patches with [distance myself <= size-of-environment]
      [
        ;set pcolor green
        set optimum initial-env-optimum
        set optimum-0 optimum
        set noise 0
        ;set degree-maladaptation compute-degree-maladaptation
        set time-to-extinction 0
        ;set mean-env-optimum initial-env-optimum
        set carrying-capacity k-capacity
      ]
    ]
  ]; else: spatial-explicit:
  [
    ; change topology to a box
    ; possibilities: torus (default), box, vertical cylinder, horizontal cylinder
    __change-topology false false ; two inputs world wraps horizontally? wraps vertically?

    ; add a latitudinal environmental gradient
    ask patches with [ pxcor = 0 ]
    [
      ; impose a gradient with an environmental value of the upper most patch
      ; differing in 0.2 units standard-deviation from the initial (bottom patch)
      ; The environment is thought at a regional scale (e.g., Europe), and that
      ; elevation plays no role on the optimum phenotype as given by the
      ; environment
      set optimum initial-env-optimum + pycor / max-pycor * GRADIENT-SD
      set optimum-0 optimum
      set noise 0
      set time-to-extinction 0
      set carrying-capacity k-capacity
      ; there is no longitudinal gradient
      ask patches with [ pycor = [pycor] of myself ]
      [
        set optimum [optimum] of myself
        set optimum-0 optimum
        set noise 0
        set time-to-extinction 0
        set carrying-capacity k-capacity
      ]
    ]
    ; bad - unsuitable - habitats cannot support any turtles / individuals
    ask patches with [ pcolor = HABITAT-BAD ][ set carrying-capacity 0 ]
  ]

  ;***************************************************************************


  ; create a population of N individuals,
  ; initialize individuals (the turtles) and
  ; calculate phenotype (update phenotype)
  create-turtles population-size
  [
    set color CUSTOM-TURTLE-COLOR
    set shape CUSTOM-TURTLE-SHAPE
    set size CUSTOM-TURTLE-SIZE
    ; evolution of dispersal
    set dispersal-propensity 0 ; no dispersal
    ; evolution of plasticity
    ;set sensibility-factor 0   ; no plasticity

    ;update stage and sex
    update-stage-sex

    if-else (type-of-model? != "spatial-explicit")
    [
      move-to one-of patches with [ pcolor = HABITAT-GOOD ]
    ]
    [; else, if spatially explicit
      if(START-LOCATION = "NORTH")
      [
        move-to one-of patches with [ pcolor = HABITAT-GOOD and
                                      pycor > max-pycor - (max-pycor / 3)
                                    ]
      ]
      if(START-LOCATION = "SOUTH")
      [
        move-to one-of patches with [ pcolor = HABITAT-GOOD    and
                                      pycor < max-pycor / 3
                                    ]
      ]
      if(START-LOCATION = "CENTER")
      [
        move-to one-of patches with [ pcolor = HABITAT-GOOD    and
                                      pycor > max-pycor / 3 and
                                      pycor < max-pycor - (max-pycor / 3) and
                                      pxcor > max-pxcor / 3 and
                                      pxcor < max-pxcor - (max-pxcor / 3)
                                    ]
      ]
    ]

  ]

  ; set initial conditions according to the method for simulating genetics
  if-else (how-genetics? = "implicit")
  [ ask turtles [ set-initial-conditions-implicit-genetics ] ]
  [
    if-else (how-genetics? = "explicit")
    [ ask turtles [ set-initial-conditions-explicit-genetics ] ]
    [ error "undefined method for modelling genetics" stop ] ;catch exception
  ]

  ; set degree-maladaptation: used to calculate the genetic-load
  if(type-of-model? = "non-spatial")
  [
    ask patch 0 0 [ set degree-maladaptation 0 ]
  ]

  ;update phenotypic response of turtles
  ;ask turtles
  ;[
  ;  if-else (how-plasticity? = "standard-model")
  ;  [ update-phenotype           ] ; standard model
  ;  [ update-phenotypic-response ] ; account for plasticity
  ;]

  ; make sure that the no plasticity scenario will not pay plasticity costs
  ; due to a degree of plasticity value greater than zero (i.e., slope > 0)
  if(how-plasticity? = "standard-model" and slope > 0)
  [ set slope 0 ]


  ; update mean of the distribution of effect size
  ;if-else (beneficial-mutations != 0.5)
  ;[set mu-dist-mut-effects mean-dist (beneficial-mutations) (mut-effect-size)]
  ;[set mu-dist-mut-effects 0]
  set mu-dist-mut-effects 0

  ; if the model is spatially explicit then constant population is not allowed
  if(type-of-model? = "spatial-explicit" and N-constant? = true)
  [
    set N-constant? false
    print "Constant population is prohibited when type-of-model? = spatial-explicit"
  ]

  ;update-output:
  ;update-output

end






;**************************************************************************************
;        TO GO
;**************************************************************************************
to go

  ; interations counter
  ; Each iteration means one generation. The model
  ; uses non-overlapping generations:
  tick

  ; record video
  if (crt-video)
  [
    if (vid:recorder-status = "inactive" )
    [
      vid:start-recorder
    ]
    vid:record-view
  ]


  ; update fitness and fecundity of turtles:
  ; DEBUG
  if(DEBUG?)
  [
    print "update phenotypic response, fitness and fecundity ... "
    reset-timer
  ]

  ;*************************************************************
  ; TO GENETIC LOAD
  ; done for the population at tau = 0; no epigenetic mutations
  ; Thus, only the genetics is considered
  if (type-of-model? = "non-spatial")
  [
    ask turtles
    [
      update-phenotype
      check-fitness-negative-exponential
      set-fecundity-Bjoerklund
    ]
    ask patch 0 0
    [
      set degree-maladaptation compute-degree-maladaptation
    ]
  ]
  ; END OF GENETIC LOAD
  ;*************************************************************

  ask turtles
  [
    ; update phenotypic response at this time, when assumming no time lag
    ; in the development (maturation) of the trait. This is,
    ; the environment of development and fitness evaluation are the same
    if(time-lag-dev? = false)
    [
      ;update phenotypic response
      ;lazy coding: fixed plasticity
      ;set sensibility-factor 0.03
      ;end of lazy coding

      ; update the phenotypic response of turtles
      if-else (how-plasticity? = "standard-model")
      [ update-phenotype           ] ; standard model
      [ update-phenotypic-response ] ; account for plasticity
    ]; end of time-lag-dev?

    ; calculate fitness of each individual
    if-else( fitness-function = "Bjoerklund2009")
    [
      check-fitness-Bjoerklund2009
      ; DEBUG:
      ;print "fitness according to Bjoerklund"
    ]
    [ ; else, negative-exponential
       check-fitness-negative-exponential
    ]
    ; set individual stress
    set stress 1 - fitness

    ; calcualte fecundity of each individual
    set-fecundity-Bjoerklund
  ]

  ; DEBUG
  if(DEBUG?)
  [ type "done, time in seconds: " print timer ]

  ; update-output:
  ;DEBUG
  if(DEBUG?)
  [ print "update-output ... " reset-timer]

  if(plot-output? = true)
  [
    update-output
    plot-mean-fitness
  ]

  ;DEBUG
  if(DEBUG?)
  [ type "done, time in seconds: " print timer ]

  ; reproduction
  if(DEBUG?)
  [ print "reproduction ... " reset-timer ]

  if-else (evolution? = true)
  [
      ; store population-level additive genetic variance
      if (how-genetic-variance = "population-level")
      [
        if-else ( count turtles > 1)
        [
          set pop-gv variance [genetic-component] of turtles with [stage = "adult"]
        ]
        [ set pop-gv 0 ] ; else, if only one turtle
      ]

    ; reproduction
    if-else(N-constant?)
    [
      ask turtles [ do-reproduction ]
      ;if(count turtles with [stage = "juvenile"] < k-capacity ) ;
      while[count turtles with [stage = "juvenile"] < k-capacity ]
      [
        ask turtles [ do-reproduction ]
      ]

      if( count turtles with [stage = "juvenile" ] > k-capacity )
      [ ask n-of (count turtles with [stage = "juvenile" ] - k-capacity) turtles with [stage = "juvenile" ] [ die ] ]
    ]
    [; else, variable size

      ask turtles [ do-reproduction ]
    ]

  ]
  [;else, create turtles (every turtle replace itself)

    ask turtles
    [
      hatch random-poisson fecundity
      [
        ; set stage as juvenile
        set stage "juvenile"
        ; set random orientation of the newborn and move to a suitable patch
        set heading random 360

        ; inheritance of dispersal
        ;set dispersal-propensity random-normal d sqrt mutational-variance
      ]
    ]
  ]
  ; DEBUG
  if(DEBUG?)
  [ type "done, time in seconds: " print timer ]

  ; mortality of adults
  ;DEBUG
  if(DEBUG?)
  [ print "mortality-of-adults ... " reset-timer ]

  ask turtles with [ stage = "adult" ] [ die ]

  ;DEBUG
  if(DEBUG?)
  [ type "done, time in seconds: " print timer ]

  ; update stage and age, and phenotypic response (if time-lag-dev is active)
  ; evolution of dispersal
  ; check value range of dispersal-propensity
  ;DEBUG
  if(DEBUG?)
  [ print "update stage-sex and check value-range of dispersal-propensity" reset-timer ]
  ask turtles
  [
    ; dispersal propensity
    if (dispersal-propensity < 0 )[ set dispersal-propensity 0 ]
    if (dispersal-propensity > 1 )[ set dispersal-propensity 1 ]

    ;update stage and sex
    update-stage-sex

    ; update phenotypic response at this time, to account for a time lag
    ; in the development (maturation) of the trait. This is,
    ; the environment of development Et differs from the environment where
    ; fitness is evaluated (Et+1)
    if(time-lag-dev?)
    [
      ;update phenotypic response
      ;lazy coding: fixed plasticity
      ;set sensibility-factor 0.03
      ;end of lazy coding

      ; update the phenotypic response of turtles
      if-else (how-plasticity? = "standard-model")
      [ update-phenotype           ] ; standard model
      [ update-phenotypic-response ] ; account for plasticity
    ]; end of time-lag-dev?
  ]

  ;DEBUG
  if(DEBUG?)
  [ type "done, time in seconds: " print timer ]

  ; check time-to-extinction:
  ;DEBUG
  if(DEBUG?)
  [ print "check-extinction ... " reset-timer ]

  if (count turtles < 1 and extinction? = false)
  [
    ask patches with [ pcolor = HABITAT-GOOD ] [ set time-to-extinction (ticks - time-to-balance) ]
    set extinction? true
  ]

  ; end simulation:
  if (extinction? = true or (ticks - time-to-balance) > time-limit)
  [
    if(crt-video)
    [
      vid:save-recording "PanModel33.mp4"
      print vid:recorder-status
    ]
    stop
  ]

  ;DEBUG
  if(DEBUG?)
  [ type "done, time in seconds: " print timer ]

  ;**********************************************************************
  ; DISPERSAL HETEROGENEOUS / FRAGMENTED LANDSCAPES
  ;**********************************************************************
  ;DEBUG
  if(DEBUG?)
  [ print "Dispersal procedure ... " reset-timer ]

  if-else(type-of-model? = "non-spatial")
  [ ask turtles [move-to one-of patches with [ pcolor = HABITAT-GOOD ] ] ]
  [; add here the dispersal strategies
    ask turtles
    [
       ;if ( count neighbors with [pcolor = HABITAT-GOOD] > 0 )
       ;[ move-to one-of neighbors with [ pcolor = HABITAT-GOOD ] ]

       ; dispersal to patches within certain distance
       if( random-float 1 < dispersal-propensity )
       [
          move-to one-of patches with [pcolor = HABITAT-GOOD and
                                       distance myself < max-dispersal-distance ]
       ]
    ]
  ]
  ;DEBUG
  if(DEBUG?)
  [ type "done, time in seconds: " print timer ]
  ;**********************************************************************

  ; Environmental scenario:
  ; update-environment: updates optimum phenotype as given by the-environment according
  ; to the scenario of environmental change.
  ; If specified, the environmental change starts after a given number of
  ; iterations/generations (mutation selection balance/mutation-selection drift balance)
  ;**********************************************************************
  ; HETEROGENEOUS / FRAGMENTED LANDSCAPES
  ;**********************************************************************
  ;DEBUG
  if(DEBUG?)
  [ print "update-environment-before-rep-loop ... " reset-timer ]

  if-else( type-of-model? = "non-spatial" )
  [ ask patch 0 0 [ update-environment ] ]
  [ ask patches with [pcolor = HABITAT-GOOD] [ update-environment ] ]

  ;DEBUG
  if(DEBUG?)
  [ type "done, time in seconds: " print timer ]

  ;type "N: " print count turtles
  ;type "K: " print [carrying-capacity] of patch 0 0
  ;print "ratio: " print count turtles / [carrying-capacity] of patch 0 0

  ;**********************************************************************

  ; interations counter
  ; Each iteration means one generation. The model
  ; uses non-overlapping generations:
  ;tick


end






;**************************************************************************************
;        TO SET-GLOBAL-PARAMETERS
;**************************************************************************************
to set-global-parameters

  set DEBUG?               false   ; false: no DEBUG

  set CUSTOM-TURTLE-COLOR  orange ; (see netlogo color palette)
  set CUSTOM-TURTLE-SHAPE "circle"
  set CUSTOM-TURTLE-SIZE 1

  set genetic-mean         0.0
  set env-effect-mean      0.0 ; standard model and implicit genetics only

  set initial-env-optimum  0.0
  set female-male-ratio    0.5 ; 0.5 means random set of sex
  set extinction? false
  ; set the variance of the random environmental effect "ve" on the development of the
  ; trait. If how-plasticity? = "standard-model", then ve is computed according to given
  ; h2 and gv as in Bjoerklund et al (2009).
  if (how-plasticity? = "standard-model")
  [ set env-effect-variance   compute-env-effect-variance (heritability) (genetic-variance)]

  ; set level of density compensation, which is simulated as in Bjoerklund et al (2009),
  ; according to the density dependence effeect selected by the user: this impacts the
  ; population dynamics. For example, monotonic or oscillatory population dynamics.

  set density-compensation density-dependence-effect

  ; DEBUG:
  ;write "density dependence effect: " print density-dependence-effect
  ;write "value of density compensation: " print density-compensation

  ; maximum distance a turtle can disperse
  set max-dispersal-distance 3


  ; set the strength of selection according to fitness function and type of organism
  ; whether specialist, moderate, generalist
  if-else (fitness-function = "Bjoerklund2009")
  [
    ; values are set as in Bjoerklund et al (2009)
    if-else (type-organism = "specialist") [ set strength-selection 10 ]
    [ ;else
      if-else (type-organism = "moderate") [ set strength-selection 20 ]
      [;else
        if-else (type-organism = "generalist") [ set strength-selection 40 ]
        [;else, catch exception
          error "Unidentified type of organism ..." stop
        ]
      ]
    ]
  ]
  [;else: exponential function
    ; values are set as in Burger and Lynch (1995): strength of selection
    if-else (fitness-function = "negative-exponential")
    [
      if-else (type-organism = "specialist") [ set strength-selection 1 ]
      [ ;else
        if-else (type-organism = "moderate") [ set strength-selection 2.2 ]
        [;else
          if-else (type-organism = "generalist") [ set strength-selection 3.2 ]
          [;else, catch exception
            error "Unidentified type of organism ..." stop
          ]
        ]
      ]
    ]
    [;else; catch exception
      error "Unidentified fitness function" stop
    ]
  ]

  ; DEBUG
  ;write "strength of selection: " print strength-selection
  ;write "env-effect-variance: " print env-effect-variance

  ; end of parameter values
end







;**************************************************************************************
;        TO SETUP-LANDSCAPE
;**************************************************************************************
to setup-landscape
  ;**********************************************************************
  ; HETEROGENEOUS / FRAGMENTED LANDSCAPES
  ;**********************************************************************
  ; set parameter values
  ; roaming range
  ; Intraspecific interactions within distance size-of-environment are perceived for
  ; each turtle / individual. This affects its fecundity, for example
  set size-of-environment 10

  ; env-optimum expected deviation along the latitudinal range
  set GRADIENT-SD 1
  ; color code:
  ;(see netlogo documentation for more information on permitted values)
  set HABITAT-GOOD 0;52
  set HABITAT-BAD  19.9;48
  set HABITAT-NEUTRAL white

  ; proportion of habitat-good (random uniform generated maps only)
  set PROP-HABITAT-GOOD 0.5

  ; set lower, upper bound for the random generator
  set RN-MIN -1
  set RN-MAX  1

  if-else (type-of-model? = "non-spatial")
  [
    ; custom patch size
    set CUSTOM-PATCH-SIZE 3

    ; resize world to netlogo default values
    resize-world (- size-of-environment) size-of-environment (- size-of-environment) size-of-environment
    set-patch-size CUSTOM-PATCH-SIZE
    ; turn suitable the patches of the local environment
    ask patch 0 0
    [
      set pcolor HABITAT-GOOD
      ask patches with [distance myself <= size-of-environment]
      [ set pcolor HABITAT-GOOD ]
      ; turn HABITAT-BAD those outside radius given by size-of-environment
      ask patches with [distance myself > size-of-environment] [ set pcolor HABITAT-BAD ]
    ]
  ]
  [ ; spatially explicit:

    ; since carrying-capacity is modelled differently under explicit space
    ; make sure it has a reasonable value in favor of model performance
    ; Advance users may modify this at own risk
    if(k-capacity > 200)
    [
      set k-capacity 200
      print("warning! k-capacity adjusted in favor of model performance")
    ]

    ; set the starting location of the population
    ; one of "NORTH" "SOUTH" "CENTER"
    set START-LOCATION "SOUTH"

    ; custom patch size
    set CUSTOM-PATCH-SIZE 1

    ; Resize the world-map and its patches
    set size-of-environment 20
    set HEIGHT-MAP-SIZE (2 ^ map-scaling-factor) + 1
    set-patch-size CUSTOM-PATCH-SIZE
    resize-world 0 HEIGHT-MAP-SIZE 0 HEIGHT-MAP-SIZE
    ; set initial conditions for the map area
    ask patches [ set pcolor HABITAT-NEUTRAL ]

    ; generate the map according to selected method
    ;if-else (landscape-type? = "random")
    ;[ set-random-map  ]
    ;[ ; at the moment only works if the world is toroid
    ;  __change-topology true true
    ;  set-fractal-map
    ;]
  ]
  ;**********************************************************************
end






to crt-map
  ; generate the map according to selected method
  if-else (landscape-type? = "random")
  [ set-random-map  ]
  [ ; at the moment only works if the world is toroid
    __change-topology true true
    set-fractal-map
  ]
end
;**************************************************************************************
;        TO SET-INITIAL-CONDITIONS-IMPLICIT-GENETICS
;**************************************************************************************
; This function set the initial conditions for both, the genetic and environmental
; components of the phenotype for each individual (i.e., turtle).
; Here the genetics is implicitly modelled.
; Important!: random-normal function uses the standard deviation = sqrt(variance) of the
; distribution.
; The function works in a turtle context, example: ask turtles [ set-initial-cond... ]
to set-initial-conditions-implicit-genetics

  ; genetic component 'a':
  set genetic-component random-normal genetic-mean sqrt (genetic-variance)

  ; environmental component 'e':
  set environmental-effect random-normal env-effect-mean sqrt (env-effect-variance)

  ; DEBUG:
  ;write "turtle_" write who print ": "
  ;write "genetic component: " print genetic-component
  ;write "environ component: " print environmental-effect

end







;**************************************************************************************
;        TO SET-INITIAL-CONDITIONS-EXPLICIT-GENETICS
;**************************************************************************************
; this function sets the initial conditions for the genetic and environmental component
; of each individual turtle.
; The genetics is modelled explicitly, considering the number of loci.
; Three methods for the initialization of allele values were implemented. Only one must
; be active.
; alleles can be either:
; - real allele values set according to normal distribution as in Vincenzi (2014)
; - real allele values randomly set
; - random set of ones and zeros (on / off)
; (loci effect on phenotype is assumed to be additive).
; This function works in a turtle context, example: ask turtles [set-initial-condi... ]
to set-initial-conditions-explicit-genetics

  ; set allele values for each locus. Notice that the two dna strains are simulated
  ; separately

  ; population is assumed locally adapted. Allele values from a normal distribution
  ; as in Vincenzi (2014)
  ; N(0, va), where mean = 0 = initial environmental optimum, and
  ; va is the additive genetic variance per locus at the start of the simulation, given by
  ; VG (genetic-variance) / 2*L, L is the number-of-loci, and VG is the additive genetic
  ; variance of the quantitative trait at the start of the simulation
  ; (THIS METHOD is based on the continuum-of-alleles model of Kimura 1970: an introduction
  ; to population genetics theory (cited in Vincenzi 2014) )
  ; Vincenzi (2014) simulate additive effect through the addition of allele values

  ; factor to correct the allelic variance
  let f 2

  ;if-else (initial-dist? = "normal")
  ;[
    set dna-strain1 n-values number-of-loci [random-normal 0
                                           sqrt (genetic-variance / (f * number-of-loci)) ]
    set dna-strain2 n-values number-of-loci [random-normal 0
                                           sqrt (genetic-variance / (f * number-of-loci)) ]
  ;]
  ;[; else: uniform
    ; alleles take values from uniform distribution in range (-genetic-variance, genetic-variance)
    ;set dna-strain1 n-values number-of-loci [ genetic-variance - random-float (f * genetic-variance) ]
    ;set dna-strain2 n-values number-of-loci [ genetic-variance - random-float (f * genetic-variance) ]
  ;]


  ; another method: alleles take values of 1 or 0 (on / off)
  ; this method might not be appropriated for small number of loci
  ; set dna-strain1 n-values number-of-loci [ random 2 ]
  ; set dna-strain2 n-values number-of-loci [ random 2 ]

  ; Evolution of plasticity: Diallelic model
  ; pL loci of additive effects coding for the sensibility factor
  ; L: number-of-loci coding for the phenotypic trait (see above)
  ; p: proportion of L that define the total number of loci (pL)
  ; coding for the sensibility factor
  ; The value of each locus can be, either 0 or 1
  ; The results is the addition of all pL loci divided by a scaling parameter
  set reg-strain1 n-values (loci-reg-genome) [ 0 ]
  set reg-strain2 n-values (loci-reg-genome) [ 0 ]
  ; set the sensibility factor according to whether the organism is haploid or diploid
  if-else (haploid?)
  [ set sensibility-factor  abs (sum reg-strain1)  / reg-scaling-factor ]
  [ set sensibility-factor (abs (sum reg-strain1) + abs (sum reg-strain2))  / reg-scaling-factor ]

  ; neutral marker
  ; Note that it is initialized and scaled similarly as for the regulatory loci
  if(neutral-marker?)
  [
    set neutral-strain1 n-values (loci-reg-genome) [ 0 ]
    set neutral-strain2 n-values (loci-reg-genome) [ 0 ]

    ; set the value of the neutral marker accordingo the whether the orgnaism is
    ; haploid or diploid
    if-else (haploid?)
    [ set neutral-marker  abs (sum neutral-strain1) / reg-scaling-factor ]
    [ set neutral-marker (abs (sum neutral-strain1) + abs (sum neutral-strain2))  / reg-scaling-factor ]

  ]

  ; DEBUG
  ;type "dna-strain1: " print dna-strain1
  ;type "dna-strain2: " print dna-strain2
  ;show item (number-of-loci - 1) dna-strain1

  ; set the value for the genetic component 'a' according to the genome
  ; and to whether the organism is haploid or diploid
  ; loci effect is assumed to be additive
  if-else (haploid?)
  [ set genetic-component (sum dna-strain1) ] ; haploid
  [ set genetic-component (sum dna-strain1) + (sum dna-strain2) ] ; diploid

  ; set environmental component 'e':
  set environmental-effect random-normal env-effect-mean sqrt (env-effect-variance)

  ; reset the every 2nd dna strain for haploid organisms
  if( haploid?)
  [
    set dna-strain2     0
    set reg-strain2     0
    set neutral-strain2 0
  ]

  ; DEBUG:
  ;write "turtle_" write who print ": "
  ;write "genetic component: " print genetic-component
  ;write "environ component: " print environmental-effect
end






;**************************************************************************************
;        TO UPDATE-PHENOTYPE
;**************************************************************************************
; This function calculates (updates) the phenotype 'z' of the turtle following
; z = a + e
; where 'a' is the genetic-component and 'e' the environmental-effect on the phenotype
; The function works in a turtle context, example: ask turtles [ update-phenotype ]
to update-phenotype

  set phenotype genetic-component + environmental-effect

  ;DEBUG:
    ;write "phenotype z of turtle " write who print ":"
    ;print phenotype
    ;write "genetics: " print genetic-component
    ;write "env-effect: " print environmental-effect

end






;**************************************************************************************
;        TO UPDATE-PHENOTYPIC-RESPONSE
;**************************************************************************************
; updates the phenotype z according to the plastic response (selected method).
; The methods are similar, differing only in the parameter b, which is a sinusoidal
; function in method 2 and 3. This imposes limits to the plastic response of z depending on
; the amount of change that the organism is currently experimenting.
; Mehod 1 is based on the typical reaction norm approach, with a small modification.
; Methods assume a reference environment Q* where no plasticity is needed, and
; the phenotype z develops according to the genetic component a.
; 1) method 1: linear reaction norm
; z = a + bQt; where a is the genetic component (also breeding value); Qt the optimum
; phenotype as given by the environment; and b governs the plastic response (the slope)
;
; 2) method 2 and 3: sinusoidal reaction norms
; z = a + bQt; b = sin(abs Qt - a) These methods accounts for adaptive phenotypic
; plasticity, and consider that the plastic response has its limits
;
; 3) method 4: random plasticity
; common approach (e.g., Vincenzi 2014), where a value is draw randomly from a normal
; distribution
to update-phenotypic-response

  let a genetic-component
  let Qt [optimum] of patch-here
  let b 0
  let _e 0

   if-else (how-plasticity? = "linear-RN" ) ; method 1
   [
     ;set b 1
     set b slope
   ]
   [ ;else method 2
     if-else (how-plasticity? = "adaptive-sinusoidal")
     [
       ; this method uses the amount of change defined as the departure
       ; from the reference environment Q* (Q* = Qt, such that Qt - a = 0)
       set Qt Qt - a
       ; in this method, 1 / slope affects the plastic ability of the organism,
       ; such that small slope less plasticity, large slope, high plasticity,
       ; just as for the linear reaction norm
       let omega 1 / slope
       ;
       ; In Netlogo sin function assumes angle is given in degrees
       if (abs (omega * Qt) < pi)
       [
         let angle (abs (omega * Qt) * 180) / pi
         set b sin (angle)
       ]
     ] ; no negative effect
     [ ; else method 3
      if-else (how-plasticity? = "adaptive-logistic")
      [
        ; this method is similar to method 2, but assumes that the plastic
        ; response saturates after certain amount of change (similar to
        ; logistic function
        set Qt Qt - a
        ; in this method, 1 / slope affects the plastic ability of the organism,
        ; such that small slope less plasticity, large slope, high plasticity,
        ; just as for the linear reaction norm
        let omega 1 / slope
        ; check condition for saturation of the plastic response
        if (abs (omega * Qt) > pi / 2)
        [
          set Qt (pi / (2 * omega)) * (abs Qt) / Qt ; keep the sign
        ]
        ; In Netlogo sin function assumes angle is given in degrees
        let angle (abs (omega * Qt) * 180) / pi
        set b sin (angle)
      ]
      [ ;else method 4
        if-else (how-plasticity? = "random-noise")
        [
          set b random-normal 0 sqrt 1
          ; random plasticity is assumed noise around the breeding value a
          ; thus
          set Qt 1
        ]
        [; method 5
         if-else (how-plasticity? = "epigenetics" and how-genetics? = "explicit")
          [; evolution of plasticity
            ; based on the method demonstrated in the Rscript: studying plasticity
            ; with limits
            let epimutation-rate 0
            let tau 33           ; developmental time
            let l1 dna-strain1   ; original values in haplotypes are not modified
            let l2 dna-strain2   ; i.e., no epigenetic inheritance
            let var_e 0.0 ; variance of the residual environmental effect
                          ; independent of Qt
            let perturbation 0

            ; simulation of epigenetic regulations
            while [ tau > 0 ]
            [
              ; assess the perturbation
              set _e random-normal 0 sqrt(var_e)

              if-else(haploid?)
              [ set a (sum l1)            ]
              [ set a (sum l1) + (sum l2) ]

              set perturbation (a + _e) - Qt
              ; update value of epimutation-rate
              set epimutation-rate 1 - exp( - (perturbation ^ 2) * sensibility-factor)

              let loop-control genetic-bias; genetic-bias loci always active
              ;loop-control is an array index and in Netlogo indices begin from 0 not 1
              ; thus, genetic-bias = 1 means locus index 0 is always active
              ;       genetic-bias = 9, that loci index 0:8 are always active.
              ;       index 9 => locus 10 is under regulation
              while [loop-control < number-of-loci ]
              [
                ; locus acctivation / deactivation
                ; haplotype l1
                if ( random-float 1 < epimutation-rate)
                [
                  if-else (item loop-control l1 = 0)
                  [ set l1 replace-item loop-control l1 item loop-control dna-strain1 ]
                  [ set l1 replace-item loop-control l1 0 ] ; else, deactivated
                ]

                ; haplotype l2
                ; only for diploid organisms
                if (haploid? = false)
                [
                  if ( random-float 1 < epimutation-rate)
                  [
                    if-else (item loop-control l2 = 0)
                    [ set l2 replace-item loop-control l2 item loop-control dna-strain2 ]
                    [ set l2 replace-item loop-control l2 0 ] ; else, deactivated
                  ]
                ]

                ; update loop-control variable
                set loop-control loop-control + 1
              ] ; end of loci loop (epigenetic regulations)
              ; update developmental time
              set tau tau - 1; tau decreases in one unit every loop step
            ] ; end of outer most loop: developmental time
          ] ; end of method 5: epigenetics
          [; else catch exception
            print("epigenetics compatible with explicit genetics only")
            error "Unidentified method of plasticity ..." stop
          ]
        ]
      ]
     ]
   ]

  ; update plastic response
  set phenotype (a + _e) + (b * Qt)

  ; DEBUG:
   ;write "phenotype of" write who write ": " print phenotype
   ;write "a: " print a
   ;write "b: " print b
   ;write "Qt: " print Qt
   ;write "sigmoid(Qt): " print sigmoid (Qt)


end






;**************************************************************************************
;        TO UPDATE-STAGE-SEX
;**************************************************************************************
; update the stage of the individuals to adulthood, and set sex randomly.
; The function works in a turtle context, example: ask turtles [ update-stage-sex ]
to update-stage-sex

  ; set stage to adulthood:
  set stage "adult"
  ; set reproduced? to false
  set reproduced? false
  ; set sex randomly:
  if-else (random-float 1 <= female-male-ratio)
  [ set sex "female" ]
  [ set sex  "male"  ] ; else, male

end






;**************************************************************************************
;        TO CHECK-FITNESS-BJOERKLUND2009
;**************************************************************************************
; Fitness wi of individual i is calculated according to Björklund et al. (2009):
; wi = 1 - [(zi - tita(t))^2]/gamma
; where tita(t) is the optimum phenotype as given by the environment at time t, and
; gamma is the strength of selection
; the maximum fitness is 1
; The function works in a turtle context
to check-fitness-Bjoerklund2009

  if-else([pcolor] of patch-here = HABITAT-GOOD)
  [
    set fitness (1 - (((phenotype - [optimum] of patch-here) ^ 2) / strength-selection))
  ]
  [ error "turtle outside suitable patches" stop ] ; change this for heterogeneous landscapes

  ; DEBUG
  ;write "fitness of turtle " write who write ": " print fitness
  ;write "optimum of patch-here: " print [optimum] of patch-here
  if (fitness > 1)
  [
    write "phenotype of turtle : " write who write ": " print phenotype
    write "optimum: " print [optimum] of patch-here
    write "gamma: " print strength-selection
    write "fitness: " print fitness
    error "warning! fitness value is greater than 1" stop
  ]

end






;**************************************************************************************
;        TO CHECK-FITNESS-NEGATIVE-EXPONENTIAL
;**************************************************************************************
; Fitness wi of individual i is calculated according to a negative exponential function
; as in Burger & Lynch (1995):
; wi = exp [ - (zi - tita(t))^2]/2*gamma^2
; where tita(t) is the optimum phenotype as given by the environment at time t, and
; gamma is the strength of selection.
; wi is in range (0, 1) => maximum fitness of 1
; The function works in a turtle context
to check-fitness-negative-exponential

  let tita 0

  if-else([pcolor] of patch-here = HABITAT-GOOD)
  [ set tita [optimum] of patch-here ]
  [ error "turtle outside suitable patches!" stop]

  let gamma strength-selection

  set fitness exp ( - ( (phenotype - tita) ^ 2) / (2 * (gamma) ^ 2) )

  ; DEBUG
  ;write "fitness of turtle " write who write ": " print fitness
  ;write "optimum of patch-here: " print [optimum] of patch-here

  ; catch error:
  if (fitness > 1)
  [
    write "phenotype of turtle : " write who write ": " print phenotype
    write "optimum: " print [optimum] of patch-here
    write "gamma: " print strength-selection
    write "fitness: " print fitness
    error "warning! fitness value is greater than 1" stop
  ]

end






;**************************************************************************************
;        TO SET-FECUNDITY-BJOERKLUND
;**************************************************************************************
; Fecundity lambda is calculated according to Björklund et al. (2009) for each
; individual:
; lambda = wi*exp[alpha(1 - N/K)]
; where wi is the fitness of individual i, alpha is the density-compensation, N the
; population size, K the carrying capacity and exp the exponential function
; The function works in a turtle context
to set-fecundity-bjoerklund

  let alpha density-compensation
  let N 0
  if-else(type-of-model? = "non-spatial")
  [ set N count turtles with [ stage = "adult" ] ]
  [ set N count turtles with [ stage = "adult" and distance myself <= size-of-environment ] ]
  let K [carrying-capacity] of patch-here
  let wi fitness

  if-else(N-constant?)
  [
    set fecundity wi * alpha
  ]
  [; else, variable
    set fecundity wi * ( exp ( alpha * (1 - N / K ) ) )
  ]

  ; evolution of plasticity
  ; adjust fecundity based on plasticity costs
  ; Individuals pay a cost proportional to their degree of plasticity
  ; this is, the sensibility-factor (mechanistic plasticity based on epigenetics),
  ; the slope (reaction-norm approaches)
  if-else(how-plasticity? = "epigenetics")
  [
    set fecundity fecundity - sensibility-factor * cost
  ]
  [;else
    set fecundity fecundity - slope * cost
  ]


  ; DEBUG
  ;write "fecundity of turtle " write who write ": " print fecundity
  ;write "population size: " print N
  ;write "K: " print K

end






;**************************************************************************************
;        DO-REPRODUCTION
;**************************************************************************************
; this function handles the reproduction of turtles according to the selected
; reproductive strategy (or mode):
; - "lottery-polygyny"
; - "hermaphrodite"
to do-reproduction

  if(reproductive-mode = "lottery-polygyny")
    [
      ; lottery polygyny: females reproduce only once, but males can repeat
      if( stage = "adult" and sex = "female" )
      [
        if-else (how-genetics? = "implicit")  [ reproduce-implicit-genetics ]
        [ if-else (how-genetics? = "explicit")[ reproduce-explicit-genetics ]
          [ error "undefined method of how-genetics?" stop ] ; catch exception
        ]
      ]
    ]

    if(reproductive-mode = "hermaphrodite")
    [
      if (stage = "adult")
      [
        if-else (how-genetics? = "implicit")  [ reproduce-implicit-genetics ]
        [ if-else (how-genetics? = "explicit")[ reproduce-explicit-genetics ]
          [ error "undefined method of how-genetics?" stop ] ; catch exception
        ]
      ]
    ]

end







;**************************************************************************************
;        TO REPRODUCE-IMPLICIT-GENETICS
;**************************************************************************************
; This function implements sexual reproduction according to Björklund et al. (2009).
; Female individuals randomly select a partner of opposite sex to mate. The fitness of
; a pair is the sum of the fitness value of the two parents. Then fecundity is calculated
; considering density dependence. The number of offspring is drawn from a poisson
; distribution centered on the value of fecundity.
;
; Genetics is implicit according to the infinitesimal model of quantitative genetics which
; assumes that traits are affected by a large number of loci of additive effects.
; Therefore trait inheritance can be approximated using a normal distribution with mean
; parents trait value, and variance: half the genetic variance.
; Mutations and recombination are implicit.
;
; The function works in a turtle context, example: ask turtles [ reproduce ]
to reproduce-implicit-genetics

  ; SEXUAL REPRODUCTION:
  let my-partner 0

  ; lottery polygyny
  ; pick a random partner of opposite sex:
  if(reproductive-mode = "lottery-polygyny")
  [
    set my-partner one-of turtles with [ stage = "adult" and
                                       sex = "male"    and
                                       distance myself <= size-of-environment ]
  ]

  ; hermaphrodite
  if(reproductive-mode = "hermaphrodite")
  [
    set my-partner one-of turtles with [ stage = "adult" and
                                       ;sex = "male"    and
                                       distance myself <= size-of-environment ]
  ]

  ; else, ASEXUAL REPRODUCTION
  if(asexual? = true)
  [ set my-partner self ]
  ; test whether there is a partner available to reproduce:
  if ( my-partner != nobody )
  [
    ; DEBUG:
    ;write "me "print who
    ;write "partner " print [who] of my-partner
    ;type "my-sex: " print sex
    ;type "partner sex: " print [sex] of my-partner
    ;if(reproduced? = true) [ print "me"]
    ;if([reproduced?] of my-partner = true) [print "partner"]

    ;create a list containing the genetic-component of both parents:
    let list-genetic-parents list (genetic-component) ([genetic-component] of my-partner)

    ;DEBUG:
    ;write "my-genetic-component: " print genetic-component
    ;write "my-partner-genetic-comp: " print [genetic-component] of my-partner
    ;write "list-genetic-components-of-parents: " print list-genetic-parents

    ; calculate the mean parental genetic-component which will be inherited by offspring
    let genetic-mean-of-parents mean ( list-genetic-parents )

    ;DEBUG:
    ;write "my genetic-component: " print genetic-component
    ;write "partner genetic-component: " print [genetic-component] of my-partner
    ;write "genetic-mean: " print genetic-mean-of-parents

    ; calculate genetic variance among parents:
    let genetic-variance-of-parents variance ( list-genetic-parents )
    ; DEBUG:
    ;write "genetic-variance-of-parents: " print genetic-variance-of-parents
    ;type "adults "print count turtles with [stage = "adult"]
    ;type "total " print count turtles
    ; set the genetic-variance for the offspring according to the selected method
    let my-genetic-variance 0

    if-else ( how-genetic-variance = "parental-level")
    [
      set my-genetic-variance (genetic-variance-of-parents);parental-level
                                                           ;(Bjoerklund2009)
    ]
    [; else,
      if-else (how-genetic-variance = "parameter")
      [ set my-genetic-variance genetic-variance ] ; as in Reed et al. (2011):
      ;the population-level additive genetic variance is an input parameter (set by user)
      [; else
        if-else (how-genetic-variance = "population-level")
        [
        ; method as in Vincenzi & Piotti (2014): genetic variance is the total
        ; additive genetic variance for the trait at the population level
          set my-genetic-variance pop-gv
        ]
        ; else, catch exception
        [ error "undefined method for genetic variance" stop ]
      ]
    ]
    ; the genetic variance is simulated as in Vincenzi et al (2012). The infinitesimal
    ; model is modified to account for the decline in additive genetic variance and the
    ; new input of variation through mutation (mutational-variance).
    ; This method also accounts for offspring additive genetic variance equals half the
    ; additive genetic variance
    set my-genetic-variance (1 / 2) * (my-genetic-variance + mutational-variance)

    ;-------------------------------------------------------------------------------------
    ; INHERITANCE OF DISPERSAL
    let d mean list (dispersal-propensity) ([dispersal-propensity] of my-partner)
    ;-------------------------------------------------------------------------------------

    ; calculate the fecundity of the breeding pair pair-fecundity as in Björklund et al.
    ; Björklund(2009) assumed that the fitness of the pair was the sum of the fitness
    ; values of the two parents (wsum) which is equivalent to the sum of their fecundities.
    ; Given that wsum is used, the population is saturated above the carrying capacity.
    let fecundity-of-pair 0
    if(reproductive-mode = "lottery-polygyny")
    [
      set fecundity-of-pair ( fecundity + [ fecundity ] of my-partner )

      ; adjust the fecundity if using constant population size (N)
      if(N-constant?)
      [
        set fecundity-of-pair fecundity-of-pair / (2 * mean [fecundity] of turtles)
      ]
    ]

    if(reproductive-mode = "hermaphrodite")
    [
      set fecundity-of-pair ( fecundity )

      ; adjust the fecundity if using constant population size (N)
      if(N-constant?)
      [
        set fecundity-of-pair fecundity-of-pair / (mean [fecundity] of turtles)
      ]
    ]

    ;DEBUG:
    ;write "my-fecundity: " print fecundity
    ;write "fecundity-of-partner: " print [fecundity] of my-partner
    ;write "fecundity-of-pair: " print fecundity-of-pair

    let number-of-offspring random-poisson fecundity-of-pair
    ; DEBUG
    ;write "nr. offspring: " print number-of-offspring

    if (number-of-offspring >= 1) ; reproduction occurs
    [
      hatch number-of-offspring
      [ ; according to netlogo dictionary, each new turtle inherits all its
        ; variables, including its location, from its parent, except [who]

        ; set stage as juvenile
        set stage "juvenile"
        ; set random orientation of the newborn and move to a suitable patch
        set heading random 360
        ;move-to one-of neighbors with [ pcolor = HABITAT-GOOD ]

        set genetic-component random-normal genetic-mean-of-parents ; mean
                                sqrt my-genetic-variance

        ; DEBUG
        ;type "genetic-component of offspring: " print genetic-component
        ;type "genetic-variance: "  print my-genetic-variance

         ; environmental component 'e':
        ; Set the variance of environmental effect ve according to method of heritability
        ; if h2 = fixed, compute the variance of environmental effect according to the
        ; genetic variance and level of heritability
        ; else, ve is a constant parameter
        let ve 0 ; initialization
        if (how-plasticity? = "standard-model")
        [ set ve compute-env-effect-variance (heritability) (my-genetic-variance)]

        set environmental-effect random-normal env-effect-mean sqrt ve

        ; inheritance of dispersal
        set dispersal-propensity random-normal d sqrt mutational-variance

        ; DEBUG:
        ;write "genetic-variance: " print my-genetic-variance
        ;write "env-effect-variance: " print ve
      ] ; end of inheritance

   ] ; end of hatching n offspring ( reproduction)

    ; update reproduced? status of the partner to true
    ask my-partner [ set reproduced? true ]
    ;type "genetic-component of partner " print [genetic-component] of my-partner
  ]; end of if statement: whether a partner is available for reproduction

  ; update reproduced? of this turtle or individual to true:
    set reproduced? true
    ;type "my-genetic-component " print genetic-component

end






;**************************************************************************************
;        TO REPRODUCE-EXPLICIT-GENETICS
;**************************************************************************************
; This function implements sexual reproduction and lottery polygyny.
; Females randomly select a partner of opposite sex to mate. The fitness of a pair
; is the sum of the fitness value of the two parents. Then fecundity is calculated
; considering density dependence. The number of offspring is drawn from a poisson
; distribution centered on the value of fecundity.
;
; Genetics is explicit. This means that loci are explicitly modelled. Loci effect on
; trait value is assumed to be additive.
; Recombination is simulated through the random selection of parental allele values.
; Then, mutations take place with probability mut-rate.
;
; The function works in a turtle context, example: ask turtles [ reproduce ]
to reproduce-explicit-genetics

  ; store the loci (dna strains) of this turtle (parent 1)
  let my-loci1 dna-strain1
  let my-loci2 dna-strain2
  ; evolution of plasticity
  let my-reg-strain1 reg-strain1
  let my-reg-strain2 reg-strain2
  ; neutral marker
  let my-neutral-s1 0
  let my-neutral-s2 0
  if (neutral-marker?)
  [
    set my-neutral-s1 neutral-strain1
    set my-neutral-s2 neutral-strain2
  ]

  ; SEXUAL REPRODUCTION:
  let my-partner 0

  ; lottery polygyny
  ; pick a random partner of opposite sex:
  if(reproductive-mode = "lottery-polygyny")
  [
    set my-partner one-of turtles with [ stage = "adult" and
                                       sex = "male"    and
                                       distance myself <= size-of-environment ]
  ]

  ; hermaphrodite
  if(reproductive-mode = "hermaphrodite")
  [
    set my-partner one-of turtles with [ stage = "adult" and
                                       ;sex = "male"    and
                                       distance myself <= size-of-environment ]
  ]

  ; else, ASEXUAL REPRODUCTION
  if(asexual? = true)
  [ set my-partner self ]

  ; test whether there is a partner available to reproduce:
  if ( my-partner != nobody )
  [
    ; DEBUG:
    ;write "me "print who
    ;write "partner " print [who] of my-partner
    ;type "my-sex: " print sex
    ;type "partner sex: " print [sex] of my-partner
    ;if(reproduced? = true) [ print "me"]
    ;if([reproduced?] of my-partner = true) [print "partner"]

    ;-------------------------------------------------------------------------------------
    ; INHERITANCE OF DISPERSAL
    ; evolution of dispersal
    let d mean list (dispersal-propensity) ([dispersal-propensity] of my-partner)

    ; INHERITANCE OF PLASTICITY
    ; evolution of plasticity
    ;let mean-sensibility mean list (sensibility-factor) ([sensibility-factor] of my-partner)
    ;let mu-rate 0.001
    ;-------------------------------------------------------------------------------------

    ; store loci of my-partner (parent 2)
    let partner-loci1 [dna-strain1] of my-partner
    let partner-loci2 [dna-strain2] of my-partner
    ; evolution of plasticity (parent 2)
    let partner-reg-strain1 [reg-strain1] of my-partner
    let partner-reg-strain2 [reg-strain2] of my-partner
    ; neutral marker (parent 2)
    let partner-neutral-s1 0
    let partner-neutral-s2 0
    if(neutral-marker?)
    [
      set partner-neutral-s1 [neutral-strain1] of my-partner
      set partner-neutral-s2 [neutral-strain2] of my-partner
    ]

    ;DEBUG:
    ;write "Loci parent 1: " write my-loci1 write " and " print my-loci2
    ;write "Loci parent 2: " write partner-loci1 write " and " print partner-loci2

    ; calculate the fecundity of the breeding pair pair-fecundity as in Björklund et al.
    ; Björklund(2009) assumed the fitness of the pair as the sum of the fitness
    ; values of the two parents (wsum) which is equivalent to the sum of their fecundities.

    ; fecundity adjustments depending on the reproductive mode such that the population does
    ; not exceed the carrying capacity of the environment
    let fecundity-of-pair 0
    if(reproductive-mode = "lottery-polygyny")
    [
      set fecundity-of-pair ( fecundity + [ fecundity ] of my-partner )

      ; adjust the fecundity if using constant population size (N)
      if(N-constant?)
      [
        set fecundity-of-pair fecundity-of-pair / (2 * mean [fecundity] of turtles)
      ]
    ]

    if(reproductive-mode = "hermaphrodite")
    [
      set fecundity-of-pair ( fecundity )

      ; adjust the fecundity if using constant population size (N)
      if(N-constant?)
      [
        set fecundity-of-pair fecundity-of-pair / (mean [fecundity] of turtles)
      ]
    ]

    ;DEBUG:
    ;write "my-fecundity: " print fecundity
    ;write "fecundity-of-partner: " print [fecundity] of my-partner
    ;write "fecundity-of-pair: " print fecundity-of-pair

    ; set the genetic variance. Only for the calculation of environmental variance
    ; according to the heritability
    let genetic-variance-of-parents variance list (genetic-component)
                                                  ([genetic-component] of my-partner)
    let my-genetic-variance 0

    if-else ( how-genetic-variance = "parental-level")
    [
      set my-genetic-variance (genetic-variance-of-parents); parental-level
                                                           ; (Bjoerklund2009)
    ]
    [; else,
      if-else (how-genetic-variance = "parameter")
      [ set my-genetic-variance genetic-variance ] ; as in Reed et al. (2011):
        ;the population-level additive genetic variance is an input parameter (set by user)
      [; else
        if-else (how-genetic-variance = "population-level")
        [
        ; method as in Vincenzi & Piotti (2014): genetic variance is the total
        ; additive genetic variance for the trait at the population level
          set my-genetic-variance pop-gv
        ]
        ; else, catch exception
        [ error "undefined method for genetic variance" stop ]
        ]
    ]
    ; additive genetic variance of offspring equals half the additive genetic variance:
    set my-genetic-variance (1 / 2) * (my-genetic-variance)

    let number-of-offspring random-poisson fecundity-of-pair
    ; DEBUG
    ;write "nr. offspring: " print number-of-offspring

    if (number-of-offspring >= 1) ; reproduction occurs
    [
      hatch number-of-offspring
      [
        ; set stage as juvenile
        set stage "juvenile"
        ;set color blue
        ; set random orientation of the newborn and move to a suitable patch
        set heading random 360
        ;move-to one-of neighbors with [ pcolor = HABITAT-GOOD ]

        ; evolution of dispersal
        ; inheritance of dispersal
        set dispersal-propensity random-normal d sqrt ((1 / 2) * (genetic-variance + mutational-variance))

        ; evolution of plasticity
        ; inheritance of plasticity
        ;set sensibility-factor random-normal mean-sensibility sqrt (mutational-variance)
        ;if(random-float 1 < mu-rate)
        ;[ set sensibility-factor one-of [0 0.1 0.5 1 2 5 10 20] ]


        ;*********************************************************************************
        ; INHERITANCE GENETIC COMPONENT OF PHENOTYPE
        ;*********************************************************************************
        ; inheritance with recombination and mutations, explicit genetics:
        ; this routine considers the ploidy level of the organism

        ;1st mutation takes place during the creation of gametes inside each parent
        ; Each locus is checked for mutations according to the given mutation rate.
        ; The effect-size of the mutation is implemeted using a normal distribution
        ; N(0, variance), where the variance controls for the effect-size (Vincenzi 2014).
        ;
        ; A mutation can be explicitly simulated per locus or as a single mutation per
        ; strain with probability number-of-loci * mut-rate-per-locus, as in Vincenzi (2014)
        ; the effect-size of mutations can be modified in the future. For example,
        ; Vincenzi (2014) based the effect-size of mutations on a normal distribution
        ; N(0, effect-size). Through modifications of the variance of the distribution
        ; one can control the effect size of the mutation. Note that using the normal
        ; distribution, the expected value is the mean (in the example, 0 effect).

        ;DEBUG:
        ;write "parent1 strain1: " print my-loci1
        ;write "parent1 strain2: " print my-loci2

        ;parent1
        let strain1-parent1 evaluate-mutations(my-loci1)(number-of-loci)("infinite")(mut-rate-per-locus)(mut-effect-size)
        let strain2-parent1 0
        if(haploid? = false)
        [set strain2-parent1 evaluate-mutations(my-loci2)(number-of-loci)("infinite")(mut-rate-per-locus)(mut-effect-size)]

        ;parent2
        let strain1-parent2 evaluate-mutations(partner-loci1)(number-of-loci)("infinite")(mut-rate-per-locus)(mut-effect-size)
        let strain2-parent2 0
        if(haploid? = false)
        [set strain2-parent2 evaluate-mutations(partner-loci2)(number-of-loci)("infinite")(mut-rate-per-locus)(mut-effect-size)]

        ; 2nd inheritance with recombination:

        ; in the model, haploid organisms that undergo sexual reproduction are assumed to form a
        ; diploid zygote that split into haploid cells by meiosis that then make the haploid organism
        ; Thus, recombination could still take place
        if-else (haploid?)
        [
          ;offspring strain1
          set dna-strain1 do-recombination (strain1-parent1) (strain1-parent2) (number-of-loci)
        ]
        [; else, diploid
          ;offspring strain1
          set dna-strain1 do-recombination (strain1-parent1) (strain2-parent1) (number-of-loci)
          ;offsptring strain2
          set dna-strain2 do-recombination (strain1-parent2) (strain2-parent2) (number-of-loci)
        ]


        ; DEBUG
        ;write "offspring strain1: " print dna-strain1
        ;write "offspring strain2: " print dna-strain2

        ; evolution of plasticity
        ; mutation and recombination of the regulatory loci
        let genome-size loci-reg-genome
        let s1p1 0
        let s2p1 0
        let s1p2 0
        let s2p2 0
        let mu mut-rate-per-locus * (10 ^ reg-mu-scaling)

        if(how-plasticity? = "epigenetics")
        [
          ;parent1
          set s1p1 evaluate-mutations(my-reg-strain1)(genome-size)(allele-model)(mu)(mut-effect-size)
          if(haploid? = false)
          [set s2p1 evaluate-mutations(my-reg-strain2)(genome-size)(allele-model)(mu)(mut-effect-size)]

          ;parent2
          set s1p2 evaluate-mutations(partner-reg-strain1)(genome-size)(allele-model)(mu)(mut-effect-size)
          if(haploid? = false)
          [set s2p2 evaluate-mutations(partner-reg-strain2)(genome-size)(allele-model)(mu)(mut-effect-size)]

          ; inheritance with recombination
          ; in the model, haploid organisms that undergo sexual reproduction are assumed to form a
          ; diploid zygote that split into haploid cells by meiosis that then make the haploid organism
          ; Thus, recombination could still take place
          if-else (haploid?)
          [
            ;offspring strain1
            set reg-strain1 do-recombination (s1p1) (s1p2) (genome-size)
          ]
          [; else, diploid
            ;offspring strain1
            set reg-strain1 do-recombination (s1p1) (s2p1) (genome-size)
            ;offsptring strain2
            set reg-strain2 do-recombination (s1p2) (s2p2) (genome-size)
          ]

        ]
        ; neutral marker (similar to regulatory region)
        ; using same genome-size
        if(neutral-marker?)
        [
          ;parent1
          set s1p1 evaluate-mutations(my-neutral-s1)(genome-size)(allele-model)(mu)(mut-effect-size)
          if(haploid? = false)
          [set s2p1 evaluate-mutations(my-neutral-s2)(genome-size)(allele-model)(mu)(mut-effect-size)]

          ;parent2
          set s1p2 evaluate-mutations(partner-neutral-s1)(genome-size)(allele-model)(mu)(mut-effect-size)
          if(haploid? = false)
          [set s2p2 evaluate-mutations(partner-neutral-s2)(genome-size)(allele-model)(mu)(mut-effect-size)]
          ;recombination
          if-else (haploid?)
          [
            ;offspring strain1
            set neutral-strain1 do-recombination (s1p1) (s1p2) (genome-size)
          ]
          [; else, diploid
            ;offspring strain1
            set neutral-strain1 do-recombination (s1p1) (s2p1) (genome-size)
            ;offsptring strain2
            set neutral-strain2 do-recombination (s1p2) (s2p2) (genome-size)
          ]
        ]

        ; end of inheritance explicit genetics
        ;***********************************************************************************

        ;***********************************************************************************

        ; evolution of plasticity
        ; update sensibility factor of the offspring
        ; according to whether the organism is haploid or diploid
        ; sensibility factor: explicitly result from additive effets of regulatory genome
        ;                     and is always >= 0
        if-else (haploid?)
        [
           set sensibility-factor  abs (sum reg-strain1)  / reg-scaling-factor
        ]
        [; else, diploid
           set sensibility-factor (abs (sum reg-strain1) + abs (sum reg-strain2))  / reg-scaling-factor
        ]


        ; neutral marker
        ; similar scaling as for the regulatory region
        ; consider whether the organism is haploid or diploid
        if(neutral-marker?)
        [
          if-else (haploid?)
          [ set neutral-marker  abs (sum neutral-strain1)  / reg-scaling-factor                              ]
          [ set neutral-marker (abs (sum neutral-strain1) + abs (sum neutral-strain2))  / reg-scaling-factor ]

        ]

        ; update the genetic component of the offspring
        ; accordingo to the ploidy level of the organism
        if-else (haploid?)
        [ set genetic-component (sum dna-strain1)                     ]
        [ set genetic-component (sum dna-strain1) + (sum dna-strain2) ]


        ; DEBUG
        ;type "genetic-component: " print genetic-component

        ; environmental component 'e':
        ; Set the variance of environmental effect ve according to method of heritability
        ; if h2 = fixed, compute the variance of environmental effect according to the
        ; genetic variance and level of heritability
        ; else, ve is a constant parameter
        let ve 0 ; initialization
        if (how-plasticity? = "standard-model")
        [ set ve compute-env-effect-variance (heritability) (my-genetic-variance)]

        set environmental-effect random-normal env-effect-mean sqrt ve

        ; DEBUG:
        ;write "genetic-variance: " print my-genetic-variance
        ;write "env-effect-variance: " print ve
      ] ; end of inheritance

   ] ; end of hatching n offspring ( reproduction)

    ; update reproduced? status of the partner to true
    ask my-partner [ set reproduced? true ]

  ]; end of if statement: whether a partner is available for reproduction

  ; update reproduced? of this turtle or individual to true:
    set reproduced? true

end





;**************************************************************************************
;        TO-REPORT EVALUATE-MUTATIONS
;**************************************************************************************
; This method returns the initial genome strain after checking for mutations
; according to given mutation-rate
; Mutated values in the strain are set based on the selected method:
; "diallelic" values are either 0 or 1
; "infinite" values are continuos real numbers from a normal distribution centered
; on zero and standard deviation given by the effect-size of mutations (input parameter)
; genome-size denotes the number of loci in the strain
to-report evaluate-mutations [strain genome-size method mut-rate effect-size]

  let l 0
  if-else (method = "diallelic")
  [
    while[l < genome-size]
    [
      if(random-float 1 < mut-rate)
      [
        set strain replace-item l strain one-of [0 1]
      ]
      set l l + 1
    ]
  ]
  [; else, method 2:
    if-else(method = "infinite")
    [
      while[l < genome-size]
      [
        if(random-float 1 < mut-rate)
        [
          let dummy (item l strain) + (random-normal 0 sqrt effect-size)
          set strain replace-item l strain dummy
        ]
        set l l + 1
      ]

    ]
    [
      error "unidentified method of evaluate-mutations"
    ]
  ]

  report strain

end





;**************************************************************************************
;        TO-REPORT DO-RECOMBINATION
;**************************************************************************************
; This method returns the result s of recombining the two input strains (s1 and s2)
; Loci are assumed unlinked
; genome-size denotes the number of loci in each strain
; Strains should have the same length!
to-report do-recombination [s1 s2 genome-size]

  let l 0
  ; initialization of the resulting strain
  let s n-values genome-size [0]

  while[l < genome-size ]
  [
    set s replace-item l s (one-of (list item l s1
                                         item l s2))
    set l l + 1
  ]

  report s

end





;**************************************************************************************
;        TO UPDATE-ENVIRONMENT
;**************************************************************************************
; Currently, there are two main environmental scenarios: climate-change and cyclic,
; respectively.
; In climate change scenario,
; the mean environmental optimum changes at rate r every iteration. For the model, each
; iteration is equivalent to a generation
; The optimum Q updates according to:
; 1) Qt = Q0 + rt (Directional deterministic)
; 2) Qt = Q0 + rt; Qt* = Qt + E, (Directional stochastic)
; where E is the type of noise defined as:
; Et = aE(t-1) + b*Dt, where D = N(0, 1); here, a is not the genetic component!
; Here b is assumed to be b = sqrt VEt*(1 - a^2), as in Schwager et al (2006), where
; VEt is the environmental variance at time t; and the
; parameter a (set by the user), determines the strength of the autocorrelation
; (based on Ripa and Lundberg 1996, and Schwager et al. 2006)
;      a = 0 (white noise)
;  0 < a < 1 (red noise) Björklund et al used a = 0.7
; -1 < a < 0 (blue noise) Björklund et al used a = -0.7
;
; The parameter VEt can change or not in time, depending on the rate of change k of
; the variance. Thus, VEt = VE + kt, where VE is the inital environmental variance.
; This consideration allows to simulate the increased probability of extreme events
; as predicted by climatic IPCC scenarios, as in Vincenzi et al (2012).
; The function works in a patch context, example: ask patches [ update-optimum ]
;
; In the cyclic environment,
; The optimum Qt changes according to a sinusoidal function that considers two
; parameters:
; the amplitude A, and the period T, thus:
; Qt = A.sin(2.pi.t / T), according to Burger and Krall (2004)
to update-environment

  ; initial optimum of the environment Q0 = tita-0 according to Björklund2009
  let tita-0 optimum-0
  let new-optimum 0
  let new-noise 0
  let VE env-variance
  let k rate-change-of-env-variance
  let autocorr level-autocorr
  let b sqrt (1 - (autocorr ^ 2)) ; scaling factor of the variance according to
                           ; Schwager et al (2006)

  ; first catch exception:
  if (scenario? != "climate-change" and scenario? != "cyclic")
  [ error "update-environment: undefined scenario of environment" stop ]

  ; update the optimum according to the scenario of environment
  if (scenario? = "climate-change")
  [
    ; set parameter values for the scenario climate-change
    let r rate-change-of-optimum

    ; update moving mean envirnomental optimum
    if (ticks > time-to-balance)
    [ set new-optimum tita-0 + (r * (ticks - time-to-balance) ) ]; 1) directional deterministic
  ]

  if (scenario? = "cyclic") ; according to Burger and Krall 2004
  [
    ; set parameter values for the cyclic scenario
    let A amplitude ; amplitude of the wave
    let T period ; period of the wave T = 1 / fr; fr is the frequency

    ; update the environmental optimum
    if(ticks > time-to-balance)
    [
      let arg (2 * pi * (ticks - time-to-balance)) / T ; radians
      ; Netlogo sin function takes degrees
      set arg arg * 180 / pi
      set new-optimum ( A * sin ( arg ) )
    ]
  ]

  ; check for increasing environmetal variance and update the environmental variance
  if (ticks > time-to-balance)
  [set VE VE + (k * (ticks - time-to-balance) )]

  ; update the value of the mean-env-optimum of the patch before adding stochastic noise
  ;set mean-env-optimum new-optimum

  ; apply stochasticity around the environmental optimum
  ; stochastic environment with potentially different colour of noise
  ; depending on the level of outocorrelation a (level-autocorr)
  ;type "noise before " print noise
  set new-noise (autocorr * noise) + (b * (sqrt VE) * (random-normal 0 sqrt 1))
  set new-optimum new-optimum + new-noise
  ;type "noise after " print noise

  ; update optimum and noise of the patch
  set optimum new-optimum
  set noise new-noise

  ;**********************************************************************
  ; HETEROGENEOUS / FRAGMENTED LANDSCAPES
  ;**********************************************************************
  if(type-of-model? = "non-spatial")
  [
    ask patches with [distance myself <= size-of-environment]
    [
      set optimum new-optimum
      set noise new-noise
    ]
  ]
  ;**********************************************************************

  ; DEBUG
  ;type "new optimum: " print optimum
  ;type "new noise: "   print noise
  ;type "patch optimum " print optimum
  ;type "patch noise: " print noise
  ;type "env-variance: " print env-variance
  ;type "new env-variance: " print VE

end






;**************************************************************************************
;        TO-REPORT SET-STOCHASTIC-K-CAPACITY
;**************************************************************************************
; Here the carrying capacity K is considered a state variable subject to stochastic
; change in time.
; This function adds stochastic white noise around the current carrying capacity of the
; patch (passed as an argument).
; The noise is given by a normal distribution with mean 0 and standard deviation
; defined as a percent of the mean k-capacity which is user defined
;to-report set-stochastic-k-capacity
;
;  let var k-capacity * std-of-k / 100
;
;  let new-K  random-normal k-capacity sqrt (var ^ 2)
;
; ;Debug:
; ;type "carrying-capacity: " print new-K
;
; ; since negative K makes no sense,
; report (abs new-K)
;
;end






;**************************************************************************************
;        TO-REPORT COMPUTE-ENV-EFFECT-VARIANCE
;**************************************************************************************
; This function reports the corresponding value of environmental effect variance given the
; specified values of heritability and genetic variance
; h2: narrow-sense heritability
; gv: additive genetic variance
to-report compute-env-effect-variance [h2 gv]

  ; catch exception h2 <= 0:
  if-else (h2 <= 0)
  [ error "heritability muss be greater than 0 (h2 > 0)" ]
  [ report (gv / h2 ) - gv ] ; else

end






;**************************************************************************************
;        TO-REPORT COMPUTE-DEGREE-MALADAPTATION
;**************************************************************************************
; this function compute the degree of maladaptation based on the genetic load
; defined as
; genetic-load = abs(W0max - W0_avg) / W0max
; where W0avg and W0max are, respectively, the average and highest fitness in a
; population without epigenetic alterations (i.e., at developmental tau = 0; all
; loci active)
to-report compute-degree-maladaptation

 if-else ( count turtles > 0 )
 [
    let max_wi max   [ fecundity ] of turtles
    let mean_w mean [ fecundity ] of turtles

   ; compute the degree-maladaptation
   report ( max_wi - mean_w ) / max_wi
 ]
 [ report 0 ] ; else (extinct population)

end





;**************************************************************************************
;        TO IMPORT-MAP
;**************************************************************************************
; import pcolors from file-name
; if file-name is empty / not valid, print error message
; user may want to adjust the colors manually, e.g.,
; ask patches with [pcolor > 41 and pcolor < 71] [ set pcolor HABITAT-GOOD]
to import-map

  ; select directory
  set-current-directory user-directory

  ; if file-name does not exist pop-up error message
  if-else ( not empty? file-name and file-exists? file-name )
  [
    ; import map
    import-pcolors file-name
  ]
  [; else
    type "unable to find file " type file-name print " in current directory"
  ]


end






;**************************************************************************************
;        TO SET-RANDOM-MAP
;**************************************************************************************
; randomly (Bernoulli distribution) create a proportion p of suitable habitats
to set-random-map

  ;let p 0.5 ; proportion suitable habitats
  let habitats PROP-HABITAT-GOOD * (max-pxcor + 1) * (max-pycor + 1)
  ask n-of habitats patches [ set pcolor HABITAT-GOOD ]

end

;**************************************************************************************
;        TO SET-FRACTAL-MAP
;**************************************************************************************
; generate a random map based on the diamond-square algorithm (fractal model). The
; method was developed following the youtube videos by Klayton Kowalski and
; Mathematics of Computer Graphics and Virtual Environments
to set-fractal-map

  ; number of iteration: max_i ~ HEIGHT-MAP-SIZE
  let i 0
  ; habitat contagion
  let habitat-contagion 7

  ; h: ROUGHNESS ; the parameter H in With et al. (1997)
  let h 2 ^ (-2 * i * ROUGHNESS)

  ;** INITIALIZATION   ***********
  ; set a random value to each corner of the map
  ask patches with [(pxcor = 0 and pycor = 0) or
                    (pxcor = 0 and pycor = max-pycor) or
                    (pxcor = max-pxcor and pycor = max-pycor) or
                    (pxcor = max-pxcor and pycor = 0)
                   ]
  [
    set rank h * random-uniform RN-MIN RN-MAX
  ]

  ; set initial chunk size: height_map_size - 1
  let chunk-size HEIGHT-MAP-SIZE - 1
  ; set the step (called half) to walk each square / diamond
  let half 0
  ;** END OF INITIALIZATION ******

  ;** UPDATE PATCH RANK **********
  ; based on diamond-square method
  while [chunk-size > 1]
  [
    set i i + 1 ; update iterations
    set h 2 ^ (-2 * i * ROUGHNESS)
    set half chunk-size / 2

    ; square-step
    square-step h half chunk-size
    ; diamond-step
    diamond-step h half chunk-size

    set chunk-size chunk-size / 2
  ]
  ;** END OF UPDATE PATCH RANK ***

  ;** DRAW MAP *******************

  let interval count patches / N-HABITATS
  ; loop control variable
  let k 1
  ; color variable for habitat contagion
  let current-color 131
  ; steps used to walk the rank values and set the corresponding habitat type (color)
  let lower-bound 0 ;min [rank] of patches
  let upper-bound interval ;lower-bound + interval
  while [k <= N-HABITATS ]
  [
    if-else ( k = 1 )
    [
      foreach sublist sort-on [rank] patches lower-bound upper-bound
      [ the-patch -> ask the-patch [set pcolor HABITAT-GOOD ] ]
      ; update current-color
      set current-color HABITAT-GOOD
    ]
    [ if-else ( k = N-HABITATS )
      [
         foreach sublist sort-on [rank] patches lower-bound upper-bound
        [ the-patch -> ask the-patch [set pcolor HABITAT-BAD ] ]
         ;set current-color HABITAT-BAD ; drop error due to agentset with only one element :S
      ]
      [;else
        foreach sublist sort-on [rank] patches lower-bound upper-bound
        [ the-patch -> ask the-patch [set pcolor HABITAT-GOOD + k * 1] ] ;+10,new colorname
        set current-color HABITAT-GOOD + k * 1
      ]
    ]

    ; simulate spatial habitat contagion to smooth the landscape
    ; patches surounded by suitable neighbors turn suitable as well
    ask patches with
    [count neighbors with [pcolor = current-color] >= habitat-contagion]
    [ set pcolor 131 ]
    ; turn black all pink patches
    ask patches with [pcolor = 131] [ set pcolor current-color ]

    ; update walking steps and loop-control variable before repeatin the loop
    set lower-bound upper-bound
    set upper-bound upper-bound + interval
    set k k + 1
  ]

  ; habitat contagion to smooth the landscape of the HABITAT-BAD
  ask patches with
  [count neighbors with [pcolor = HABITAT-BAD] >= habitat-contagion]
  [ set pcolor 131 ]
  ; turn black all pink patches
  ask patches with [pcolor = 131] [ set pcolor HABITAT-BAD ]
  ;** END OF DRAW MAP ************
end

;**************************************************************************************
;        TO SQUARE-STEP
;**************************************************************************************
; perform the square step of the diamond-square algorithm
; arguments to the function:
; h: range of the random generator. Already accounts for the roughness
; half: current step value used to walk each square
; chunk-size: current chunk-size value (a chunk is the size of a square)
to square-step [ h half chunk-size ]

  ; initialization at position (1,1)
  let x0 1
  let y0 1
  ;print ("in")
  while [ x0 < max-pxcor - 1 ]
  [
    ; set the rank of the middle point patch (of the current square)
    ask patch (x0 + half) (y0 + half)
    [
      set rank (h * random-uniform -1 1) +
     ([rank] of patch (x0) (y0) +
      [rank] of patch (x0) (y0 + chunk-size) +
      [rank] of patch (x0 + chunk-size) (y0 + chunk-size) +
      [rank] of patch (x0 + chunk-size) (y0)) / 4
    ]
    ; move to next square
    set y0 y0 + chunk-size
    if( y0 > max-pycor - 1 )
    [
      set y0 0
      set x0 x0 + chunk-size
    ]
    ;type "chunksize" print chunk-size
    ;type "y:" print y0
    ;type "x:" print x0
  ]
  ;print("out")

end

;**************************************************************************************
;        TO DIAMOND-STEP
;**************************************************************************************
; perform the diamond step of the diamond-square algorithm
; arguments to the function:
; h: range of the random generator. Already accounts for the roughness
; half: current step value used to walk each square
; chunk-size: current chunk-size value (a chunk is the size of a square)
to diamond-step [ h half chunk-size ]

  ; initialization at position (1, half step)
  let x0 1
  let y0 (x0 + half) mod chunk-size
  ; result and n-count are used such that the mean rank of neighboring patches
  ; does not include patches outside map range (i.e., undefined patches)
  let result 0
  let n-count 0

  ;print("in")
  while [ x0 < max-pxcor ]
  [
    if (y0 = 0)[ set y0 chunk-size ]

    ; set the rank value of the middle point of the diamond
    ; omit patches outside map range, and
    ; update the sum result (and count, for average)
    if(x0 - half > 0)
    [
      set result result + [rank] of patch (x0 - half) (y0)
      set n-count n-count + 1
    ]
    if(y0 - half > 0)
    [
      set result result + [rank] of patch (x0) (y0 - half)
      set n-count n-count + 1
    ]
    if(y0 + half <= HEIGHT-MAP-SIZE)
    [
      set result result + [rank] of patch (x0) (y0 + half)
      set n-count n-count + 1
    ]
    if(x0 + half <= HEIGHT-MAP-SIZE)
    [
      set result result + [rank] of patch (x0 + half) (y0)
      set n-count n-count + 1
    ]

    ask patch (x0) (y0)
    [
      set rank (h * random-uniform -1 1) + (result / n-count)
    ]

    ; find the next diamond
    set y0 y0 + chunk-size
    if( y0 > max-pycor )
    [
      set y0 (x0 + half) mod chunk-size
      set x0 x0 + half
    ]

    ; reset result and n-count
    set result 0
    set n-count 0
  ]

end

;**************************************************************************************
;        TO-REPORT RANDOM-UNIFORM
;**************************************************************************************
; custom random generator in range (min-value, max-value)
; random generator according to uniform distribution
to-report random-uniform [min-value max-value]

  report min-value + random (max-value + 1 - min-value)

end





;**************************************************************************************
;        TO UPDATE-OUTPUT
;**************************************************************************************
; this function plots the distribution of phenotypes in the population
to update-output

  ;type "total " print count turtles
  ;type "adults " print count turtles with [stage = "adult"]
  ;type "juveniles " print count turtles with [stage = "juvenile"]
  let dummy 0 ; control count turtles = 0

  ; plot time series:
  ;let threshold 0.01
  set-current-plot "Time series"
  set-current-plot-pen "turtles"
  plot count turtles
  ;set-current-plot-pen "pen-1"
  ;plot count turtles with [sensibility-factor >  threshold]
  ;set-current-plot-pen "pen-2"
  ;plot count turtles with [sensibility-factor <= threshold]

  ; evolution of plasticity
  ; sensibility-factor omega
  set-current-plot "genetic-load"
  if (count turtles > 1 and type-of-model? = "non-spatial")
  [
    set-current-plot-pen "load"
    plot [degree-maladaptation] of patch 0 0
  ]

  set-current-plot "trait-evolution"
  if (count turtles > 1)
  [
    set-current-plot-pen "plasticity"
    plot mean [sensibility-factor] of turtles
    set-current-plot-pen "plasticity_median"
    plot median [sensibility-factor] of turtles

    ; neutral marker
    if(neutral-marker?)
    [
      set-current-plot-pen "neutral-marker"
      plot mean [neutral-marker] of turtles
     ; set-current-plot-pen "neutral-mode"
     ; plot item 0 modes [neutral-marker] of turtles
    ]
  ]

  ; plot environmental optimum vs mean phenotypic response
  set-current-plot "environmental optimum vs mean phenotypic response"
  set-current-plot-pen "optimum"
  plot mean [optimum] of patches with [ pcolor = HABITAT-GOOD ]
  set-current-plot-pen "mean z"
  if (count turtles > 0)
  [ plot mean [phenotype] of turtles ]

  ; plot phenotypic variance
  set-current-plot "geno/phenotypic variance"
  ; VP
  set-current-plot-pen "VP"
  set dummy 0
  if (count turtles > 1)
  [ set dummy variance [phenotype] of turtles ]
  plot dummy
  ; VG
  set-current-plot-pen "VG"
  set dummy 0
  if (count turtles > 1)
  [ set dummy variance [genetic-component] of turtles ]
  plot dummy

  ; store the highest current phenotipic value in the population
  ; this value will be used below to set the range of the x axis
  ;let highest-value [phenotype] of max-one-of turtles [phenotype]
  set-current-plot "frequency distribution"


  ;set-plot-x-range -1 end-point-of-environment
  set-plot-pen-mode 1        ; (0 for lines, 1 for bars)
  set-histogram-num-bars 10
  ; plot the phenotypic frequency in the population:
  set-current-plot-pen "phenotype"
  histogram [phenotype] of turtles
  ; plot the genotypic frequency in the population:
  set-current-plot-pen "genotype"
  if-else(how-genetics? = "explicit")
  [
    if-else(haploid?)
    [histogram [(sum dna-strain1)] of turtles]
    [histogram [(sum dna-strain1) + (sum dna-strain2)] of turtles]
  ]
  [histogram [genetic-component] of turtles]
  ; plot the environmental optimum
  set-current-plot-pen "env-optimum"
  plot-pen-reset
  set-plot-pen-mode 2
  ; plot a vertical line to show the current optimum phenotype in the environment
  let n 0
  while [n < k-capacity]
  [
    plotxy mean [optimum] of patches with [pcolor = HABITAT-GOOD ] n
    set n n + 0.1
  ]

end






;**************************************************************************************
;        TO PLOT-MEAN-FITNESS
;**************************************************************************************
; plot mean fitness of turtles
to plot-mean-fitness

  ;type "total " print count turtles
  ;type "adults " print count turtles with [stage = "adult"]
  ;type "juveniles " print count turtles with [stage = "juvenile"]

  set-current-plot "average stress"
  set-current-plot-pen "stress"
  let dummy 0
  if (count turtles > 0 )
  [ set dummy mean [stress] of turtles ]
  plot dummy

end




;**************************************************************************************
;        CHANGE LOG
;**************************************************************************************
; V3334:
; Extends the model to
; - haploid / diploid organisms
; - reproduction modes: lottery-polygyny / hermaphrodite
; - scaling factor for the mutation rate of the regulatory / neutral loci
; - constant / variable population size
;
; - Implementation of constant population size based on Botero et al 2015, PNAS
;   Applies to implicit and explicit genetics
; - Created: function do-reproduction
;
; - one locus for the plasticity trait, subjecto to the same mutation rate and
;   mutation effect size as for the response trait
; - The plasticity alleles are continuous real values
; - The sensitivity trait omega is real number > 0
; - Fixed: there was an error by the process evaluate-mutations. Only the first locus
;   was mutating under the infinite method. Also, the effect was wrong, since it needs
;   to add a random normaly distributed value to the current allele value
; - set genotypic variance as the variance contained in the allele values of turtles
;
; - implementation of time-lag-dev? which indicates whether the environment of trait
;   development and fitness evaluation are the same
; - code maintenance
; - check whether the high VG means coexistence of genetic polymorphisms
;   using
;   extensions [ table ]
;   length table:keys table:counts [genetic-component] of turtles
;
; - new parameter: genetic-bias (q)
;   q bias: L - q; 0 < q < L; q belongs to natural numbers
;   The condition q = L => genetic determinism, i.e., no epigenetic regulation
;
; - added calculation and plotting of the genetic-load (see also degree-maladaptation)
@#$#@#$#@
GRAPHICS-WINDOW
964
17
1035
89
-1
-1
3.0
1
10
1
1
1
0
0
0
1
-10
10
-10
10
0
0
1
ticks
30.0

BUTTON
100
10
169
43
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
7
44
398
191
frequency distribution
NIL
freq
-15.0
15.0
0.0
10.0
true
true
"" ""
PENS
"phenotype" 1.0 0 -13791810 true "" ""
"env-optimum" 1.0 2 -2674135 true "" ""
"genotype" 1.0 0 -16777216 true "" ""

BUTTON
7
10
99
43
redo-map
clear-all\nreset-ticks\nask patches [set pcolor white]\nsetup-landscape
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
170
10
296
43
Run the model!
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
7
191
398
338
time series
time in generations
N
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"turtles" 1.0 0 -13791810 true "" ""

SLIDER
400
81
603
114
k-capacity
k-capacity
10
1000
1000.0
10
1
NIL
HORIZONTAL

SLIDER
793
144
949
177
population-size
population-size
10
1000
1000.0
10
1
NIL
HORIZONTAL

CHOOSER
793
178
949
223
type-organism
type-organism
"specialist" "moderate" "generalist"
1

INPUTBOX
400
460
510
520
genetic-variance
1.0
1
0
Number

SLIDER
400
215
601
248
rate-change-of-optimum
rate-change-of-optimum
0
1
0.0
0.01
1
NIL
HORIZONTAL

CHOOSER
399
141
603
186
scenario?
scenario?
"climate-change" "cyclic"
1

SLIDER
779
740
938
773
heritability
heritability
0.1
1
1.0
0.1
1
NIL
HORIZONTAL

PLOT
7
631
216
773
genetic-load
time in generations
variance
0.0
2.0
0.0
0.05
true
true
"" ""
PENS
"load" 1.0 0 -16777216 true "" ""

CHOOSER
793
258
949
303
fitness-function
fitness-function
"Bjoerklund2009" "negative-exponential"
1

INPUTBOX
400
567
510
627
time-to-balance
0.0
1
0
Number

CHOOSER
654
728
777
773
how-genetic-variance
how-genetic-variance
"parameter" "parental-level" "population-level"
1

SLIDER
606
85
791
118
env-variance
env-variance
0
2
0.0
0.1
1
NIL
HORIZONTAL

PLOT
7
487
398
631
average stress
time in generations
stress
0.0
5.0
0.0
1.0
true
false
"" ""
PENS
"stress" 1.0 0 -16777216 true "" ""

CHOOSER
400
521
510
566
how-genetics?
how-genetics?
"implicit" "explicit"
1

INPUTBOX
511
506
646
566
mut-rate-per-locus
1.0E-4
1
0
Number

INPUTBOX
511
567
646
627
mut-effect-size
1.0
1
0
Number

SLIDER
511
472
646
505
number-of-loci
number-of-loci
1
50
10.0
1
1
NIL
HORIZONTAL

TEXTBOX
402
46
793
81
....................................... ECOLOGY .......................................
14
54.0
1

TEXTBOX
515
416
793
450
...................... EVOLUTION ........................
14
93.0
1

SLIDER
297
10
398
43
time-limit
time-limit
100
2000
2000.0
100
1
NIL
HORIZONTAL

SLIDER
606
119
791
152
level-autocorr
level-autocorr
-0.9
0.9
0.0
0.1
1
NIL
HORIZONTAL

SLIDER
606
153
791
186
rate-change-of-env-variance
rate-change-of-env-variance
0
0.1
0.0
0.01
1
NIL
HORIZONTAL

TEXTBOX
463
119
536
137
Scenario:
14
54.0
1

TEXTBOX
807
10
947
28
Type of Organism
14
14.0
1

TEXTBOX
535
446
779
464
.............. Explicit genetics .................
14
93.0
1

TEXTBOX
658
646
770
667
Implicit genetics
14
93.0
1

INPUTBOX
654
667
777
727
mutational-variance
0.1
1
0
Number

PLOT
7
337
398
487
environmental optimum vs mean phenotypic response
time
optimum
0.0
1.0
0.0
1.0
true
true
"" ""
PENS
"optimum" 1.0 0 -10899396 true "" ""
"mean z" 1.0 0 -14454117 true "" ""

CHOOSER
777
533
936
578
how-plasticity?
how-plasticity?
"standard-model" "random-noise" "linear-RN" "adaptive-sinusoidal" "adaptive-logistic" "epigenetics"
5

SLIDER
777
579
936
612
slope
slope
0.5
2
0.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
666
64
765
82
Stochasticity
14
54.0
1

TEXTBOX
405
192
791
210
Directional / Climate change\t        Cyclic environmental change
14
54.0
1

SLIDER
614
214
791
247
amplitude
amplitude
0
5
2.0
1
1
NIL
HORIZONTAL

SLIDER
614
249
791
282
period
period
1
200
200.0
1
1
NIL
HORIZONTAL

TEXTBOX
404
64
604
85
....... The Environment .............
14
54.0
1

TEXTBOX
782
718
931
736
if standard model,
14
93.0
1

SLIDER
793
224
949
257
density-dependence-effect
density-dependence-effect
0
5
1.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
783
444
930
462
            Plasticity
14
53.0
1

CHOOSER
400
293
601
338
type-of-model?
type-of-model?
"spatial-explicit" "non-spatial"
1

SLIDER
400
338
601
371
map-scaling-factor
map-scaling-factor
1
9
8.0
1
1
NIL
HORIZONTAL

CHOOSER
717
345
837
390
landscape-type?
landscape-type?
"random" "fractal"
1

SLIDER
717
390
837
423
ROUGHNESS
ROUGHNESS
0
1
0.69
0.01
1
NIL
HORIZONTAL

TEXTBOX
717
307
784
338
or create \none
12
0.0
1

TEXTBOX
407
270
593
294
Neutral Landscape Model
14
54.0
1

PLOT
397
631
651
773
trait-evolution
time
omega
0.0
10.0
0.0
0.1
true
true
"" ""
PENS
"plasticity" 1.0 0 -955883 true "" ""
"neutral-marker" 1.0 0 -7500403 true "" ""
"plasticity_median" 1.0 0 -12087248 true "" ""

INPUTBOX
601
344
715
404
file-name
switzerland.png
1
0
String

SLIDER
400
371
601
404
N-HABITATS
N-HABITATS
2
5
2.0
1
1
NIL
HORIZONTAL

BUTTON
602
311
715
344
NIL
import-map
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
635
10
791
43
Run-Model-Once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
400
426
510
459
Evolution?
Evolution?
0
1
-1000

SWITCH
399
10
507
43
plot-output?
plot-output?
0
1
-1000

SLIDER
647
560
776
593
loci-reg-genome
loci-reg-genome
1
20
1.0
1
1
NIL
HORIZONTAL

INPUTBOX
647
465
776
525
reg-scaling-factor
42.0
1
0
Number

SLIDER
777
499
936
532
cost
cost
0
1
0.0
0.1
1
NIL
HORIZONTAL

SWITCH
647
594
777
627
neutral-marker?
neutral-marker?
0
1
-1000

SWITCH
793
30
949
63
asexual?
asexual?
1
1
-1000

SWITCH
793
64
949
97
haploid?
haploid?
1
1
-1000

CHOOSER
793
98
949
143
reproductive-mode
reproductive-mode
"lottery-polygyny" "hermaphrodite"
0

SLIDER
647
526
776
559
reg-mu-scaling
reg-mu-scaling
0
3
0.0
1
1
NIL
HORIZONTAL

SWITCH
508
10
634
43
N-constant?
N-constant?
1
1
-1000

CHOOSER
780
668
936
713
allele-model
allele-model
"infinite" "diallelic"
0

PLOT
216
631
397
773
geno/phenotypic variance
NIL
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"VP" 1.0 0 -16777216 true "" ""
"VG" 1.0 0 -13791810 true "" ""

SWITCH
777
465
936
498
time-lag-dev?
time-lag-dev?
1
1
-1000

SLIDER
780
634
936
667
genetic-bias
genetic-bias
1
number-of-loci - 1
2.0
1
1
NIL
HORIZONTAL

TEXTBOX
782
615
892
633
if epigenetics:
12
53.0
1

BUTTON
793
311
881
344
crt-map
crt-map\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
838
390
943
423
stop-recording
vid:save-recording \"PanModel33.mp4\"\nprint vid:recorder-status\nstop
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
838
356
943
389
crt-video
crt-video
1
1
-1000

@#$#@#$#@
## WHAT IS IT?

This model study the ability of a sexual population to adapt to its local environment. The population can adapt by means of genetic changes and phenotypic plasticity. In addition, there are many scenarios for the simulation of environmental change (e.g., climate change, seasonal changes)

## HOW IT WORKS

Set the parameter values defining the
Ecology: characteristics of the environment, carrying capacity, type of organisms (r or k strategist)

Evolution: genetic properties and phenotypic plasticity

Then, set the time that determines for how many generations the model will be run

Setup and go to run the model

The figure panels will follow in real time the eco-evolutionary dynamic. 

## HOW TO USE IT

Please refer to the manual and documentation located in the githup repository (link below)

## THINGS TO NOTICE

There are two different approaches to simulate the underlying genetics of the population. Implicit and explicit genetics. By default the model uses the explicit genetics method.

## THINGS TO TRY

It is recommended to start with the simplest as possible scenario (e.g., no stochasticity in environmental fluctuations, and no plasticity, i.e., plasticity = standard-model). Then, the complexity can be progressively increased.

For example, one could try testing the consequences of changing genetic parameters of the model. Then, the model complexity can be increased to evaluate, for instance, the effect of environmental stochasticity on the persistence of the population

## EXTENDING THE MODEL

It would be interesting to extend the model to account for environmental heterogeneity. The code already provides hints on where to add the new lines in this respect

## NETLOGO FEATURES

The model uses native netlogo functions only

## RELATED MODELS

There is another simpler version of the model available on GitHub

## CREDITS AND REFERENCES

For more information, please visit:
https://github.com/danielrm84/PanModel33
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.3.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="evolution_plasticity_periodic_diploid_revision" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>clear-all
reset-ticks
setup-landscape
setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean     [sensibility-factor] of turtles</metric>
    <metric>mean     [neutral-marker] of turtles</metric>
    <metric>mean     [fitness] of turtles</metric>
    <metric>variance [fitness] of turtles</metric>
    <metric>mean     [phenotype] of turtles</metric>
    <metric>variance [phenotype] of turtles</metric>
    <metric>mean     [genetic-component] of turtles</metric>
    <metric>variance [genetic-component] of turtles</metric>
    <enumeratedValueSet variable="plot-output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;cyclic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="amplitude">
      <value value="1"/>
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ROUGHNESS">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="allele-model">
      <value value="&quot;infinite&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reproductive-mode">
      <value value="&quot;lottery-polygyny&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="period">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-type?">
      <value value="&quot;fractal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-constant?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haploid?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file-name">
      <value value="&quot;map.png&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="slope">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="asexual?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neutral-marker?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-mu-scaling">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-model?">
      <value value="&quot;non-spatial&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-HABITATS">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-scaling-factor">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="map-scaling-factor">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loci-reg-genome">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-dist?">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Evolution?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-plasticity?">
      <value value="&quot;epigenetics&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-lag-dev?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-bias">
      <value value="5"/>
      <value value="7"/>
      <value value="9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="evolution_plasticity_uniform_diploid_revision" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>clear-all
reset-ticks
setup-landscape
setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean     [sensibility-factor] of turtles</metric>
    <metric>median   [sensibility-factor] of turtles</metric>
    <metric>variance [sensibility-factor] of turtles</metric>
    <metric>mean     [neutral-marker] of turtles</metric>
    <metric>median   [neutral-marker] of turtles</metric>
    <metric>variance [neutral-marker] of turtles</metric>
    <metric>mean     [fitness] of turtles</metric>
    <metric>variance [fitness] of turtles</metric>
    <metric>mean     [phenotype] of turtles</metric>
    <metric>variance [phenotype] of turtles</metric>
    <metric>mean     [genetic-component] of turtles</metric>
    <metric>variance [genetic-component] of turtles</metric>
    <metric>mean     reduce sentence [reduce sentence (list dna-strain1 dna-strain2)] of turtles</metric>
    <metric>variance reduce sentence [reduce sentence (list dna-strain1 dna-strain2)] of turtles</metric>
    <enumeratedValueSet variable="plot-output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;climate-change&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="amplitude">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ROUGHNESS">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="allele-model">
      <value value="&quot;infinite&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reproductive-mode">
      <value value="&quot;lottery-polygyny&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="period">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-type?">
      <value value="&quot;fractal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-constant?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haploid?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file-name">
      <value value="&quot;map.png&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="slope">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="asexual?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neutral-marker?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-mu-scaling">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-model?">
      <value value="&quot;non-spatial&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-HABITATS">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-scaling-factor">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="map-scaling-factor">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loci-reg-genome">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-dist?">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Evolution?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-plasticity?">
      <value value="&quot;epigenetics&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-lag-dev?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="evolution_plasticity_periodic_diploid_T1_3_10_revision" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>clear-all
reset-ticks
setup-landscape
setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean     [sensibility-factor] of turtles</metric>
    <metric>mean     [neutral-marker] of turtles</metric>
    <metric>mean     [fitness] of turtles</metric>
    <metric>variance [fitness] of turtles</metric>
    <metric>mean     [phenotype] of turtles</metric>
    <metric>variance [phenotype] of turtles</metric>
    <metric>mean     [genetic-component] of turtles</metric>
    <metric>variance [genetic-component] of turtles</metric>
    <enumeratedValueSet variable="plot-output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;cyclic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="amplitude">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ROUGHNESS">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="allele-model">
      <value value="&quot;infinite&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reproductive-mode">
      <value value="&quot;lottery-polygyny&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="period">
      <value value="1"/>
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-type?">
      <value value="&quot;fractal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-constant?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haploid?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file-name">
      <value value="&quot;map.png&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="slope">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="asexual?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neutral-marker?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-mu-scaling">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-model?">
      <value value="&quot;non-spatial&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-HABITATS">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-scaling-factor">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="map-scaling-factor">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loci-reg-genome">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-dist?">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Evolution?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-plasticity?">
      <value value="&quot;epigenetics&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-lag-dev?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-bias">
      <value value="5"/>
      <value value="7"/>
      <value value="9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="evolution_plasticity_time_series" repetitions="30" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>clear-all
reset-ticks
setup-landscape
setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean     [sensibility-factor] of turtles</metric>
    <metric>median   [sensibility-factor] of turtles</metric>
    <metric>variance [sensibility-factor] of turtles</metric>
    <metric>mean     [neutral-marker] of turtles</metric>
    <metric>median   [neutral-marker] of turtles</metric>
    <metric>variance [neutral-marker] of turtles</metric>
    <enumeratedValueSet variable="plot-output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;cyclic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost">
      <value value="0"/>
      <value value="0.5"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="amplitude">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ROUGHNESS">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="allele-model">
      <value value="&quot;infinite&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reproductive-mode">
      <value value="&quot;lottery-polygyny&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="period">
      <value value="1"/>
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-type?">
      <value value="&quot;fractal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-constant?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haploid?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file-name">
      <value value="&quot;map.png&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="slope">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="asexual?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neutral-marker?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-mu-scaling">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-model?">
      <value value="&quot;non-spatial&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-HABITATS">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-scaling-factor">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="map-scaling-factor">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loci-reg-genome">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-dist?">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Evolution?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-plasticity?">
      <value value="&quot;epigenetics&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-lag-dev?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="evolution_plasticity_uniform_diploid_bias_revision" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>clear-all
reset-ticks
setup-landscape
setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean     [sensibility-factor] of turtles</metric>
    <metric>mean     [neutral-marker] of turtles</metric>
    <metric>mean     [fitness] of turtles</metric>
    <metric>variance [fitness] of turtles</metric>
    <metric>mean     [phenotype] of turtles</metric>
    <metric>variance [phenotype] of turtles</metric>
    <metric>mean     [genetic-component] of turtles</metric>
    <metric>variance [genetic-component] of turtles</metric>
    <enumeratedValueSet variable="plot-output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;climate-change&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="amplitude">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ROUGHNESS">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="allele-model">
      <value value="&quot;infinite&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reproductive-mode">
      <value value="&quot;lottery-polygyny&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="period">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-type?">
      <value value="&quot;fractal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-constant?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haploid?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file-name">
      <value value="&quot;map.png&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="slope">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="asexual?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neutral-marker?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-mu-scaling">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-model?">
      <value value="&quot;non-spatial&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-HABITATS">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-scaling-factor">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="map-scaling-factor">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loci-reg-genome">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-dist?">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Evolution?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-plasticity?">
      <value value="&quot;epigenetics&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-lag-dev?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-bias">
      <value value="5"/>
      <value value="7"/>
      <value value="9"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="evolution_plasticity_periodic_diploid_revision_extreme" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>clear-all
reset-ticks
setup-landscape
setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean     [sensibility-factor] of turtles</metric>
    <metric>mean     [neutral-marker] of turtles</metric>
    <metric>mean     [fitness] of turtles</metric>
    <metric>variance [fitness] of turtles</metric>
    <metric>mean     [phenotype] of turtles</metric>
    <metric>variance [phenotype] of turtles</metric>
    <metric>mean     [genetic-component] of turtles</metric>
    <metric>variance [genetic-component] of turtles</metric>
    <enumeratedValueSet variable="plot-output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;cyclic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="amplitude">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ROUGHNESS">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="allele-model">
      <value value="&quot;infinite&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reproductive-mode">
      <value value="&quot;lottery-polygyny&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="period">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-type?">
      <value value="&quot;fractal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-constant?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haploid?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file-name">
      <value value="&quot;map.png&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="slope">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="asexual?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neutral-marker?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-mu-scaling">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-model?">
      <value value="&quot;non-spatial&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-HABITATS">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-scaling-factor">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="map-scaling-factor">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loci-reg-genome">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-dist?">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Evolution?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-plasticity?">
      <value value="&quot;epigenetics&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-lag-dev?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-bias">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="evolution_plasticity_periodic_diploid_T1_3_10_revision_extreme" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>clear-all
reset-ticks
setup-landscape
setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean     [sensibility-factor] of turtles</metric>
    <metric>mean     [neutral-marker] of turtles</metric>
    <metric>mean     [fitness] of turtles</metric>
    <metric>variance [fitness] of turtles</metric>
    <metric>mean     [phenotype] of turtles</metric>
    <metric>variance [phenotype] of turtles</metric>
    <metric>mean     [genetic-component] of turtles</metric>
    <metric>variance [genetic-component] of turtles</metric>
    <enumeratedValueSet variable="plot-output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;cyclic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="amplitude">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ROUGHNESS">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="allele-model">
      <value value="&quot;infinite&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reproductive-mode">
      <value value="&quot;lottery-polygyny&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="period">
      <value value="1"/>
      <value value="3"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-type?">
      <value value="&quot;fractal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-constant?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haploid?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file-name">
      <value value="&quot;map.png&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="slope">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="asexual?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neutral-marker?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-mu-scaling">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-model?">
      <value value="&quot;non-spatial&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-HABITATS">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-scaling-factor">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="map-scaling-factor">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loci-reg-genome">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-dist?">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Evolution?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-plasticity?">
      <value value="&quot;epigenetics&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-lag-dev?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-bias">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="evolution_plasticity_uniform_diploid_revision_extreme" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>clear-all
reset-ticks
setup-landscape
setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean     [sensibility-factor] of turtles</metric>
    <metric>mean     [neutral-marker] of turtles</metric>
    <metric>mean     [fitness] of turtles</metric>
    <metric>variance [fitness] of turtles</metric>
    <metric>mean     [phenotype] of turtles</metric>
    <metric>variance [phenotype] of turtles</metric>
    <metric>mean     [genetic-component] of turtles</metric>
    <metric>variance [genetic-component] of turtles</metric>
    <enumeratedValueSet variable="plot-output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;climate-change&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="amplitude">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ROUGHNESS">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="allele-model">
      <value value="&quot;infinite&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reproductive-mode">
      <value value="&quot;lottery-polygyny&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="period">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-type?">
      <value value="&quot;fractal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0.004"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-constant?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haploid?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file-name">
      <value value="&quot;map.png&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="slope">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="asexual?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neutral-marker?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-mu-scaling">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-model?">
      <value value="&quot;non-spatial&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-HABITATS">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-scaling-factor">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="map-scaling-factor">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loci-reg-genome">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-dist?">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Evolution?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-plasticity?">
      <value value="&quot;epigenetics&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-lag-dev?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-bias">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="evolution_plasticity_uniform_haploid_revision" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>clear-all
reset-ticks
setup-landscape
setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean     [sensibility-factor] of turtles</metric>
    <metric>mean     [neutral-marker] of turtles</metric>
    <metric>mean     [fitness] of turtles</metric>
    <metric>variance [fitness] of turtles</metric>
    <metric>mean     [phenotype] of turtles</metric>
    <metric>variance [phenotype] of turtles</metric>
    <metric>mean     [genetic-component] of turtles</metric>
    <metric>variance [genetic-component] of turtles</metric>
    <enumeratedValueSet variable="plot-output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;climate-change&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="amplitude">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ROUGHNESS">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="allele-model">
      <value value="&quot;infinite&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reproductive-mode">
      <value value="&quot;lottery-polygyny&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="period">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-type?">
      <value value="&quot;fractal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0.001"/>
      <value value="0.004"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-constant?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haploid?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file-name">
      <value value="&quot;map.png&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="slope">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="asexual?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neutral-marker?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-mu-scaling">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-model?">
      <value value="&quot;non-spatial&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-HABITATS">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-scaling-factor">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="map-scaling-factor">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loci-reg-genome">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-dist?">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Evolution?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-plasticity?">
      <value value="&quot;epigenetics&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-lag-dev?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-bias">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="evolution_plasticity_periodic_haploid_revision" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>clear-all
reset-ticks
setup-landscape
setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean     [sensibility-factor] of turtles</metric>
    <metric>mean     [neutral-marker] of turtles</metric>
    <metric>mean     [fitness] of turtles</metric>
    <metric>variance [fitness] of turtles</metric>
    <metric>mean     [phenotype] of turtles</metric>
    <metric>variance [phenotype] of turtles</metric>
    <metric>mean     [genetic-component] of turtles</metric>
    <metric>variance [genetic-component] of turtles</metric>
    <enumeratedValueSet variable="plot-output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;cyclic&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="amplitude">
      <value value="1"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ROUGHNESS">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="allele-model">
      <value value="&quot;infinite&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reproductive-mode">
      <value value="&quot;lottery-polygyny&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="period">
      <value value="1"/>
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-type?">
      <value value="&quot;fractal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-constant?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haploid?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file-name">
      <value value="&quot;map.png&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="slope">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="asexual?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neutral-marker?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-mu-scaling">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-model?">
      <value value="&quot;non-spatial&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-HABITATS">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-scaling-factor">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="map-scaling-factor">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loci-reg-genome">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-dist?">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Evolution?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-plasticity?">
      <value value="&quot;epigenetics&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-lag-dev?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-bias">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="evolution_plasticity_uniform_diploid_revision_timeseries" repetitions="30" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>clear-all
reset-ticks
setup-landscape
setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>mean     [sensibility-factor] of turtles</metric>
    <enumeratedValueSet variable="plot-output?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scenario?">
      <value value="&quot;climate-change&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-rate-per-locus">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mut-effect-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-organism">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cost">
      <value value="0"/>
      <value value="0.5"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-capacity">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fitness-function">
      <value value="&quot;negative-exponential&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="amplitude">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ROUGHNESS">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="allele-model">
      <value value="&quot;infinite&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="level-autocorr">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-dependence-effect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reproductive-mode">
      <value value="&quot;lottery-polygyny&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="period">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="landscape-type?">
      <value value="&quot;fractal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-env-variance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rate-change-of-optimum">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-constant?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="haploid?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-to-balance">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="file-name">
      <value value="&quot;map.png&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="slope">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="asexual?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mutational-variance">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neutral-marker?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-loci">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="heritability">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-mu-scaling">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="type-of-model?">
      <value value="&quot;non-spatial&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-HABITATS">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="reg-scaling-factor">
      <value value="42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetics?">
      <value value="&quot;explicit&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-genetic-variance">
      <value value="&quot;parental-level&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="map-scaling-factor">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loci-reg-genome">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-variance">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-dist?">
      <value value="&quot;normal&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-size">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-limit">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Evolution?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="how-plasticity?">
      <value value="&quot;epigenetics&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-lag-dev?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="genetic-bias">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
