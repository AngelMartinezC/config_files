# SOLAR ARRAYS

- <h3>SM0:</h3>
 
 Solar model folder `x=100 y=200`


- <h3>SM1:</h3>

 Solar model folder `x=100 y=200` but for new depth at half

 - **Attempt 2:** I have to correct radtopix function because has not a true limit.


- <h3>SM2:</h3>
 
 Solar model folder `x=150 y=300` as original 
 

</br></br>


----

# SIMULATION SCHEMES


<!--# solar\_model\_0:
+ `X2BEG` conditions userdef.
# solar\_model\_1:
* `X2BEG` conditions reflective: 
# _solar\_model\_2:_
* If variables are not set, then are `False`. -->

* <span style="color:red">``side == 0``:</span> With `INTERNAL_BOUNDARY` flag enabled, there is no density change along the array.

## No magnetic field:

- <span style="color:blue">_SCHEME A_:</span> No magnetic field tube. $B_e=0.1$.
    
  |   |  <span style="color:green">`X2_BEG`</span> | <span style="color:green">`X2_END`</span>  |  <span style="color:violet">`InitDomain`</span> |
  |:-:|:-:|:-:|:-:|
  | $\rho$   | $\rho_e$  |  $=$ | $=$ | 
  | $P$      |  $P_{e}$ | $=$  | $=$  |
  | $B_{X2}$ | $B_e$  | $B_e$ | $B_e$ | 
  | $V_{X2}$ | `outflow`   | $-0.05~V_{X2}$   |  |
  
  - From `6 h`, the array shows to have convective-like currents at the bottom of the box.
  
  - From `1 d: 3 h` there is a great horizontal component of velocity with increases at `1 d: 10 h` in [0,1] Mm becoming unstable.
  
- <span style="color:blue">_SCHEME B_:</span> No magnetic field tube. $B_e=X2[j]/1\times 10^3 + 0.1$.
    
  |   |  <span style="color:green">`X2_BEG`</span> | <span style="color:green">`X2_END`</span>  |  <span style="color:violet">`InitDomain`</span> |
  |:-:|:-:|:-:|:-:|
  | $\rho$   | $\rho_e$  |  $=$ | $=$ | 
  | $P$      |  $P_{e}$ | $=$  | $=$  |
  | $B_{X2}$ | $B_e$  | $B_e$ | $B_e$ | 
  | $V_{X2}$ | `outflow`   | $-0.05~V_{X2}$   |  |

  - From `6 h`, the array shows to have convective-like currents at the bottom of the box.
  
  - Great changes at `21 h` there are upper temporary changes. At `23 h; 40 m` there are great changes. 


  
---







##  Magnetic field:
 
  
- <span style="color:red">_SCHEME 0_</span>: `InitDomain` magnetic field tube at $t=0$.
  
 - $B_e=0.1$
  
 - $B_i=1000$
  
 - $T_f=1.5$ 
  
 - $\rho_e,~P_e:$ Solar Model
  
  |   |  <span style="color:orange">`X2_BEG`</span>  | <span style="color:orange">`X2_END`</span> |  <span style="color:violet">`InitDomain`</span>  |  <span style="color:green">`X2_BEG`</span> | <span style="color:green">`X2_END`</span>  |  <span style="color:violet">`InitDomain`</span> |
  |:-:|:-:|:-:|:-:|:-:|:-:|:-:|
  | $\rho$   | $T_f ~\rho_e$  |  $=$ | $=$ | $\rho_e$ | $=$ | $=$|
  | $P$      |  $P_{e} + \frac{1}{2}\big(B_{e}^2-B_{i}^2 \big)$ | $=$  | $=$  | $P_{e}$ | $=$ | $=$|
  | $B_{X2}$ | $B_i$  | $B_i$ | $B_i$ | $B_e$ | $B_e$ | $B_e$ |
  | $V_{X2}$ | `outflow`   | $-0.05~V_{X2}$   | | `outflow`  | $-0.07~V_{X2}$  | |
  
 Magnetic structure breaks at 21/45 minutes. 
 
 ---
 
 
 
 
- <span style="color:red">_SCHEME 1_</span>: `InitDomain` magnetic field tube at $t=0$.
  
 - $B_e=0.1$
  
 - $B_i=1000$
  
 - $T_f=1.5$ 
  
 - $\rho_e,~P_e:$ Solar Model
  
  |   |  <span style="color:orange">`X2_BEG`</span>  | <span style="color:orange">`X2_END`</span> |  <span style="color:violet">`InitDomain`</span>  |  <span style="color:green">`X2_BEG`</span> | <span style="color:green">`X2_END`</span>  |  <span style="color:violet">`InitDomain`</span> |
  |:-:|:-:|:-:|:-:|:-:|:-:|:-:|
  | $\rho$   | $T_f ~\rho_e$  |  $=$ | $=$ | $\rho_e$ | $=$ | $=$|
  | $P$      |  $P_{e} + \frac{1}{2}\big(B_{e}^2-B_{i}^2 \big)$ | $=$  | $=$  | $P_{e}$ | $=$ | $=$|
  | $B_{X2}$ | $B_i$  | $B_i$ | $B_i$ | <span style="color:red">`outflow`</span> | $B_e$ | $B_e$ |
  | $V_{X2}$ | `outflow`   | $-0.05~V_{X2}$   | | `outflow`  | $-0.07~V_{X2}$  | |
  
 Magnetic structure breaks at 20 minutes. 
  
  
  ---
  
- <span style="color:red">_SCHEME 2_</span>: Scheme following Bgradient as SCHEME B, setting the magnetic field from $t=0$ in `side == 0`.

  `dt` is too small at `7h 22m`. A huge $x-$velocity component is seen to the left of the tube. Apparently I have to set $B_{X2}$ after a "setting time" of `~6 h` according to former schemes.
  
  
  ---
  
- ~~<span style="color:red">_SCHEME 3_</span>: Scheme as SCHEME 2 but with `B_EXT` fixed in `side == 0`~~:

 ```c 
    TOT_LOOP(k,j,i){
     Bext = x2[j]/1.0e3 + Bext0;
     if (x1[i]>-radius && x1[i] < radius){ 
      d->Vc[BX2][k][j][i] = Bint/press_unit;
      d->Vc[RHO][k][j][i] = (T_f)*InputDataInterpolate(id1,x1[i],x2[j],x3[k]);
      d->Vc[PRS][k][j][i] = InputDataInterpolate (id2,x1[i],x2[j],x3[k]) +
        (1.0/(2.0*pow(press_unit,2)))*(pow(Bext,2)-pow(Bint,2));
      }
     d->Vc[BX2][k][j][i] = Bext/press_unit;
    }
 ```
 Bad magnetic field.
    
  ---
    
- <span style="color:red">_SCHEME 4_</span>: Scheme as SCHEME 2  but with else `B_EXT` fixed in `side == 0`

 ```c
    TOT_LOOP(k,j,i){
     Bext = x2[j]/1.0e3 + Bext0;
      if (x1[i]>-radius && x1[i]<radius){ 
       d->Vc[BX2][k][j][i] = Bint/press_unit;
       d->Vc[RHO][k][j][i] = (T_f)*InputDataInterpolate(id1,x1[i],x2[j],x3[k]);
       d->Vc[PRS][k][j][i] = InputDataInterpolate (id2,x1[i],x2[j],x3[k]) +
        (1.0/(2.0*pow(press_unit,2)))*(pow(Bext,2)-pow(Bint,2));
     }
     else {
      d->Vc[BX2][k][j][i] = Bext/press_unit;
     }
    }
    
 ```
 Great changes near the magnetic field. Small `dt ~ 1e-3` at 20300 s. 
 
  ---

	
- <span style="color:red">_SCHEME 5_</span>: SCHEME 0  but with else `B_EXT = 0.1` fixed in `side == 0`

 Magnetic field is going to be set at frame 348 (`5h 48 min`) where the convection is begining to be significant.
 
  1. ~~As "initial" condition (setting all the magnetic field struture suddenly in `1 min`) setting boundary conditions for 3 minutes~~: The magnetic field suppresed by the environment. 
 
  2. ~~As "initial" condition (just temporary magnetic field)~~: Sudden loss of structure.

  3. ~~From `5 h`~~: For the setting time of mangetic field (`5 min`), all is good, then the structure vanishes as there is some greater outside total pressure.
 
  ---

	
- <span style="color:red">_SCHEME 6_</span>: Two magnetic field lines.
   
  ---





 </br>

---


</br>

# <span style="color:green">  PLUTO V4.4 </span>

</br> </br>
	
- <span style="color:blue">_SCHEME 7_</span>: `solar_model_6 folder` of new PLUTO version. Reading model from SM6.

  First reading the gravity force in runtime.
  
  - What was I thinking! I can save precious run time setting a constant value in boundary instead of reading and interpolating the solar model file. I can store boundary values in an array. :(.
    
     - I just did. Running time improves by more than $33\%$. The time for `~1d 4h` is `1h 24min`.
     
  ---
  
- ~~<span style="color:red">_SCHEME 8_</span>:~~ `solar_model_8 folder` of new PLUTO version. Same as SCHEME 7, but with magnetic field.

  -  `movie_8` improving pressure for `side==0` does not improve too much. There is the same small dt.
  
  - Improve gravity by reading to a constant value but the hole array.
  
  - Make with MHD continua as Nakariakov, then test the horizontal condition for magnetic field

  ---
      
- <span style="color:blue">_**SCHEME 9, 11**_</span>: `solar_model_9 folder` of new PLUTO version. No magnetic field but boundary conditions on both X-walls (either periodic of reflective, the former seems to be better).

  - Seems to behave propertly. At least the huge velocities in x-direction are not noticeable and convective-like currents seem to be present. At `~12 days` of code time, the convective currents are still present.
   
 ---

- <span style="color:red">_SCHEME 10_</span>: `solar_model_10 folder`. The mgnetic field is set at `InitDomain`. I should apply this "initial condition" on a stable time.

 ---

- <span style="color:red">_SCHEME 12_</span>: `solar_model_12 folder`. Same as SCHEME 11 with outflow conditon for both X walls.

 ---

- <span style="color:purple">_SCHEME 13_</span>: `solar_model_13 folder`. I am going to set a background magnetic field for the same scheme of scheme 11 with straggled components of the magnetic field.


 ---

- <span style="color:purple">_SCHEME 14_</span>: `solar_model_10 folder`. The mgnetic field is set at `InitDomain`. I should apply this "initial condition" on a stable time.

 ---

- <span style="color:red">_SCHEME 15_</span>: Smaller grid.
 
  - I should use the bigger grid possible

 ---

 ---

 CHANGE PLUTO ENTRIES

 ---
 
 ---

- <span style="color:purple">_SCHEME 16_</span>: Test the setting of the magnetic field. I am going to try several schemes to get the most stable tube. The start SCHEME is 9, from about 5 days where there is a window of almost steady convective currents (velocity vector field).


- <span style="color:red">_SCHEME 17_</span>: Here I set a great InputDataFile in order to perform a smooth interpolation to be read on Init function. The y_range is set manually because of dynamic variables to 302~304 (in order to preserve ghosts zones). This is based on a new `solar_model_test_2.3` which, with the same grid, is seen to have a stability window beyond `0d 20h` up to `2h`.

  - Setting a mangetic field from stability (from a down convective cell, no matter the width) crashes the program.

  - Removing magnetic field but mantaining boundary magneticnetic tube conditions does not improve results, ending up with great velocities in x-direction.


- <span style="color:red">_SCHEME 18_</span>: As SCHEME 17 but with different velocity reflection on the `X2_END` $v_y$ component (reflection at 80%).

  - It shows stability. It's like the former one but stretched.

</br> </br> </br> </br> </br> </br> 


<!-- blank line --> 
