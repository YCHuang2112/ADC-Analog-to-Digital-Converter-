* IDC_mod simulation
* for modulator simulation only 
.option runlvl=6
* simulation option
.options post_version=9601 
.options method=gear lvltim=2 vntol=1e-6 reltol=1e-6 rmin=1e-15 
+ delmax=1n 
+ post probe INTERP 
.options post


**** netlist
.protect
.lib 'rf018.l' tt
.unprotect
.param nnoiseflag=1 pnoiseflag=1
.inc 'CMP5_mis.cir'




.param pwra=1.8 pwrd=1.8
* power supply
vavdd  AVDD  0 dc pwra
vagnd  AGND  0 dc 0
vavdd2  AVDD2  0 dc pwra
vsub   SUB   0 dc 0


.param fclk=100e6 tp='1/fclk' td1=0.1n td2=0.1n tr=0.1n tf=0.1n 
vclk  CLK 0 pulse ( 0 pwra '0.6*tp' tr tf 'tp/2-td1-td2-tr' tp ) 
.param cycle=10

* signal source
evcm1 vcm1 0 vcm 0 1
*.param lsb='pwra/2^10'
.param vp='1.8'

.param nlsb=4 v_delta='2*vp/(2**10)' n_part='(2**10-nlsb)/2' n_parts_in_lsb="100*cycle" ramp_time="(nlsb*n_parts_in_lsb+1)*tp"
v_inp INP vcm1 pulse('-vp+v_delta*n_part' 'vp-v_delta*n_part'  '0.6*tp'  'ramp_time-tp' tf  'ramp_time' 'ramp_time*2')
v_inn INN vcm1 pulse('vp-v_delta*n_part' '-vp+v_delta*n_part'  '0.6*tp'  'ramp_time-tp' tf  'ramp_time' 'ramp_time*2')

v_refp  VREFP 0  'pwra'
v_refn  VREFN 0  '0'
v_vcm   VCM   0  'pwra/2'

* temperature
.temp 60 

* analysis
evcvs_cmp ideal_cmp_out 0 sp1 sn1  1000  MAX=pwra MIN=0
.tran step='tp' stop='ramp_time' 
.probe v(clk) v(inp) v(inn) v(stage1_p) v(stage1_n) v(outp) v(outn) v(ideal_cmp_out)
.meas tran pwr avg power
.trannoise v(outp,outn) method=mc samples=50 seed=2 autocorrelation=0 fmin=1 fmax='100/tp' scale=1.0

favdd meas_avdd 0 vavdd 1
cavdd meas_avdd 0 1u
favdd2 meas_avdd2 0 vavdd2 1
cavdd2 meas_avdd2 0 1u

.meas tran t2_pwr_avdd find par('v(meas_avdd)*pwra*1e-6/ramp_time') at='0*tp+ramp_time'
.meas tran t1_pwr_avdd find par('v(meas_avdd)*pwra*1e-6/ramp_time') at='0*tp'
.meas tran pwr_avdd  param="t2_pwr_avdd-t1_pwr_avdd"
.meas tran t2_pwr_avdd2 find par('v(meas_avdd2)*pwra*1e-6/ramp_time') at='0*tp+ramp_time'
.meas tran t1_pwr_avdd2 find par('v(meas_avdd2)*pwra*1e-6/ramp_time') at='0*tp'
.meas tran pwr_avdd2  param="t2_pwr_avdd2-t1_pwr_avdd2"

**** probing
*.probe vin=v(vinp,vinn)  

.end
