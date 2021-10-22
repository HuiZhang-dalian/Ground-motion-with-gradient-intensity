proc PeakResponse {period Ry} {
	wipe
	wipeAnalysis
	global dt
	global nPts
	global record
	global Fact
	global dataDir
	global F	
	model BasicBuilder -ndm 1 -ndf 1
	node 0 0. 
	fix 0 1 
	node 1 0.
	mass 1 1. 
	set stiffness [expr pow(2 * 3.141593 / $period, 2)]
	# uniaxialMaterial Elastic $matTag $E <$eta> || $eta:  damping tangent (optional, default=0.0)
	uniaxialMaterial Elastic 1 $stiffness 
	# element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2 ... <-orient $x1 $x2 $x3 $yp1 $yp2 $yp3>
	element zeroLength 1 0 1 -mat 1 -dir 1 
	# 这里没有找到对应的命令，而这个命令很相似, Series -dt $dt -filePath $fileName <-factor $cFactor> || while the values are taken from a file specified by $fileName
	timeSeries Path 1 -dt $dt -filePath $record
	set dampingRatio 0.05
	set alpha [expr 4*3.141593 / $period * $dampingRatio]
	rayleigh $alpha 0 0 0
	# pattern UniformExcitation $patternTag $dir -accel (TimeSeriesType arguments)<-vel0 $ver0> 
	# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	# $patternTag 			unique pattern object tag
	# $dir 							direction of excitation (1, 2 or 3) used in formulating the inertial loads for the transient analysis
	#TimeSeriesType arguments	TimeSeries object associated with the acceleration record used in determining the inertial loads
	#$vel0							initial velocity to be assigned to each node (optional, default: zero)
	pattern UniformExcitation 1 1 -accel 1 -fact 1
	# this command is used to construct a TransformationConstrainHandler which will cause the constraints to be enforced using the transformation method
	constraints Transformation
	# This command is used to construct a RCMUnumberer object; The RCM numberer uses the reverse Cuthill McKee (REF?) algorithm to number th degrees of freedom.
	numberer RCM
	# This command is used to construct a general sparese system of equations object which will be factored and solved during the analysis using the UMFPACK solver. (REF?)
	system UmfPack
	# This command is used to construct a CTestEnergyIncr object which tests positive force convergence if one half of the inner-product 
	# of the x and b vectors (displacement and increment and unbalance) in the LinearSOE object less than the sspecified tolerance.
	# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	# test EnergyIncr $tol $maxNumIter <$printFlag>
	# $tol								convergence tolerance
	# $maxNumIter			maximum number of iterations that will be performed before
	# $printFlag					flag used to print information on convergence (optional):
	#											0: no print output (default);	
	#											1: print infromation on each step;
	#											2: print information when convergence has been achieved;	
	# 										4: print norm, dU and dR vectors;	
	#											5: at convergence failure, carry on, print error message, but do not stop analysis.
	test EnergyIncr 1 200
	# This command is used to construct a Newton Raphson algorithm object which used the Newton-Raphson method to advance to the next time step.
	# Note: The tangent is updated at each iteration.
	algorithm Newton
	# This command is used to construct a TransientIntegrator object of type Newmark
	# integrator Newmark $gamma $beta
	# $gamma				Newmark parameter γ
	# #beta 					Newmark parameter β
	# The damping matrix D is specified as a combination of stiffness and mass-proportional damping matrices:
	# D = $alphaM * M + $betaK * Kcurrent +$betaKinit * Kinit + $betaKcomm * KlastCommit 
	# The mass and stiffness matrices are defined as: 
	# M 								mass matrix used to calculate Rayleigh Damping 
	# Kcurrent 					stiffness matrix at current state determination used to calculate Rayleigh Damping 
	# Kinit 							stiffness matrix at initial state determination used to calculate Rayleigh Damping 
	# KlastCommit 		stiffness matrix at last-committed state determination used to calculate Rayleigh Damping 
	# Note: Damping should be specified using the Rayleigh (page 342) command. 
	integrator Newmark 0.5 0.25
	# This command is used to construct a DirectIntegrationAnalysis object. 
	# This analysis object is constructed with the component objects previously created by the analyst. If none has been created, default objects are constructed and used: 
	#Component Object Default object 
	# SolutionAlgorithm (page 329), TransientIntegrator (page 333):
	# 	NewtonRaphson (page 329) EquiSolnAlgo with a CTestNormUnbalance (page 326) with a tolerance of 1e-6 and a maximum of 25 iterations 
	# ConstraintHandler (page 316) :
	# 	PlainHandler (page 318) ConstraintHandler 
	# DOF_Numberer (page 321) :
	# 	RCM (page 322) DOF_Numberer 
	# LinearSOE (page 323), LinearSolver (page 323) :
	# 	profiled symmetric positive definite  (page 324) LinearSOE and Linear Solver 
	# Chapter 36    analysis Command 341 :
 	# 	Integrator Newmark TransientIntegrator (page 337) with γ=0.5 and
	analysis Transient
	set maxAccel 0.0
	set maxDisp 0.0 
	set maxVel 0.0
	for {set i 0} {$i < $nPts} {incr i} {
		analyze 1 $dt
 		if {$maxAccel < [expr abs([lindex [eleResponse 1 forces] 0])]} {
 			set maxAccel [expr abs([lindex [eleResponse 1 forces] 0])]
 		}
 		if {$maxDisp < [expr abs([lindex [nodeDisp 1]])]} {
 			set maxDisp [expr abs([lindex [nodeDisp 1]])]
 		}
 		if {$maxVel < [expr abs([lindex [nodeVel 1]])]} {
 			set maxVel [expr abs([lindex [nodeVel 1]])]
 		}
	}
	wipe
	wipeAnalysis
	model BasicBuilder -ndm 1 -ndf 1
	node 0 0.
	fix 0 1 
	node 1 0.
	mass 1 1. 
	set stiffness [expr pow(2 * 3.141593 / $period, 2)]
	set moe [expr 1.0/$Ry]
	set s1p $moe;set e1p [expr $s1p/$stiffness];set e2p [expr 4*$e1p];set s2p [expr $s1p*1.06];set s1n -$s1p;set e1n -$e1p;set s2n -$s2p;set e2n -$e2p;
	set e3p [expr (4+1.06/0.3)*$e1p];set s3p 1e-10;set e3n -$e3p;set s3n -$s3p;set pinchX 1;set pinchY 1; set damage1 0;set damage2 0.1;set beta 0.5;
	# This command is used to construct a uniaxial biliear hysteretic material object with pinching of force and deformation, damage due to ductility and energy,
	# and degraded unloading stiffness based on ductility
	# uniaxialMaterial Hysteretic $matTag $s1p $e1p $s2p $e2p <$s3p $e3p> $s1n $e1n $s2n $e2n <$s3n $e3n> $pinchX $pinchY $damage1 $damage2 <$beta>
	# $matTag unique material object integer tag 
	# $s1p $e1p stress and strain (or force & deformation) at first point of the envelope in the positive direction 
	# $s2p $e2p stress and strain (or force & deformation) at second point of the envelope in the positive direction 
	# $s3p $e3p stress and strain (or force & deformation) at third point of the envelope in the positive direction (optional) 
	# $s1n $e1n stress and strain (or force & deformation) at first point of the envelope in the negative direction* 
	# $s2n $e2n stress and strain (or force & deformation) at second point of the envelope in the negative direction* 
	# $s3n $e3n stress and strain (or force & deformation) at third point of the envelope in the negative direction (optional)* 
	# $pinchX pinching factor for strain (or deformation) during reloading $pinchY pinching factor for stress (or force) during reloading 
	# $damage1 damage due to ductility: D1(mu-1) $damage2 damage due to energy: D2(Eii/Eult) 
	# $beta power used to determine the degraded unloading stiffness based on ductility, mu-beta (optional, default=0.0) 
	# *NOTE: negative backbone points should be entered as negative numeric values 
	
	uniaxialMaterial Hysteretic 1 $s1p $e1p $s2p $e2p $s3p $e3p $s1n $e1n $s2n $e2n $s3n $e3n $pinchX $pinchY $damage1 $damage2 $beta
	#uniaxialMaterial Hysteretic 1 $s1p $e1p $s2p $e2p $s1n $e1n $s2n $e2n $pinchX $pinchY $damage1 $damage2 $beta
	# uniaxialMaterial Steel01 1 $moe $stiffness 0.02
	
	# This command is used to construct a zeroLength element object, which is defined by two nodes at the same location. 
	# The nodes are connected by multiple UniaxialMaterial (page 43) objects to represent the force-deformation relationship for the element. 
	# element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2 ... <-orient $x1 $x2 $x3 $yp1 $yp2 $yp3> 
	# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	# $eleTag 										unique element object tag $iNode  $jNode end nodes 
	# $matTag1 $matTag2 ... 		tags associated with  previously-defined UniaxialMaterials (page 43) 
	# $dir1 $dir2 ...  							material directions:   1,2,3  translation along local x,y,z axes, respectively   4,5,6 rotation about local x,y,z axes, respectively 
	# 
	# the orientation vectors can be specified for the element (optional): 
	# $x1 $x2 $x3 vector components in global coordinates defining local x-axis (vector x) 
	# $yp1 $yp2 $yp3 						vector components in global coordinates defining vector yp which lies in the local x-y plane for the element:   the local z-axis is defined by the cross product between the vectors x and yp  
	# If the optional orientation vectors are not specified, the local element axes coincide with the global axes. 
	# The valid queries to a zero-length element when creating an ElementRecorder (page 307) object are 'force,' 'deformation,' 'stiff,' and 'material $matNum matArg1 matArg2 ...' Where $matNum is the tag associated with the material whose data is to be output. 
	element zeroLength 1 0 1 -mat 1 -dir 1 
	timeSeries Path 1 -dt $dt -filePath $record
	set dampingRatio 0.05
	set alpha [expr 4*3.141593 / $period * $dampingRatio]
	rayleigh $alpha 0 0 0
	pattern UniformExcitation 1 1 -accel 1 -fact  1
	constraints Transformation
	numberer RCM
	system UmfPack
	test EnergyIncr 1 200  
	algorithm Newton
	integrator Newmark 0.5 0.25
	analysis Transient
	set maxAccel2 0.0
	set maxDisp2 0.0 
	set maxVel2 0.0
	#recorder Node -file $dataDir/$record$Ry$period-accel.txt -time -node 1 -dof 1  accel;			# displacements of free node
	#recorder Node -file $dataDir/$record$Ry$period -node 1 -dof 1  reaction;
	#recorder Node -file $dataDir/$record$Ry$period-vel.txt -node 1 -dof 1  vel;
	#recorder Node -file $dataDir/$record$Ry$period-disp.txt -node 1 -dof 1  disp;
	
	# recorder Node -file $dataDir/accel.txt -time -node 1 -dof 1  accel;			# displacements of free node
	# recorder Node -file $dataDir/reaction.txt -node 1 -dof 1  reaction;
	# recorder Node -file $dataDir/vel.txt -node 1 -dof 1  vel;
	# recorder Node -file $dataDir/disp.txt -node 1 -dof 1  disp;
	# analyze $nPts $dt
	
	
	# This part is added by ZZF to induce the accelergration data, in order to compare the value of the absolut acceleration and acceleration and the eleResponse.
	set nop [open $record r]
	set data [read $nop]
	set maxAbsAccel 0.0
	# This part is finised
	set D0 0
	set F0 0
	set E0 0
	set EH0 0
	set EC0 0
	set EK0 0
	set Ein0 0
	set VV0 0
	for {set i 0} {$i < $nPts} {incr i} {
		analyze 1 $dt
		set D [nodeDisp 1]
		set F [nodeReaction 1] 
		set VV [ nodeVel 1]
		# to calculate the hysteretic energy，here is the elastic force work
		 set E [expr (($F+$F0)*($D0-$D)*0.5+$E0)]
		# to calculate the input energy, modified by ZZF
		if {$i< $nPts-1} {
			set Ein [expr (-(([lindex $data $i])+([lindex $data $i+1]))*($D-$D0)*0.5+$Ein0)]
		} else {
			set Ein [expr (-([lindex $data $i]+0)*($D-$D0)*0.5+$Ein0)]
		}
		# here, is the damp energy
		set EC [expr ($alpha*($VV+$VV0)*($D-$D0)*0.5+$EC0)]
		# here, is the dynamic energy
		set EK [expr (0.5*pow($VV,2))]
		set D0 $D
		set F0 $F
		set E0 $E
		set Ein0 $Ein
		set EC0 $EC
		set VV0 $VV
		#set EK0 $EK
		
 		if {$maxAccel2 < [expr abs([lindex [eleResponse 1 forces] 0])]} {
 			set maxAccel2 [expr abs([lindex [eleResponse 1 forces] 0])]
 		}
 		if {$maxDisp2 < [expr abs([lindex [nodeDisp 1]])]} {
 			set maxDisp2 [expr abs([lindex [nodeDisp 1]])]
 		}
 		if {$maxVel2 < [expr abs([lindex [nodeVel 1]])]} {
 			set maxVel2 [expr abs([lindex [nodeVel 1]])]
 		}
 		if {$maxAbsAccel < [expr abs([lindex $data $i]+[lindex [eleResponse 1 forces] 0])]} {
 			set maxAbsAccel [expr abs([lindex $data $i]+[lindex [eleResponse 1 forces] 0])]
 		}
	}
	#puts [expr max(abs([nodeDisp 1]))]
	#puts $maxDisp2
	#puts $EH
	#puts $E0
	#puts $EC
	#puts $EK
	#puts $Ein
	# here, is the hysteristic energy
   	set EH [expr $Ein-$EC-$EK-$E]
	#set E [expr abs($E)]
	set  E [expr abs($E)]
	set  energy [expr $E/$moe/$maxDisp/$Ry]
	set u [expr $maxDisp2/$maxDisp*$Ry ]
	#set maxfile [open $dataDir/$record$Ry-maxvalues.txt w+]
	#puts $maxfile [format "%.2f\t%.8f\t%.8f\t%.8f\n"  $period $maxAccel2 $maxVel2 $maxDisp2]
	stop
	# return [list $maxDisp2 $E]
	# return [list $maxDisp2 $E]
	return [list $E $energy $u $maxAccel2 $maxVel2 $maxDisp2 $maxAbsAccel];
}
source NationCA_ListName_Dt_Dn.tcl
set dataDir model2Data;
file mkdir $dataDir;
for {set i 0} {$i<=1000} {incr i} {
	set record1 [lindex $Name $i];
	# puts $Name
	#set record [string trimright $record1 .bx]
	set record [lindex $Name $i]
	puts $record
	set dt [lindex $Dt $i];
	set nPts [lindex $Npts $i];
	foreach Ry { 10 } {
 	#set spectraFile [open $dataDir/$record$Ry.txt w]
	set energyFile [open $dataDir/$record$Ry.txt w]
	foreach j { 10 } {
		set period [expr $j*0.1]
		set response [PeakResponse $period $Ry]
		# puts $spectraFile [format "%.4f\t%.4f" $period   [expr [lindex $response 1] / [lindex $response 0]] ]
		puts $energyFile [format "%.2f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n" $period  [lindex $response 0] [lindex $response 1] [lindex $response 2] [lindex $response 3] [lindex $response 4] [lindex $response 5] [lindex $response 6] ];
		puts " $j analyze done"
}
# puts "$Ry spectrum done"

# close $spectraFile
puts "$Ry energy done"
close $energyFile
}
}