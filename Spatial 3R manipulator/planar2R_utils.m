(* ::Package:: *)

BeginPackage["planar2Rutils`",{"mathutils`"}]



planar2Rutils::usage = " 'planar2R_utils' package provides simple utilities for finding forward & inverse kinematics, workspace boundaries & singularity condition of a 2-R manipulator.

It contains the following functions: 

1) forwardKinematics2R[{x0,y0},l1,l2,{\[Theta]1,\[Theta]2},{x,y}] 
          Solves the forward kinematics of 2-R planar manipulator whose base joint is at {x0,y0} given, link lengths 'l1','l2' and joint angles '\[Theta]1','\[Theta]2' along with the varibles for X & Y coordinates of the end-effector.       

2) inverseKinematics2R[{x0,y0},l1,l2,{x,y},{\[Theta]1,\[Theta]2}] 
          Solves the inverse kinematics of 2-R planar manipulator whose base joint is at {x0,y0} given, link lengths 'l1','l2' and end-effector coordinates 'x','y' along with the variables for the joint angles.

3) workspaceBoundary2R[{x0,y0},l1,l2,{x,y}] 
          Returns the work-space boundary equations of 2-R planar manipulator given the base joint coordinates(x0,y0) & the link lengths l1 and l2.
             
4) singularity2R[l1,l2,{\[Theta]1,\[Theta]2}] 
          Returns the singularity condition for a 2-R planar manipulator given the variables representing link lengths-'l1','l2' & joint angles-'\[Theta]1','\[Theta]2'. 

P.S : Type '?<function_name>' for more details.   
" 


forwardKinematics2R::usage = "Returns a list of rules for X & Y coordinates of the end-effector in the specified variables."
inverseKinematics2R::usage = "Returns a list of rules for the possible values of joint angles in the specified variables."
workspaceBoundary2R::usage = "Shows the equations of work-space boundaries in cartesian coordinates 'x' & 'y'."
singularity2R::usage = "Shows the condition for singularity in the specified configuration variables."


Begin["`Private`"]; 


forwardKinematics2R[{x0_,y0_},l1_, l2_, {\[Theta]1_, \[Theta]2_}, {x_, y_}] := 
  Module[{sol},
	sol = {x -> x0 + l1*Cos[\[Theta]1] + l2*Cos[\[Theta]1 + \[Theta]2], 
			y -> y0 + l1*Sin[\[Theta]1] + l2*Sin[\[Theta]1 + \[Theta]2]};
	Return[sol];
];

inverseKinematics2R[{x0_,y0_},l1_, l2_, {x_,y_}, {\[Theta]1_, \[Theta]2_}] := 
  Module[{d, p, eqFK, eqn, \[Theta]12, elim\[Theta]1, sol\[Theta]12, sol\[Theta]},
	p={x,y};
	d = Norm[p-{x0,y0}];
   (* Check if the end effector point is within the workspace *)
   
	If[NumericQ[l1] && NumericQ[l2] && Map[negativeQ, (l1 + l2) - d ||  (d - Abs[l1 - l2])],
        Print["inverseKinematics2R::error: The points x= ", x, " and y= ", y, " are beyond the workspace of the 2-R serial chain!"];
    Return[{{}}];
    ];
   
   (* Check if l1 and l2 are non-zero positive numbers *)

	If[NumericQ[l1] && NumericQ[l2] && ((Map[negativeQ, l1 || l2]) || Map[zeroQ, l1 || l2]),
		Print["inverseKinematics2R::error: The link lengths should be positive, whereas 
				l1 = ", l1, " l2 = ", l2];
		Return[{{}}];
    ];
   
	eqFK = forwardKinematics2R[{x0,y0},l1, l2, {\[Theta]1,\[Theta]2}, {x,y}] /. {\[Theta]1 + \[Theta]2 -> \[Theta]12};
	eqn = ((p /. eqFK) - p);
   
	elim\[Theta]1 = solveLinTrig2[eqn, \[Theta]1];
	sol\[Theta]12 = solveLinTrig1[elim\[Theta]1[[2]], \[Theta]12];
	sol\[Theta] = {elim\[Theta]1[[1,1]]/.#,\[Theta]2->(\[Theta]12-(\[Theta]1/.elim\[Theta]1[[1]]))/.#}&/@sol\[Theta]12;

	Return[sol\[Theta]];
];

workspaceBoundary2R[{x0_,y0_},l1_,l2_,{x_,y_}]:=
   Module[{eq1,eq2},
    eq1=(x-x0)^2+(y-y0)^2-(l1+l2)^2;
	eq2=(x-x0)^2+(y-y0)^2-(l1-l2)^2;
	Return[{eq1*eq2}];
];

singularity2R[l1_,l2_,{\[Theta]1_,\[Theta]2_}]:=
  Module[{},
	Return[l1*l2*Sin[\[Theta]2]];
];


End[]
EndPackage[]



