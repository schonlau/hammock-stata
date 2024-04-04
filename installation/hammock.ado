program define hammock
*! 1.4.2   Apr 03, 2024: removed "set matsize" (obsolete in modern Stata) and corresponding error msg
	version 14.2
	syntax varlist [if] [in], [  Missing BARwidth(real 1) MINBARfreq(int 1) /// 
		hivar(str) HIVALues(string) SPAce(real 0.0) ///
		LABel labelopt(str) label_min_dist(real 3.0) ///
		SAMEscale(varlist) ///
		ASPECTratio(real 0.72727) COLorlist(str) shape(str) no_outline  * ]

	confirm numeric variable `varlist'

	local missing = "`missing'" != ""
	local addlabel= "`label'" != ""
	local same = "`samescale'" !=""
	local outline = 1- ("`_outline'" == "no_outline") //options that start with "no" have different syntax 

	if ("`shape'"=="") local shape="rectangle"
	local rangeexpansion=0.1  /*fraction space reserved for missing values*/
	
	local varnamewidth=`space' /*=percentage of space given to text as opposed to the graph*/
	if `addlabel'!=0 & `space'==float(0) {
		local varnamewidth=0.3   // if addlabel, change the default space to 0.3
	}

	check_sufficient_colors ,  hivar("`hivar'") hivalues("`hivalues'") colorlist(`"`colorlist'"') 
	if ("`colorlist'"=="")	 local colorlist="black red  blue teal  yellow sand maroon orange olive magenta"


	* observations to use 
	marksample touse , novarlist
	qui count if `touse' 
	if r(N) == 0 { 
		error 2000 
	} 

	preserve 
	qui keep if `touse' 
	tempvar id std_ylag width graphxlag colorgroup
	tempname label_coord 
	qui gen long `id'  = _n 
	local N = _N
	local separator "@"  /*one letter only, separates variable names */
	
	// Aspectratio
	// Stata define aspectratio as height/width of the inner plotregion (see "help aspect"). I chose as default: aspect=4/5.5=0.7272
	// Because it refers to the *inner* plot region, it is not affected by the xlabels at the bottom.
	// aspectratio is explicitly specified in the syntax, but can be overwritten in "*"	
	// Outside of Stata, the aspectratio is often defined the other way around (including in the initial drafts of this program)
	// Therefore defining a second quantity:
	local ar_x= 1/`aspectratio'
	// `xsize' and `ysize' define the whole-graph aspectratio.  `aspect' defines the plot aspect ratio.
	// Stata's default is ysize(4) xsize(5.5)
	// see "help region_options"
	
	
	// if "missing" not specified, drop observations with missing values
	// do this before range transformation; or else obs may generate labels for non-missing x-vars, but obs are later removed
	if (`missing'==0) {
		foreach v of var `varlist' { 
			qui drop if `v'==.
		}
	}

	if `addlabel'==1 {
		list_labels `varlist', separator(`separator')
		// matrix has 3 cols: 1 variable values/levels, 2 variable index in `varlist',3  ./0 to mark start of new var
		matrix `label_coord'= r(label_coord)  
		local label_text  "`r(label_text)'"
	}

	local k : word count `varlist' 
	local max=`k'
	tokenize `varlist' 
	

	/*generate colorgroup variable for highlighting*/
	/*later, missval replaces "." , so gen_colorgroup needs to go before*/
	if "`hivar'"!=""    qui gen_colorgroup , hivar(`hivar') hivalues(`hivalues')  
	else 				qui gen colorgroup=1
	qui tab colorgroup
	local ncolors=r(r) // number of colors
	
	* transform variables' range to be between 0 and 100, also adjust matrix w label coordinates
	// Note this computes the midpoints (between 0 and 100), the upper/lower points may be a little wider
	// Missing values are initially assigned -10, but after standardization (std_y) they equal 0
	if `same' local addtmp = `"samescale(`samescale')"'  // if `same' specify this option
	transform_variable_range `varlist', same(`same') addlabel(`addlabel') ///
		rangeexpansion(`rangeexpansion') mat_label_coord("`label_coord'") miss(`missing') `addtmp'
	// local yline=r(yline) //not currently using yline; 
							//such a line can be passed directly as optional argument
	if (`addlabel')  decide_label_too_close ,  mat_label_coord("`label_coord'") min_distance(`label_min_dist')
	
	* construct xlabel
	local i=1
	foreach v of var `varlist' {
		if mod(`i', 2)  local xl "`xl' `i'" 
		else 			local tl "`tl' `i'" 
		local xlabel "`xlabel' `i' ``i''"
		local i = `i' + 1
	}
	local xl = substr("`xl'", 2, .)
	local xlab_num "`xl'"
	local xlab_num "`xl'`tl'" 

	
	* transform the data from one obs per data point to one obs per graphbox 
	* using reshape and contract

	* variables yvar need be listed std_y1 std_y2 std_y3 etc.   
	// "std_y" is a stub only; cannot replace by a tempvariable
	qui reshape long std_y, i(`id') 
	keep std_y `id' _j colorgroup 

	* graphx is the variable index of std_y
	local graphx "_j"

	qui { 
		bysort `id' (`graphx') : gen `std_ylag' = std_y[_n+1] 
		drop if `std_ylag' == .  /*last variable doesn't connect to variable after it*/
	} 	

	* graphx is important for unique identification 
	* contract creates _freq variable
	contract std_y `std_ylag' `graphx' colorgroup
	* scatter std_y  `graphx'

	qui replace _freq=max(`minbarfreq',_freq)
	
	
	*** preparation for graph 
	
	* make room for labels in between rectangles
	qui gen `graphxlag'= `graphx'+ (1-`varnamewidth'/2)
	if (`varnamewidth'>0)	{ 
		qui replace `graphx'= `graphx' + `varnamewidth'/2 
	}

	* compute `width' :  refers to a percentage of `range'
	summarize `std_ylag', meanonly 
	local range = r(max) - r(min)
	// Barwidth is a multiplicative increase from 0.2
	qui gen `width' =_freq / `N' * `range' * 0.2* `barwidth' 
	
	*xrange
	local xrange= `k'-1   /* number of x variables-1==xrange */
	

	//reshape: previously each obs represented unique (box,color) combination
	//now each obs represents a unique box (with multiple colors). 
	// color variables contain the width of the color and are missing otherwise
	keep  `width' colorgroup `graphx' `graphxlag' std_y `std_ylag' 
	qui reshape wide `width', j(colorgroup) i(`graphx' `graphxlag' std_y `std_ylag')
	
	// in trivial cases there may be fewer than 4 boxes. 
	// 4 obs are needed for the rectangle vars
	if (_N<4)  set obs 4 


	// labels
	label var `graphx' " "
	label var std_y " "
	lab define myxlab `xlabel'
	lab values  `graphx' myxlab
	
	if `addlabel'==1 {
		compute_addlabeltext,  mat_label_coord(`label_coord') missing("`missing'") ///
			label_text("`label_text'") separator(`separator') labelopt(`"`labelopt'"')
		local addlabeltext=r(addlabeltext)  
	}
	else  local addlabeltext=""


	// compute ylabmin (midpoint-1/2 width) and ylabmax (midpoint+1/2 width)
	// If y_std contains missing values they are coded 0, just like any other minimum value
	/* needed to avoid that some coordinates are off the graph screen*/
	/* since def of width changes later this is only approximate */
	// accumulate the width for different colors
	tempvar width_total
	qui egen `width_total' = rsum(`width'*)
	computeylablimit std_y `std_ylag' `width_total'
	local ylabmax=r(ylabmax)  
	local ylabmin=r(ylabmin)

	cap drop `width_total'  // if parallelogram, meaning of width_total changes
	if ("`shape'"=="parallelogram") {
		// iterate through computing vertical width and  updating ylabmax, ylabmin
		iterate_width_range , xstart(`graphx') xend(`graphxlag')  ystart("std_y") yend(`std_ylag') ///
			width(`width') xrange(`xrange') ar_x(`ar_x') ///
			ncolors(`ncolors') ylabmax(`ylabmax') ylabmin(`ylabmin')
		local ylabmax=r(ylabmax)
		local ylabmin=r(ylabmin)
	}
	local yrange=`ylabmax'-`ylabmin'

	// each obs represents a unique box (with multiple colors)
	// there is one `width' variable for each color: 
	// `width' is the stub of the variable names; `width'i contains the width of color i
	GraphBoxColor , xstart(`graphx') xend(`graphxlag') ystart(std_y) yend(`std_ylag') ///
		width(`width')  ylabmax(`ylabmax') ylabmin(`ylabmin') ///
		aspectratio(`aspectratio') ar_x(`ar_x') xrange(`xrange') yrange(`yrange') ///
		 xlab_num("`xlab_num'")  graphx("`graphx'") colorlist(`"`colorlist'"') ///
		 shape("`shape'") outline(`outline') ///
		 options(`"`options'"') addlabeltext(`"`addlabeltext'"') yline(`"`yline'"')

	// matrix `label_coord' is still around but not needed in GraphBoxColor 
	
end
/**********************************************************************************/
program iterate_width_range , rclass
// for parallelogram only: 
//		iterate through computing vertical width and  updating ylabmax, ylabmin
// on exit: recompute ylabmin and ylabmax   
//			(never shrinking range, could be slightly bigger than actual minima and maxima)
// on exit: width has changed (and now refers to vertical width)
	version 17
	syntax , xstart(str) xend(str) ystart(str) yend(str)  ///
			width(str) xrange(real) ar_x(real) ncolors(int) ///
			ylabmax(real) ylabmin(real)
			
	local tol_limit=0.01  //stop iterating if increase in ylabmin/ylabmax is less than tol_limit
	local ylabmax=`ylabmax'
	local ylabmin=`ylabmin'
	
	local width_ylist= "" 
	forvalues c=1/`ncolors'  {
		tempvar width_y`c'      // need one var for each color
		qui gen `width_y`c''=.
		local width_ylist= "`width_ylist' `width_y`c''"
		// macro `width_ylist' was created before loop and persists after loop;
		// the var name `width_y`c'' was created inside of the loop and does not (but var is there)
	}

	local tol=`tol_limit'+1 // initialize current tolerance
	while (`tol'>`tol_limit') { 
		// compute width 
		local yrange=`ylabmax'-`ylabmin'  // This is a key change; I used to use predefined range=100
		forvalues c=1/`ncolors' { 
			// compute width for one color  
			// strategy: call program with the same width`c' each time, but with a different yrange
			// on exit:  have different width_y`c'
			local width_y_c = word("`width_ylist'",`c')
			compute_vertical_width `xstart' `xend' `ystart' `yend'  ///
				`width'`c' `width_y_c' `xrange' `yrange' `ar_x'
		}

		// compute ylabmax, ylabmin
		cap drop `width_total'
		tempvar width_total
		qui egen `width_total' = rowtotal(`width_ylist')   //total width of a multi-color box
						//stata syntax change: rsum changed to rowtotal, but rsum still works also; 
						//both treat missing values as zero
		computeylablimit `ystart' `yend' `width_total'
		// tolerance is the increase either at the top or at the bottom range
		// only interested in increases:   labmax new -old , labmin old -new
		local tol=max( `r(ylabmax)'-`ylabmax' , `ylabmin'-`r(ylabmin)' )
		if ((`tol'>`tol_limit')) {  	// if I update ylabmax/ylabmin, I have to recalculate width 
										// no need to update if the range is shrinking
			local ylabmax=r(ylabmax)	// update only after computing tolerance
			local ylabmin=r(ylabmin)	// update only after computing tolerance
		}
		//di "Iteration: Tolerance= `=round(`tol',0.001)' ylabmax=`ylabmax'   ylabmin=`ylabmin'"
	}
	return local ylabmax `ylabmax'
	return local ylabmin `ylabmin'

	// caution: original `width' is lost after this
	forvalues c=1/`ncolors' { 
		local width_y_c = word("`width_ylist'",`c')
		qui replace `width'`c'=`width_y_c'
	}
end 
/**********************************************************************************/
// compute_addlabeltext (needed for GraphBoxColor)
// input: mat_label_coord arg is a matrix name. Matrix is manipulated. 
//			Manipulations persist after this program closes.
//			In this program the matrix is NOT changed
// label_text: Shakespeare example: label_text="adult@adolescent@child@0@1@2@3@4@5@20@0@1@2@3@4@5@6@7@female@male@"
program  compute_addlabeltext, rclass
	version 16
	syntax , mat_label_coord(str) missing(str) label_text(str) separator(str) [ labelopt(str) ]

		* the labels are overwriting the plot. Plot before the graphboxes 
		local n_labels = rowsof(`mat_label_coord')
		tokenize `label_text', parse(`separator')
		local addlabeltext="text("   // no space between "text" and "("
		forval j=1/`n_labels' { 
			if (`mat_label_coord'[`j',2] !=0 &  `mat_label_coord'[`j',3]!=0) { 
				/* 2nd col==0 crucial if matrix has empty rows, otherwise graph disappears*/
				local pos=`j'*2-1  /* positions 1,3,5,7, ...  */
				// text to plot ="``pos''"      y = `mat_label_coord'[`j',1]        x= `mat_label_coord'[`j',2]  
				// the matrix needs to be evaluated as  `=matrixelem'  ; otherwise just the name of the matrix elem appears
				local addlabeltext= `"`addlabeltext'`=`mat_label_coord'[`j',1]' `=`mat_label_coord'[`j',2]' "``pos''" "'
			}
		}
		if (`missing'==1) local addlabeltext= `"`addlabeltext' 0 1 "missing" "'
		* add label options 
		local addlabeltext= `"`addlabeltext',`labelopt')"'	
		return local addlabeltext=`"`addlabeltext'"' // does not allow empty "" returns
end 
/**********************************************************************************/
* creates matrix with coordinates for labels 
* on exit: label_text : single string separated by `separator; of all labels
* on exit: label_coord: matrix with one row for each label and two columns containing x and y coordinates, 
*      and "." for the first label of a new variable
*      where the xcoordinates range from 1...(#variables( and the y coordinates from 1..(# labels for corresponding variable)
* on exit: label_text: a single string with "y1 x1 y2 x2 ... y_nlabel x_nlabel". For later use with tokenize
program define list_labels, rclass
	version 7
	syntax varlist , separator(string)

	tempname one_ylabel label_coord

	n_level_program `varlist'
	local n_level= r(n_level)
		
	matrix `label_coord'=J(`n_level',3,0)
	local label_text ""
	local i=0	/* the ith variable (for x-axis) */
	local offset=0  /* sum of the number of levels in previous x-variables */
	
	foreach v of var `varlist' { 
		local i= `i'+1
		qui tab `v', matrow(`one_ylabel')
	local n_one_ylabel=r(r)

	local g : value label `v'
	forval  j = 1/`n_one_ylabel' {
		local w=`one_ylabel'[`j',1]
		matrix `label_coord'[`offset'+`j',2]=`i'
		matrix `label_coord'[`offset'+`j',1]=`w'
		if (`j'==1) {  
			* This marks when a new variable starts. Needed in decide_labels_too_close
			matrix `label_coord'[`offset'+`j',3]=.
		}
		if "`g'"!="" {
			local l : label `g' `w'
			if ("`l'"=="") {
				di as error "A label value for `g' has only white space, creating problems"
				di as res ""  // any subsequent statements are in result mode (in black)
			}
			local label_text "`label_text'`l'`separator'"
		}
		else {	
			/* format numbers to display */
			local format_w=string(`w',"%6.0g") 
			local label_text "`label_text'`format_w'`separator'"
		}
	   }
	   local offset=`offset'+`n_one_ylabel'
	}
	return matrix label_coord `label_coord'
	return local label_text `"`label_text'"' 

end 

/**********************************************************************************/
program define n_level_program, rclass
	version 7
	* compute the sum of the number of levels of all variables
	* each level later corresponds to a label
	syntax varlist 

	* calc the sum of number of levels of all vars
	local n_level=0
	foreach v of var `varlist' { 
		qui tab `v'
		  local temp= r(r)
		local n_level=`n_level' + `temp'
	}
	/* I used to set matsize here,
	 but in modern flavours of stata "set matsize" is obsolete
	 Stata flavour determines matsize (see help limits) */
	return local n_level `n_level'

end 
/**********************************************************************************/
program define globalminmax, rclass
	version 7
	* compute min and max values over all variables
	syntax varlist

	local globalmax=-9999999
	local globalmin=9999999
	foreach v of var `varlist' { 
		qui sum `v'
		if (r(max)>`globalmax') {
			local globalmax=r(max)
		}
		if (r(min)<`globalmin') {
			local globalmin=r(min)
		}
	}
	return local globalmin `globalmin'
	return local globalmax `globalmax'

end 
/**********************************************************************************/
program define compute_vertical_width
   version 7
   * compute the difference of the y coordinates needed 
   * width is right-angle width (distance between two parallel lines)
   * ar_x: (aspectratio) how much longer is the x axis relative to the y axis
   * on input : width_y has already been "generate"d 
   * width_y is computed as a result  
   *	(if width has missing values, then width_y also has missing values)
   
   args xstart xend ystart yend width width_y xrange yrange ar_x

   *xdiff is the fraction of space relative to the range occupied by any particular parallelogram. 
   *  For example, with 5 variables the xrange=5-1=4.   
   *  Going from the 2nd to the third x-variable has an xdiff = 1/ 4  * ar_x .
   *ydiff is the same thing for y, i.e. the fraction that y drops/increases for one box relative to its range.
   tempvar xdiff ydiff
   qui gen `xdiff'= .
   qui gen `ydiff'= .
   qui replace `xdiff'= (`xend'-`xstart')/`xrange' * `ar_x'
   qui replace `ydiff'= (`yend'-`ystart')/ `yrange'
   
   * The centerline of the parallelogram, xdiff and ydiff from a triangle with a right angle. 
   * theta is the angle between xdiff and the centerline of the parallelogram.
   * From this big triangle:  
   *     Cos(theta) = xdiff/  sqrt(xdiff^2 +ydiff^2) 
   * The angle theta reappears for width and width_y, as part of the parallelogram: 
   *     Cos(theta)= width/ width_y  .
   * Putting these together: 
   *    Width/ width_y = xdiff/ sqrt(xdiff^2 +ydiff^2)
   * Solving for width_y : 
   *            width_y = width  / xdiff   * sqrt(xdiff^2 +ydiff^2)
   qui replace `width_y'=`width' / `xdiff' * sqrt(`xdiff'*`xdiff'+`ydiff'*`ydiff') 

end 
/**********************************************************************************/
* distinguish between colorgroup for strings and color group numeric
* string variables allowed only when highlighting a variable not in varlist
program define gen_colorgroup
	version 16
	syntax ,  hivar(varname) HIVALues(string)

		capture confirm numeric variable `hivar'
		if !_rc {
			qui gen_colorgroup_num , hivar(`hivar') hivalues(`hivalues')  //numeric variable
		}
		else {
			gen_colorgroup_str , hivar(`hivar') hivalues(`hivalues')  //string variable
		} 
end
/**********************************************************************************/
* generate the colorgroup variable,numeric
* all values not mentioned in hivalues get color=1
* colors are generated in sequence of hivalues
* if value in hivalues not contained in hivar, then the corresponding color is just skipped
* if hivar does not exist will give error message "variable `hivar' not found
* if the same value is contained multiple times, only the last one is used and colors 
     *corresponding to earlier ones are skipped
* if hivalues contains more than 8 values earlier pens are being reused
* if hivalues contains all values of hivar, then the default color is never used
* hivalues does not accept labels at this point
* if hivalues not specified will give appropriate error
program define gen_colorgroup_num
	version 7
	syntax  ,  hivar(varname) HIVALues(numlist missingokay) 

	local pen=1
	qui gen colorgroup=`pen'
	foreach v of numlist `hivalues' {
		local pen=mod(`pen',8)+1
		qui replace colorgroup=`pen' if `hivar'==`v' 
		if ("`v'"==".") {
			qui replace colorgroup=`pen' if  `hivar'==.
		}
	}
	* list  `hivar' colorgroup

end
/**********************************************************************************/
* generate the colorgroup variable, string
* same as for numeric; except the string comparison instead of ==  , and no missing values, and loop over i
program define gen_colorgroup_str
	version 16
	syntax  ,  hivar(varname) HIVALues(string) 

	local pen=1
	qui gen colorgroup=`pen'
	foreach v in `hivalues' {
		local i=1
		local pen=mod(`pen',8)+1
		while (`i'<=_N) {
			qui replace colorgroup=`pen' in `i' if  ustrregexm(`hivar'[`i'],"`v'")    // `hivar'==`v' 
			local i=`i' + 1
		}
	}
    *list  `hivar' colorgroup

end
/**********************************************************************************/
program define computeylablimit, rclass
	version 7
	* compute maximum and minimum values for y to define graph region ylab
	* ylabmax and ylabmin are rounded to full numbers
	* If I attempt to draw beyond the graph region the graph may draw triangles instead of parallelograms
	args ystart yend width 

	local ylabmin =0
	local ylabmax = 100

	tempvar endmax startmax
	qui gen `endmax' = `yend'+`width'/2
	qui gen `startmax'= `ystart' +`width'/2
	qui sum `endmax'
	local ylabmax= r(max)
	qui sum `startmax'
	local max2=r(max)
	if (`max2'>`ylabmax') {
		local ylabmax=`max2'
	} 
	local ylabmax=round(`ylabmax',1)
	
	// compute ylabmin (copied from Tiancheng's code)
	tempvar endmin startmin
	qui gen `endmin' = `yend'-`width'/2
	qui gen `startmin'= `ystart' -`width'/2
	qui sum `endmin'
	local ylabmin= r(min)
	qui sum `startmin'
	local min2=r(min)
	if (`min2'<`ylabmin') {
		local ylabmin=`min2'
	}
	local ylabmin=round(`ylabmin',1)
	
	
	return local ylabmax `ylabmax'
	return local ylabmin `ylabmin'

end

/**********************************************************************************/
* draw all boxes
* each obs has a unique combination of "xstart xend ystart yend". 
* If more than one width`i' have non missing values then the box has more than one colors 
* xstart xend ystart yend 	(vectors determining 2 midpoints of parallelogram)
* 				start and end of the line segments for the center of the parallelogram
* width 		width`i' contain the  width of the parallelogram for color i in 1/9
*				The width of a box is the sum of all the color segments width`i'
*				For parallelogram, this is the vertical width.
*				For rectangles, this is the rightangle width.
* ylabmin,ylabmax:  upper and lower ylabels for plotting. These can be slightly wider 
*					than range of y because of the size of the bars
* aspectratio:   Aspect ratio y/x (stata's definition). 
* ar_x			 Aspect ratio x/y.   (Safer to carry forward than to recompute)
* xrange : 		 #xvars -1. It's safer to carry it forward rather than to recompute it.
* xlabnum		 :   needed to construct xlabel
* graphx 		 : name of x variable for plotting
* addlabeltext   : string for adding labels:  text("ypos1 xpos1 "text1"  ypos2 xpos2 "text2" [...])
* colorlist		 : list of colors to be used
* outline		 : 1/0: whether or not boxes should have an outline (using lcolor)
*
* Strategy: 
* For each box, 3 variables (yhigh, ylow, and x) are created. 
* For a parallelogram, they have 2 observations representing beginning and end of the parallelogram
* For a rectangle, need 4 observations for the 4 different x-values for each corner
* The plotting command for each box is added to `addplot'. 
* Thoughts: 
*	-The number of variables could be  greatly reduced by using one set of 3 variables 
*		with many obs rather than many sets of variables with 4 obs. 
*		However, then may have to reset number of observations, 
*		and unclear whether reducing the number of variables improves anything.

program define GraphBoxColor 
	version 14.2
	syntax , xstart(str) xend(str) ystart(str) yend(str) width(str) graphx(str) ///
		ylabmin(real) ylabmax(real) xlab_num(str) shape(str) outline(int) ///
		aspectratio(real) ar_x(real) xrange(real) yrange(real) ///
		[  colorlist(str) addlabeltext(str) yline(str)  options(str)  ]

	tempvar w increment
	qui egen `w'= rsum(`width'*)  // total width of a multi-color box
	qui gen `increment'=0	// sum of w_k before the current color
		
	local no_outline=""
	if (`outline'==0)	local no_outline = "lcolor(black%0)"  // boxes have no, i.e. translucent outline
	//If there is an outline, my tests indicate that the outline is correctly inside the area as in "lalign(inside)" 
	//      lalign(outside) would be a problem as there would be overlap with neighboring lines
	

	// initialize, so I can use "replace" later
	tempvar yhigh ylaghigh ylaglow xlow xhigh xlaghigh xlaglow ylow
	foreach var of newlist `xlow' `xhigh' `xlaghigh' `xlaglow' `yhigh' `ylaghigh' `ylaglow' `ylow' {
		qui gen `var'=.
	}
	// when "`shape'"=="rectangle", initialize 
	if ("`shape'"=="rectangle") {
		// o_tempvars hold the outer corners of the entire rectangle
		tempvar o_ylow o_yhigh o_ylaglow o_ylaghigh o_xlow o_xhigh o_xlaglow o_xlaghigh 
		foreach var of newlist `o_ylow' `o_yhigh' `o_ylaglow' `o_ylaghigh' `o_xlow' `o_xhigh' ///
			`o_xlaglow' `o_xlaghigh'  {
				qui gen `var'=.
		}
		init_rectangle, xstart(`xstart') xend(`xend') ystart(`ystart') yend(`yend') w(`w') ///
			o_ylow(`o_ylow') o_yhigh(`o_yhigh') o_ylaglow(`o_ylaglow') o_ylaghigh(`o_ylaghigh') ///
			o_xlow(`o_xlow') o_xhigh(`o_xhigh') o_xlaglow(`o_xlaglow') o_xlaghigh(`o_xlaghigh') ///
			xrange(`xrange') ar_x(`ar_x') ylabmin(`ylabmin') ylabmax(`ylabmax')
	}

	local N=_N
	local addplot1 ""

	// for each color
	foreach k of numlist 1/9 { 
		// draw all parallelograms of one color 
		local w_k="`width'`k'"
		capture confirm variable `w_k'
		if !_rc{ 
			/* variable w_k exists*/
			if "`shape'"=="parallelogram" {
				compute_4corners_parallelogram,  increment(`increment') w_k(`w_k') w(`w') ///
					ystart(`ystart') yend(`yend')  xstart(`xstart') xend(`xend') ///
					ylow(`ylow') yhigh(`yhigh') ylaglow(`ylaglow') ylaghigh(`ylaghigh') ///
					xlow(`xlow') xhigh(`xhigh') xlaglow(`xlaglow') xlaghigh(`xlaghigh')  
			}
			else if "`shape'"=="rectangle" {
				compute_4corners_rectangle,  increment(`increment') w_k(`w_k') w(`w') ///
					o_ylow(`o_ylow') o_yhigh(`o_yhigh') o_ylaglow(`o_ylaglow') o_ylaghigh(`o_ylaghigh') ///
					o_xlow(`o_xlow') o_xhigh(`o_xhigh') o_xlaglow(`o_xlaglow') o_xlaghigh(`o_xlaghigh') ///
					ylow(`ylow') yhigh(`yhigh') ylaglow(`ylaglow') ylaghigh(`ylaghigh') ///
					xlow(`xlow') xhigh(`xhigh') xlaglow(`xlaglow') xlaghigh(`xlaghigh')
			}
			else {
				di as error "Unrecognized shape (2)"
				exit
			}

			
			// translate k into a color 
			tokenize `"`colorlist'"'
			local color= "``k''"  // need double single quotes

			// for a given color k, plot one parallelogram at a time.			
			foreach i of  numlist 1/`N' {
				
				// for each  graph box
				local yhigh1= `yhigh'[`i']
				local ylaghigh1= `ylaghigh'[`i']
				local ylow1=`ylow'[`i']
				local ylaglow1=`ylaglow'[`i']
				local xhigh1= `xhigh'[`i']
				local xlaghigh1= `xlaghigh'[`i']
				local xlow1=`xlow'[`i']
				local xlaglow1=`xlaglow'[`i']
				
				
				
				// if graph box is non-missing (checking 1 macro; they should be all missing or all not missing)
				if (`yhigh'[`i']!=. ) {
					// temp variables die at end of a program; therefore define here where the plotting occurs
					tempvar hi`i'`k' lo`i'`k' x`i'`k'
					quietly {
						gen `hi`i'`k''=.
						gen `lo`i'`k''=.
						gen `x`i'`k''=.
					}
					if ("`shape'"=="parallelogram") { 
						// still using this because it's a little faster than quadrangle
						// for entire hammock.ado speed advantage using this routine is 1%-6%
						plot_parallelogram `hi`i'`k'' `lo`i'`k'' `x`i'`k'', ///
							yhigh(`yhigh1') ylaghigh(`ylaghigh1') x1(`xlow1') x2(`xlaglow1') ///
							ylow(`ylow1')  ylaglow(`ylaglow1') color("`color'") `no_outline'
					}
					else {						
						// works for both general quadrangles and for parallelogram as a special case
						plot_quadrangle `hi`i'`k'' `lo`i'`k'' `x`i'`k'', color("`color'") `no_outline' ///
							xhigh(`xhigh1') yhigh(`yhigh1') ///	     
							xlow(`xlow1') ylow(`ylow1') ///
							xlaghigh(`xlaghigh1') ylaghigh(`ylaghigh1') ///
							xlaglow(`xlaglow1') ylaglow(`ylaglow1') 	
					}
					local temp = r(addplot)
					assert `"`temp'"'!=""
					// caution : addplot might get tooo long
					if (`"`addplot'"'=="")  local addplot  `"`temp'"'
					else 					local addplot  `"`addplot' || `temp'"'
				}
			}						
			qui replace `increment' = `increment'+ `w_k' if `w_k'!=.   // for the next color
		}
	}
	
	// *most work happens in `addplot' (all parallelograms) and `addlabeltext' (all labels)
	// the rest just sets up the coordinates but plots nothing because msymbol(none)
	// *xscale and yscale are explicitly specified because they were needed earlier to compute aspectratio
	// *The twoway command cannot be moved up to the caller program because the caller program does not have 
	// 	access to the many (e.g.200) variables used in `addplot'

	/* 
	//for debugging: This helps seeing the frame in the plot area on which aspectratio is built.
	twoway scatter std_y `graphx', ///
		ylab(`ylabmin' `ylabmax')  xlab(`xlab_num',) ylab(,valuelabel)  `yline'     ///
		legend(off) ytitle("") xtitle("") msymbol(none)  ///
		plotregion(lstyle(dashed) m(zero) ) yscale(on)  ///
		aspect(`aspectratio') `options'  ||  `addplot' `addlabeltext'
	*/
	
	twoway scatter std_y `graphx', ///
		ylab(`ylabmin' `ylabmax')  xlab(`xlab_num',valuelabel noticks) ylab(,valuelabel noticks)  `yline'     ///
		legend(off) ytitle("") xtitle("") yscale(off) xscale(noline) msymbol(none)  ///
		plotregion(style(none) m(zero)) ///
		aspect(`aspectratio') `options'  ||  `addplot' `addlabeltext'


	// These options below produce a graph without any margins except for the x-labels at the bottom
	// scatter y x,plotr(m(zero)) graphr(m(zero)) ytitle("") xtitle("") yscale(off) xscale(noline) 
	//  graphr(m(zero)) : I can't set these to zero because the numbers go off the right/left side

end 

///////////////////////////////////////////////////////////////////////////////////////////////////////
//  prepare plotting a single parallelogram: 
//  create command `addplot' and create three variables needed for the plot
//  rarea needs a different observation for each x-value. A parallelogram has 2 x-values, 
//      so there are 2 observations here.
program define plot_parallelogram, rclass
	version 14.2
	syntax varlist (min=3 max=3), yhigh(real) ylow(real) x1(real) ylaghigh(real) ylaglow(real) x2(real) [ * ]

	if (`yhigh'==. | `x1'==. | `ylaghigh'==.) {
	  di as error "Bug: Missing values in plot_parallelogram"
	  di "`yhigh';  `ylow'; `x1';  `ylaghigh';  `ylaglow';  `x2' " 
	  exit 
	}
	tokenize `varlist'
	local hi `1'
	local lo `2'
	local x `3'

	quietly{
		replace `x'= `x1' in 1
		replace `x'= `x2' in 2
		replace `hi'= `yhigh' in 1
		replace `lo'= `ylow' in 1
		replace `hi'= `ylaghigh' in 2
		replace `lo'= `ylaglow' in 2
	}

	//twoway rarea  `hi' `lo' `x', `options'
	
	// restricting to two obs "in 1/2" can reduce speed by  ~25%
	local addplot `"rarea `hi' `lo' `x' in 1/2, `options'"'
	return local addplot `"`addplot'"'

end
///////////////////////////////////////////////////////////////////////////////////////////////////////
//  prepare plotting a single quadrangle  (can be used for rectangles)
//  create command `addplot' and create three variables needed for the plot
//  varlist has only 3 variables, because rarea only needs 3 variables
//  rarea needs a different observation for each of the four different x-values;
//     the rectangle is plotted in 3 segements between the 4 x-values
program define plot_quadrangle, rclass
	version 17.0
	syntax varlist (min=3 max=3), ///
		xhigh(real) yhigh(real) ///	     
		xlow(real) ylow(real) ///
		xlaghigh(real) ylaghigh(real) ///
		xlaglow(real) ylaglow(real)  [ * ]

	if (`yhigh'==. | `xhigh'==. | `ylaghigh'==.) {
	  di as error "Bug: Missing values in plot_quadrangle"
	  di "`yhigh';  `ylow'; `xhigh';  `ylaghigh';  `ylaglow';  `xlaghigh' " 
	  exit 
	}
	tokenize `varlist'
	local yhi `1'
	local ylo `2'
	local x `3'

	if (`yhigh' <= `ylaghigh'){

		// compute y2hi, the missing upper point for the first sub-area
		slopes, xleft(`xhigh') yleft(`yhigh')  xright(`xlaghigh') yright(`ylaghigh')
		local y2hi= `=r(intercept)' + `=r(slope)' * `xlow'

		// compute y3lo, the missing lower point for the third sub-area
		slopes, xleft(`xlow') yleft(`ylow')  xright(`xlaglow') yright(`ylaglow')
		local y3lo= `=r(intercept)' + `=r(slope)' * `xlaghigh'

		// for each x-value, give the upper and lower y-values of the area (could be identical)
		quietly{
			replace `x'= `xhigh' in 1
			replace `yhi'= `yhigh' in 1
			replace `ylo'= `yhigh' in 1

			replace `x'= `xlow' in 2
			replace `yhi'= `y2hi' in 2
			replace `ylo'= `ylow' in 2

			replace `x'= `xlaghigh' in 3
			replace `yhi'= `ylaghigh' in 3
			replace `ylo'= `y3lo' in 3

			replace `x'= `xlaglow' in 4
			replace `yhi'= `ylaglow' in 4
			replace `ylo'= `ylaglow' in 4
		}
	}
	else{
		// compute y2hi, the missing upper point for the first sub-area
		slopes, xleft(`xlow') yleft(`ylow')  xright(`xlaglow') yright(`ylaglow')
		local y2low= `=r(intercept)' + `=r(slope)' * `xhigh'

		// compute y3lo, the missing lower point for the third sub-area
		slopes, xleft(`xhigh') yleft(`yhigh')  xright(`xlaghigh') yright(`ylaghigh')
		local y3high= `=r(intercept)' + `=r(slope)' * `xlaglow'

		// for each x-value, give the upper and lower y-values of the area (could be identical)
		quietly{
			replace `x'= `xlow' in 1
			replace `yhi'= `ylow' in 1
			replace `ylo'= `ylow' in 1

			replace `x'= `xhigh' in 2
			replace `yhi'= `yhigh' in 2
			replace `ylo'= `y2low' in 2

			replace `x'= `xlaglow' in 3
			replace `yhi'= `y3high' in 3
			replace `ylo'= `ylaglow' in 3

			replace `x'= `xlaghigh' in 4
			replace `yhi'= `ylaghigh' in 4
			replace `ylo'= `ylaghigh' in 4
		}
	}

	//twoway rarea  `yhi' `ylo' `x', `options'
	
	local addplot `"rarea `yhi' `ylo' `x' in 1/4, `options'"'
	return local addplot `"`addplot'"'

end
////////////////////////////////////////////////////////
// compute the "outer" corners of a multi-color rectangle
// input: 	xstart, xend, ystart, yend, w 
//			reals: xrange, ylabmin ylabmax, ar_x  (needed to compute alpha's)
// output: 	o_ylow o_yhigh o_ylaglow o_ylaghigh 
//			o_xlow o_xhigh o_xlaglow o_xlaghigh 
	// Notation as used in the paper's appendix : 
	// B' = o_ylow , o_xlow
	// A'=  o_yhigh, o_xhigh
	// M= xstart, ystart 
	// AC=w
	// A,B,C,E  do not have existing variables and need be computed
	// xend, yend have no corresponding letter, but it's the middle between Afar and Cfar
// note convert:  degrees = radians * 180/_pi 
program init_rectangle
	version 17
	syntax, xstart(str) xend(str) ystart(str) yend(str) w(str) ///
		o_ylow(str) o_yhigh(str) o_ylaglow(str) o_ylaghigh(str) ///
		o_xlow(str) o_xhigh(str) o_xlaglow(str) o_xlaghigh(str) ///
		xrange(real) ylabmin(real) ylabmax(real) ar_x(real)
			
	// y  (~0..~100) has a different scale than x (1..#xvar-1)
	// When plotting, use the variables on their usual scale
	// For computations with angles, need to account for the different scales and aspectratio
	// the angle refers to the physical distance on the plot, not to distances with different ranges in x and y
	
	// compute alpha, cosalpha, sinalpha, for each graph box
	tempvar deltax deltay alpha cosalpha sinalpha 
	qui gen `alpha'=.
	qui gen `cosalpha'=.   // cos(alpha)
	qui gen `sinalpha'=.   // sin(alpha)
	// deltax and deltay are not the physical scatter plot distance, but const multiplier cancels in ratio
	// need to calculate the y max and min difference for standardize delta y
	local ymaxmin_diff= (`ylabmax'-`ylabmin')
	qui gen `deltax'= (`xend'-`xstart')/`xrange' *`ar_x' //same calc as in compute_vertical_width
	qui gen `deltay'= (`yend'-`ystart')/`ymaxmin_diff'   			 	//same calc as in compute_vertical_width

	// Loop because trigonometric functions don't take variables as arguments. 
	// (could convert to a matrix as an alternative.)
	// Alpha appears correct; The default layout atan(5.5/4) = 54 degrees which is similar to simple example
	foreach i of numlist 1/`=_N' { 
		local ratio= `deltay'[`i']/`deltax'[`i']  
		qui replace `alpha' = atan(`ratio')  in `i' //alpha is in radians
		qui replace `cosalpha'=cos(`alpha'[`i']) in `i'
		qui replace `sinalpha'=sin(`alpha'[`i']) in `i'
	}
	
	quietly{
		tempvar vertical_w 
		gen `vertical_w'= `w' / `cosalpha'   // alpha must be in radians

		// A
		tempvar ax ay
		gen `ax' = `xstart'
		gen `ay' = `ystart' + 0.5 * `vertical_w'

		// B
		tempvar bx by
		gen `bx' = `xstart'
		gen `by' = `ystart' - 0.5 * `vertical_w'

		// C 
		tempvar cx cy

		// sinalpha was computed on scaled w, scaled cx.
		// Unscaling is achieved by replacing:  cx_scaled = cx/ xrange *ar_x ,etc.
		// for y the unscaling cancels out	
		gen `cx' = `ax'+ `w'*`sinalpha'  /`ymaxmin_diff' * `xrange'/`ar_x'
		gen `cy' = `ay'- `w'*`cosalpha'

		// E  
		tempvar ex ey 
		gen `ex' = (`ax' + `cx')/2
		gen `ey' = (`ay' + `cy')/2

		// A'  
		replace `o_xhigh' =  `ax' + `xstart' - `ex' 
		replace `o_yhigh' =  `ay' + `ystart' - `ey'

		// C'  
		replace `o_xlow' =  `cx' + `xstart' - `ex'  
		replace `o_ylow' =  `cy' + `ystart' - `ey'

		// "Afar"  Starting from the far midpoint , go up 
		replace `o_xlaghigh' = `xend' + (`o_xhigh' - `xstart')
		replace `o_ylaghigh' = `yend' + (`o_yhigh' - `ystart') 

		// "Cfar"  Starting from the far midpoint , go down
		replace `o_xlaglow' = `xend' - (`o_xhigh' - `xstart')
		replace `o_ylaglow' = `yend' - (`o_yhigh' - `ystart')
	}
	
	// for debugging: look at one box (=observation) at a time
	// this only works for default aspect ratio
	/*
	local box=1
	twoway scatteri `=`ay'[`box']' `=`ax'[`box']'  "A"  ///
					`=`by'[`box']' `=`bx'[`box']'  "B"  ///
					`=`cy'[`box']' `=`cx'[`box']'  "C"  ///
					`=`ey'[`box']' `=`ex'[`box']'  "E"  ///
					`=`o_ylow'[`box']' `=`o_xlow'[`box']'  "o_low"  ///
					`=`o_yhigh'[`box']' `=`o_xhigh'[`box']'  "o_high"  ///
					`=`o_ylaglow'[`box']' `=`o_xlaglow'[`box']' (9) "o_laglow"  ///
					`=`o_ylaghigh'[`box']' `=`o_xlaghigh'[`box']' (9) "o_laghigh" ///
					`=`yend'[`box']' `=`xend'[`box']' (9) "Mid end" ///
					, name(box`box', replace) plotr(m(zero)) ///
					ytitle("") xtitle("") 
	// view using  "graph display box<#>"
	*/
end
////////////////////////////////////////////////////////
// compute a,b in simple linear regression
// opted not to use "regress" because that requires variables 
//	and I would use preserve/ restore which is too slow
program slopes , rclass
   version 17
   syntax , xleft(real) xright(real) yleft(real) yright(real)

   local n=2 // sample size 
   local sumx= `xleft'+`xright'
   local sumy= `yleft'+`yright'
   local sumxsqr= (`xleft')^2+(`xright')^2
   local crossxy = `xleft'*`yleft'+`xright'*`yright'

   local num= `n' * `crossxy' - (`sumx' * `sumy')
   local denum=`n'*`sumxsqr' - (`sumx')^2
   local b=`num'/`denum'
   local a=(`sumy' - `b'*`sumx' )/`n'
   
   return local slope `b'
   return local intercept `a'
end
///////////////////////////////////////////////////////////////////
// compute coordinates of 4 corners for all graph boxes of the same color
// each obs has a unique combination of "xstart xend ystart yend" (multi-color box)
// but this program only deals with the color specified by w_k
// input: 
// w_k: 	name of width variable for color k. 
// w:		total width of a graph box (across all colors)
// increment: sum of w_k for all previous (k-1) colors
// ystart, yend, xstart, xend: names of input coordinates 
// other strings are names of output coordinates of 4 corners
program compute_4corners_parallelogram 
	version 17
	syntax , increment(str) w_k(str) w(str) ///
		ystart(str) yend(str) xstart(str) xend(str) ///
		ylow(str) yhigh(str) ylaglow(str) ylaghigh(str) ///
		xlow(str) xhigh(str) xlaglow(str) xlaghigh(str) 
	
	quietly {
		replace `ylaglow' =  (`yend'   -   `w' / 2)  // outer lower corner
		replace `ylow' =  (`ystart' -  `w' / 2)    // outer lower corner
		
		replace `ylaglow' = `ylaglow' + `increment' 
		replace `ylow' = `ylow'+ `increment' 
		replace `ylaglow'=. if `w_k'==.
		replace `ylow'=. if `w_k'==.
		replace `ylaghigh' = `ylaglow' +  `w_k'
		replace `yhigh' = `ylow' +  `w_k'
		replace `xhigh'=`xstart'  // for parallelogram xhigh and xlow  are the same
		replace `xlow'=`xstart'
		replace `xlaghigh'=`xend'  // for parallelogram xlaghigh and xlaglow  are the same
		replace `xlaglow'= `xend'
	}
end 
///////////////////////////////////////////////////////////////////////////////////////////
// compute corners of an inside rectangle for color k
// strategy:
// a) (previously) compute the outer 4 corners (o_*) of the whole multi-color box 
// b) compute 2 left corners for the current box
//	  by moving `increment'/`w' along the vector (o_xlow,o_ylow)->(o_xhigh,o_yhigh) 
//	  compute 2 right corners for the current box
//	  by moving (`increment'+`w_k')/`w' along the vector (o_xlow,o_ylow)->(o_xhigh,o_yhigh)
// This program is called once per color
// input : coordinates of the outer rectangle (multi-color box)
// output: coordinates of the inner rectangle for color `k' , or missing if w_k==. 
program compute_4corners_rectangle
	version 17
	syntax , increment(str) w_k(str)  w(str) ///
		ylow(str) yhigh(str) ylaglow(str) ylaghigh(str) ///
		xlow(str) xhigh(str) xlaglow(str) xlaghigh(str) ///
		o_ylow(str) o_yhigh(str) o_ylaglow(str) o_ylaghigh(str) ///
		o_xlow(str) o_xhigh(str) o_xlaglow(str) o_xlaghigh(str) 
	
	quietly {
		// lower left 
		replace `xlow' = `o_xlow'+ `increment'/`w' * (`o_xhigh'-`o_xlow')
		replace `ylow' = `o_ylow'+ `increment'/`w' * (`o_yhigh'-`o_ylow')
		replace `ylow'=. if `w_k'==.  
		// upper left 
		replace `xhigh' = `o_xlow'+ (`increment'+`w_k')/`w' * (`o_xhigh'-`o_xlow')
		replace `yhigh' = `o_ylow'+ (`increment'+`w_k')/`w' * (`o_yhigh'-`o_ylow')
		// lower right 
		replace `xlaglow' = `o_xlaglow'+ `increment'/`w' * (`o_xlaghigh'-`o_xlaglow')
		replace `ylaglow' = `o_ylaglow'+ `increment'/`w' * (`o_ylaghigh'-`o_ylaglow')			
		replace `ylaglow'=. if `w_k'==.
		// upper right 
		replace `xlaghigh' = `o_xlaglow'+ (`increment'+`w_k')/`w' * (`o_xlaghigh'-`o_xlaglow')
		replace `ylaghigh' = `o_ylaglow'+ (`increment'+`w_k')/`w' * (`o_ylaghigh'-`o_ylaglow')	
	}
end
///////////////////////////////////////////////////////////////////////////////////////////
// colorlist should not be smaller than number of values in hivalues + 1 (default color)
program define check_sufficient_colors
	version 16
	syntax,  [ hivar(str) HIVALues(string)  COLorlist(str) ] 

	if ("`hivar'"!="" & "`colorlist'"!="") {
		local hicount: word count `hivalues'
		
		// colorlist may contain strings, so "word count" does not work
		local colcount=0
		tokenize `"`colorlist'"'
		while "`*'" ~= "" {
			macro shift
			local colcount= `colcount' +1 
		}
		if (`colcount'<`hicount'+1) {
			di as error " There are not enough colors for highlighting. "
			di as error "(The first color is the default color and is not used for highlighting.)"
			error 2000
		}
	}
end
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// purpose :  	transform variables to have a range between 0 and 100, 
//				compute corresponding label positions
//				Specifically, this refers to the center positions of the parallelograms.
//				(Later, range will be extended to accommodate upper and lower points, and missing vals)
// output: variables in varlist changed (we are in "preserve" mode, so changes disappear at the end)
// output: label_coord matrix changed if labels are used
// output: std_y`i' variables generated (lowest values (e.g. missing values) are standardized to 0)
// note: miss is an int; previously option missing was present/absent
program transform_variable_range , rclass
	version 16.0
	syntax varlist ,  addlabel(int) same(int) mat_label_coord(str) miss(int) ///
		rangeexpansion(real) [ samescale(varlist) ]

	* compute max,min values over the variables specified in samescale
	if `same' {
		globalminmax `samescale' 
		local globalmax=r(globalmax)
		local globalmin=r(globalmin)
		local globalrange=`globalmax'-`globalmin'
	}

	local i=0
	foreach v of var `varlist' { 
		local i = `i' + 1   /*needed to find the right variables of label_coord */
		* new tempvar for each loop iteration 
		quietly sum `v'
		if (r(min)!=r(max)) {
			local range=r(max)-r(min)
		}
		else { 
			* var has only one value
			local range=1	
		}
		local min=r(min)
		if (`same' & ustrregexm("`samescale'","`v'")) {
			* override calculations with variables specified in samescale
			local range=`globalrange'
			local min=`globalmin'
		}
		if (`miss') {
			* missing option specified 
			local yline "yline(4, lcolor(black))" 	// horizontal line to separate missing values
			local missval=`min'-`rangeexpansion'*`range'
			qui replace `v'=`missval' if `v'==.
			local min=`missval'  /* redefine minimum value and range */
			local range=`range'*(1+`rangeexpansion')
		}	
		local my_y  "std_y`i'"
		qui gen `my_y' = (`v'-`min')/ `range' * 100
		
		if `addlabel'==1 {
			local n_labels = rowsof(`mat_label_coord') 
			forval ii = 1 / `n_labels' {
				if (`mat_label_coord'[`ii',2]==`i') { 
					/* label belongs to variable `v' */
					matrix `mat_label_coord'[`ii',1]= (`mat_label_coord'[`ii',1] -`min')/ `range' * 100
				}
			}
		}
		
	}
	return local yline `"`yline'"'
end
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// purpose :  	decide whether the  label distance on the y-axis is too close. If yes, set mat_label_coord[,3]==0
//				This has to happen *after* transforming variables to have a range between 0 and 100
// input :	  mat_label_coord[,3] has "." for the first element of each variable (where no distances are computed)
//			  min_distance: minimum distance between labels  (on a scale from 0-100)
// assume :		labels are sorted in ascending order
// output:	  mat_label_coord[,3] contains decision (0= do not plot, 1= plot, .= bottom most label (also first label of new var in sequence)
program decide_label_too_close 
	version 16.0
	syntax ,  mat_label_coord(str) min_distance(real)

	local n_labels = rowsof(`mat_label_coord') 
	local add_dist=0  // previous distance to previously plotted label
	forval j = 2 / `n_labels' {
		if (`mat_label_coord'[`j',3]!=.) {
			/* if not the first label of a new variable,  distance = y_current - y_last */
			local distance= `add_dist' + (`mat_label_coord'[`j',1]- `mat_label_coord'[`j'-1,1])
			if (`distance'<`min_distance') {
				matrix `mat_label_coord'[`j',3]=0  // do not plot
				local add_dist= `distance' // distance to previously plotted label
			}
			else {
				matrix `mat_label_coord'[`j',3]= 1
				local add_dist=0
			} 
		 }
		else {
			local add_dist=0  // reset add_dist for new variable
		}
	}
	
end 
///////////////////////////////////////////////////////////////////////////////////////////////
//Version 
//*! 1.0.0   21 February 2003 Matthias Schonlau 
//*! 1.0.1   2 March 2003: Allow variables with one value, Bug fixed when barwidth was too large \ 
//*! 1.0.2   13 March 2003: width changed to right-angle-width
//*! 1.0.3   20 Nov 2003: hivar and hival options implemented \
//*! 1.0.4   no changes \
//*! 1.0.5   4 Nov 2008: added samescale option, fixed bug related to >8 colors
//*! 1.1.0    ongoing 2017: major rewrite due to use of twoway rarea for parallelograms
//*! 1.1.1    17 May 2022: added labeltextsize option
//*! 1.1.2    19 May 2022: added argument to samescale option
//*! 1.1.3    30 June 2022: fixed bug when not enough colors specified for highlighting
//*! 1.1.4    July 7 2022: missing values can now be highlighted
//*! 1.2.0    August 11, 2022: shape rectangle added 
//*! 1.2.1    Nov 15, 2022: rectangle is default shape, space allowed w/o label, update helpfile
//*! 1.2.2    Jan 27, 2023: minbarfreq option added
//*! 1.2.3    Apr 12, 2023: Added warning if label value contains only white space
//*! 1.2.4   May 24, 2023: label_min_dist option
//*! 1.2.5   Oct 10, 2023: added option outline 
//*! 1.2.6   Nov 10, 2023: allow hivar to be a string variable 
//*! 1.2.7   Nov 13, 2023: remove display statement leftover from debugging
//*! 1.2.8   Nov 16, 2023: fixed bug related to "labelopt" 
//*! 1.2.9   Nov 17, 2023: colorlist now allows RGB values
//*! 1.3.0   Jan 26, 2024: removed duplicate code, improved documentation compute_vertical_width
//*! 1.4.0   Feb 21, 2024: For parallelogram, iterate between computing range and vertical width
//*! 1.4.1   Apr 01, 2024: fixed bug plotting labels for obs with missing values when not specifying missing
//*! 1.4.2   Apr 03, 2024: removed "set matsize" (obsolete in modern Stata) and corresponding error msg