doAxialMeasure = true;


//get some informations about the image
image_Name = getTitle();
getDimensions(image_width, image_height, image_channels, image_slices, image_frames);
getVoxelSize(voxel_width, voxel_height, voxel_depth, voxel_unit);


doFWHM();

if (doAxialMeasure){
	// define the length of the line
	lineLength = 40;
	// to find bright spot , https://imagej.nih.gov/ij/docs/guide/146-29.html
	prominence = 100 ;
	selectImage(image_Name);
	run("Reslice [/]...", "output=["+voxel_depth+"]");
	run("Find Maxima...", "prominence="+prominence+" output=[Point Selection]");
	// get coordinates of the Maxima(s)
	getSelectionCoordinates(xpoints, ypoints);
	// make a Line using the Coordinates of the First Maxima (the brightest)
	// and the line length defined
	makeLine(xpoints[0], ypoints[0] - lineLength/2,xpoints[0], ypoints[0] + lineLength/2 );
	// do the measure
	doFWHM();
}


// companion function to make the math for FWHM
function doFWHM(){
	image_Name = getTitle();
	getVoxelSize(voxel_width, voxel_height, voxel_depth, voxel_unit);

	y = getProfile();
	x = Array.getSequence(lengthOf(y));
	Fit.doFit("Gaussian", x, y) ;
	Fit.plot;
	
	// parameter d of gaussian		
	sortedParameter = Fit.p(3); 
	rSquared = Fit.rSquared ; 
	
	// http://fr.wikipedia.org/wiki/Largeur_%C3%A0_mi-hauteur		
	FWHM = (2 * sqrt( 2 * log(2) ) ) * sortedParameter ;
	setResult("FWHM ("+voxel_unit+")", nResults, FWHM * voxel_height);
	setResult("Label",nResults-1,image_Name);
	setResult("rSquared",nResults-1,rSquared);
	updateResults();
	selectWindow("Results");
}