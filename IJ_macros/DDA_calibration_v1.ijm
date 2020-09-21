run("Clear Results");
run("Set Measurements...", "integrated redirect=None decimal=3");

// 16 um / 40x / 1.6 = 0.25 um --> scale 5x = 50 nm

direction 	= -1 	// 1 -> move up (cell is above) ; -1 -> move down (cell is below)
step 		= 5 	// step in pxl
distance 	= 100 	// distance in pxl

// Create both ROIs using the manual ROI as template

roiManager("Select", 0);
Roi.getBounds(x, y, width, height);
roiManager("translate", 0, height*direction*-1); // move rectangle up
roiManager("add");

for (p = 0; p < distance; p = p+step) {
wait(100);
roiManager("Select", 1);
roiManager("translate",0 , step*direction);
run("Measure");
};

nR = nResults;
I1 = newArray(nR);

// Grab the results for I1
for (i=0; i<nR;i++) {
	I1[i] = getResult("RawIntDen", i);
}

run("Clear Results");

for (p = 0; p < distance; p = p+step) {
wait(100);
roiManager("Select", 0);
roiManager("translate", 0, step*direction);
run("Measure");
};

nR = nResults;
I2 = newArray(nR);

// Grab the results for I2

for (i=0; i<nR;i++) {
	I2[i] = getResult("RawIntDen", i);
}

run("Clear Results");

// Make the new table
for (i=0; i<nR;i++) {
	setResult("I1", i, I1[i]);
	setResult("I2", i, I2[i]);
	setResult("I1-I2/I1+I2", i, (I1[i]-I2[i])/(I1[i]+I2[i]));
}
updateResults();


// Plot the curve

Plot.create("Plot of Results", "steps", "I1-I2/I1+I2");
Plot.add("Connected Circles", Table.getColumn("I1-I2/I1+I2", "Results"));
Plot.setStyle(0, "blue,#a0a0ff,1.0,Circle");
