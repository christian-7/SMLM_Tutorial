run("Clear Results");

// Create both ROIs

roiManager("Select", 0);
Roi.getBounds(x, y, width, height);
roiManager("translate", 0, height); // move rectangle up
roiManager("add");

// Print out their positions

roiManager("Select", 0);
Roi.getBounds(x, y, width, height);
print(x, y);
roiManager("Select", 1);
Roi.getBounds(x, y, width, height);
print(x, y);

// Start measurement

roiManager("Select", 0);
roiManager("multi-measure measure_all");

nR = nResults;
I1 = newArray(nR);

// Grab the results for I1
for (i=0; i<nR;i++) {
	I1[i] = getResult("RawIntDen", i);
}

run("Clear Results");
roiManager("Select", 1);
roiManager("multi-measure measure_all");

nR = nResults;
I2 = newArray(nR);

// Grab the results for I2

for (i=0; i<nR;i++) {
	I2[i] = getResult("RawIntDen", i);
}

// Keep the old table
IJ.renameResults("Raw Results");

// Make the new table
for (i=0; i<nR;i++) {
	setResult("I1", i, I1[i]);
	setResult("I2", i, I2[i]);
	setResult("I1-I2/I1+I2", i, (I1[i]-I2[i])/(I1[i]+I2[i]));
}
updateResults();


// Plot the curve

Plot.create("Plot of Results", "frames", "I1-I2/I1+I2");
Plot.add("Circle", Table.getColumn("I1-I2/I1+I2", "Results"));
Plot.setStyle(0, "blue,#a0a0ff,1.0,Circle");