// roiManager("Select", 1); // select the line
//Roi.getCoordinates(linex, liney);
//Roi.getContainedPoints(xpoints, ypoints);
//Array.show(linex);
//Array.show(xpoints);



roiManager("Select", 1);

step = 1

for (p = 0; p < xpoints.length; p = p+step)

roiManager("Select", 1);
roiManager("translate", step, 1);
run("Measure");
