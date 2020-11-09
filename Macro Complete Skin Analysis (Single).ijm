//#####################################################################################################################################

// Title:        Skin patch analysis

// Script Name:  Macro Complete Skin Analysis (Single).ijm
// Copyright:    Johannes S. P. Doehl (NIAID) & Paul M. Kaye (University of York)

// Purpose:      Analysis of visceralizing Leishmania-infected whole-body skin images for the extraction of the following data:
//				 1 - extraction of x & y coordinates of the skin ouline and ear wholes, if present
//				 2 - detection of skin patches, measurement of patch size and density and nearest enighbour analysis
//				 3 - distance measurements of patches from the right ear

// Author:       Johannes S. P. Doehl
// Date:         March 2020

// Publication:

//#####################################################################################################################################

//ADJUSTABLE VARIABLES

//Set Image Scale in µm !!!!!!!!!!!!!!!! change appropriately !!!!!!!!!!!!!!!!
scale = 465.67; // adjust according to your image resolution by measuring scale bar
//scale = 927.02; // adjust according to your image resolution by measuring scale bar

//Marker line thickness
lw1 = 1; // for measurement
yY = 5; //yellow line; //full size = 10; half size = 5
yR = 10; //red line; //full size = 20; half size = 10
yW = 10; //white line; //full size = 20; half size = 10
yM = 5; //magenta line; //full size = 10; half size = 5
lineWidth = 40; //// line thickness of distance measurement line

//label dimesions for label in density image
zoomMag = 7; //full size = 12; half size = 7

//For Distance Measurement
//define earhole size to exclude abitrary skin tears / gaps
earHole = 30000000; //30,000,000 µm^2 = 30 mm^2 = 0.3 cm^2

//Color Threshold for patch selection
//Hue
Hmin = 0;
Hmax = 255;
//Saturation
Smin = 0;
Smax = 255;
//Brightness
Bmin = 70;
Bmax = 255;

//Mimimum size that describes a parasite patch
patchSize = 4148.65; //µm^2
//according to Guertin & Sabatini (2015) [doi: 10.1038/npg.els.0003359], the average macrophage diameter is 20-30 µm;
//if 25 µm is the mean diameter for macrophages, then one macrophage has an average 2D surface of ~490 µm^2
//a pixel is in half size about 21.47 µm x 21.47 µm = 461 µm^2
//Thus, 4 macrophages will require at least 3 x 3 = 9 pixels to be covered; Thus, 64.41 µm x 64.42 µm = 4148.65  µm^2
//Since, the "scale" is adjusted for full sized images, the pixel number will dobule but the dimensions in µm will stay the same

//Max distance allowed between patches for cluster analysis
maxDis = 3900; // µm
//Lutz & Neiva (1912) decibed Lutzomyia longipalpis female adult sand flies as about 2 mm full body length, a wingspan around 2.3 mm and a body with of about 6.5 mm
//Thus, if one assumes that one body length into any direction from where a sand fly lands on the skin can be considered a foraging area, then patches in a cluster should be max about 4 mm apart (edge to edge)
//Since patches will have to reach into this radius and the labrium is about 35 µm in diameter according to Brinson (1993)than from any side the patch needs to reach into the foragin area by at least 50 µm
//Thus, 4000 µm - 2x50 µm = 3900 µm diameter or max patch separation (edge to edge)

//define last line for the loop
loopEnd = 5;

//#####################################################################################################################################

//SCRIPT

//Define directory from where to upload images

//clear previously used RAM from previous runs
run("Collect Garbage");

//activates Bio-Formats Importer options
run("Bio-Formats Macro Extensions");

///gets image name, extention and directory
imagePath = File.openDialog(" "); //Full file path with file name and extention
imageName = File.getName(imagePath); //extracts file name with extention
dotIndex = lastIndexOf(imageName, "."); //index for the . between file name and extention
fileLength = lengthOf(imageName); //gets the total length of the file name.
fileExt = substring(imageName, (dotIndex+1), fileLength); //gets the extention
imageName = substring(imageName, 0, dotIndex); //gets the file name without extention
pathLength = lengthOf(imagePath); //gets the total length of the whole file path
dir1 = substring(imagePath, 0, (pathLength-fileLength)); //gets the directory of the image path

//opens file without opening Bio-Fomats importer window
run("Bio-Formats Windowless Importer", "open=["+dir1+"/"+imageName+"."+fileExt+"]");

//set scale
run("Set Scale...", "distance="+scale+" known=10000 pixel=1 unit=µm global");

//extracts the pixel size in µm depending on the set scale
getPixelSize(unit, pw, ph, pd);
conv = pw

//split colors to get red only
run("Split Channels");
Blue = getTitle();
close(Blue);
Green=getTitle();
run("Put Behind [tab]");
Red=getTitle();
run("Duplicate...", "title=[Area Measurements]");
image2 = getTitle();
run("Merge Channels...", "c1=["+Red+"] c2=["+Green+"] create");
rename("Red Green Merge");
run("Flatten");
image3 = getTitle();
close("Red Green Merge");

//select starting image
selectWindow(image3);

//Getting raw image dimensions
w = getWidth;
h = getHeight;

//decides whether ear holes are present and line needs to be drawn
YesEar = getNumber("Are ear holes present? Yes = 0, No = 1", 0);

if (YesEar == 0) {
	//adjusts the line to measure
	run("Line Width...", "line="+lw1+"");
	
	//define spinal line
	waitForUser("Draw middle line along spine");
	getLine(x1, y1, x2, y2, lineWidth);
	Xs1 = x1;
	run("Select None");
}

//-------------------------------------------------------------------------------------------------------------------------------------

//saving patch and density images and the results before closing
dir2 = dir1+"/../../OutputWhole"+File.separator; //this creates the new path; "/../"selects parent folder
File.makeDirectory(dir2); //this creates the new path
path = dir2+imageName+"."+fileExt;//this creates the new path
dotIndex = lastIndexOf(path, "."); 
path = substring(path, 0, dotIndex); // remove ext.

selectWindow(image2);
setBatchMode(true);
run("Duplicate...", " ");
run("RGB Color");	

//analysis for patches
run("Colors...", "foreground=white background=black selection=yellow"); //adjust marker line colour
run("Color Threshold...");
// Color Thresholder 2.0.0-rc-69/1.52n
min=newArray(3);
max=newArray(3);
filter=newArray(3);
a=getTitle();
run("HSB Stack");
run("Convert Stack to Images");
selectWindow("Hue");
rename("0");
selectWindow("Saturation");
rename("1");
selectWindow("Brightness");
rename("2");
min[0]=Hmin;
max[0]=Hmax;
filter[0]="pass";
min[1]=Smin;
max[1]=Smax;
filter[1]="pass";
min[2]=Bmin;
max[2]=Bmax;
filter[2]="pass";
for (i=0;i<3;i++){
  selectWindow(""+i);
  setThreshold(min[i], max[i]);
  run("Convert to Mask");
  if (filter[i]=="stop")  run("Invert");
}
imageCalculator("AND create", "0","1");
imageCalculator("AND create", "Result of 0","2");
for (i=0;i<3;i++) {
  selectWindow(""+i);
  close();
}
selectWindow("Result of 0");
close();
selectWindow("Result of Result of 0");
rename(a);
// Colour Thresholding-------------end
setBatchMode(false);
//Measuring patch areas mean gray values and centre of mass
run("Set Measurements...", "area mean standard min centroid center perimeter fit feret's skewness kurtosis stack display redirect=[Area Measurements] decimal=2");
run("Analyze Particles...", "size="+patchSize+"-Infinity show=[Masks] display exclude record add"); //measure selected areas by color threshold
rename("NND");
image4 = getTitle();
maxCount = nResults;
print("Patch Counts");
print(maxCount);
selectWindow("Log");
saveAs("Text", path+"_patch_counts.csv");
run("Close");
selectWindow("Threshold Color"); //select and close Color Threshold Window
run("Close");

//-------------------------------------------

//generate image density images
selectWindow(image2);
rename("Density Image");
run("royal");
run("Calibration Bar...", "location=[Upper Right] fill=White label=Black number=5 decimal=0 font=12 zoom="+zoomMag+" bold overlay");
run("Flatten");
image5 = getTitle();
save(path+"_density.tif");
close("Density Image");
run("Collect Garbage");
setBatchMode(true);
//-------------------------------------------

if (isOpen("Results")) {//the following two if-conditions are supposed to help to get around the problem when no patches are detected and thus, no selection can be created

	//Nearest Neighbour Distance (NND) and average distance to nearest 6 neichbors measurements
	//Haeri & Haeri (2015), 3, e1
	selectWindow(image4);
	run("ND ");
	selectWindow("Distance Between Neighboring Particles");
	saveAs("Results", path+"_NND.csv");
	run("Close");
	run("Collect Garbage");

	//-------------------------------------------
	
	//Clustering patches by max distance (Graph)
	//Ben Tupper
	selectWindow(image4);
	run("Duplicate...", "title=[Cluster]");
	run("Create Selection");
	run("Colors...", "foreground=white background=black selection=red"); //adjust marker line colour
	run("Line Width...", "line="+yR); //modifies boarder line width
	run("Graph ", "mask=Cluster distance=Edges neighbors="+maxDis+" lines label matrix");
	selectWindow("Distance Matrix");
	saveAs("Results", path+"_Distance_Matrix.csv");
	run("Close");
	selectWindow("Log");
	saveAs("Text", path+"_Clusters.csv");
	run("Close");
	selectWindow("Adjacencies");
	run("Flatten");
	close("Adjacencies");
	close("Cluster");
	run("Collect Garbage");
	setBatchMode(false);
	//-------------------------------------------

	//get whole skin boarders
	selectWindow(image3);
	run("Duplicate...", "title=[Cluster Boarder]");
	run("Colors...", "foreground=white background=black selection=black"); //adjust marker line colour
	run("8-bit");
	//setAutoThreshold("Default dark");
	run("Threshold...");
	setThreshold(1, 255);
	run("Create Selection");
	selectWindow("Adjacencies-1");
	run("Restore Selection");
	run("Properties... ", "  width="+yW);
	run("Flatten");
	save(path+"_distance.tif");
	close("Area Measurements-1");
	close("Cluster Boarder");
	close("Adjacencies-1");
	close("Adjacencies-2");
	run("Collect Garbage");
	
	//get complete patch coverage
	selectWindow(image4);
	run("Create Selection");
	run("Colors...", "foreground=white background=black selection=yellow"); //adjust marker line colour
	run("Set Measurements...", "area mean standard min centroid center perimeter fit skewness kurtosis stack display redirect=None decimal=2");
	selectWindow(image3);
	run("Duplicate...", "title=[Patch Boarders]");
	image6=getTitle(); //loads image with patch boarders into variable
	run("Restore Selection");
	run("Properties... ", "  width="+yY); //modifies boarder line width
	run("ROI Manager...");
	roiManager("Add");
	roiManager("Select", 0);
	roiManager("Rename", "Sum Patch Area");
	roiManager("Measure");
				
	//create image with marked boarders for measured patches
	selectWindow(image6);
	run("Flatten");
	run("Select None");
	run("Restore Selection");
	run("Flatten");
	image6 = getTitle();
	selectWindow(image5);
	run("Colors...", "foreground=white background=black selection=magenta"); //adjust marker line colour
	run("Restore Selection");
	run("Properties... ", "  width="+yM); //modifies boarder line width
	run("Flatten");
	image5 = getTitle();
	selectWindow("ROI Manager");
	run("Close");
	run("Collect Garbage");
}

//-------------------------------------------

//get whole skin area and skin boarders
selectWindow(image3);
run("Duplicate...", "title=[Whole Skin]");
image7=getTitle(); //loads image with patch boarders into variable
run("Duplicate...", "title=[Skin Boarder]");
run("Colors...", "foreground=white background=black selection=white"); //adjust marker line colour
run("RGB Color");
run("8-bit");
//setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(1, 255);
run("Create Selection");
run("ROI Manager...");
roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "Skin Area");
roiManager("Measure");
run("Close");
selectWindow(image5);
run("Restore Selection");
run("Properties... ", "  width="+yW);
run("Flatten");
image5=getTitle(); //loads image with patch boarders into variable
selectWindow(image6);
run("Restore Selection");
run("Properties... ", "  width="+yW);
run("Flatten");
image6=getTitle(); //loads image with patch boarders into variable

//close windows not needed
selectWindow("Threshold");
run("Close");
close("Whole Skin");
close("Skin Boarder");
close("Red Green Merge (RGB)-1");
close("Red Green Merge (RGB)-2");
close("NND");
close("Density Image-1");
close("Density Image-2");
close("Patch Boarders");
close("Patch Boarders-1");
close("Patch Boarders-2");
run("Collect Garbage");

//-------------------------------------------

//save
selectWindow(image5);
save(path+"_densityboarder.tif");
close();
selectWindow(image6);
save(path+"_patches.tif");
selectWindow("Results");
saveAs("Results", path+"_Results.csv");

//-------------------------------------------------------------------------------------------------------------------------------------

//Extracts x & y coordinates of the skin outline and the ear holes

//saving patch and density images and the results before closing
dir2 = dir1+"/../../OutputPolygon"+File.separator; //this creates the new path; "/../"selects parent folder
File.makeDirectory(dir2); //this creates the new path
path = dir2+imageName+"."+fileExt;//this creates the new path
dotIndex = lastIndexOf(path, "."); 
path = substring(path, 0, dotIndex); // remove ext.

//-----------------------------------------------

//To get x,y coorcinates for the skin outline
selectWindow(image3);
run("Duplicate...", "title=[Sink Outline]");
run("Colors...", "foreground=white background=black selection=yellow"); //adjust marker line colour
run("8-bit");
//setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(1, 255);
//setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(1, 255);
run("Set Measurements...", "area mean standard min centroid center perimeter fit skewness kurtosis stack display redirect=None decimal=2");
run("Analyze Particles...", "size=4148.65-Infinity show=Masks display include");
run("Create Selection");
run("Properties... ", "  width="+yW);

"List XY Coordinates" { //This will get all the xy coordinates for the outline
	requires("1.30k");
	getSelectionCoordinates(x, y);
	for (i=0; i<x.length; i++)
		print(i+" "+x[i]+" "+y[i]);
}

//closes all open windows not required anymore
close("Sink Outline");
close("Mask of Sink Outline");
selectWindow("Threshold");
run("Close");
//selectWindow("Results");
//run("Close");
selectWindow("Log");
saveAs("Text", path+"_SkinOutline.csv");
selectWindow("Log");
run("Close");
run("Collect Garbage");

//-------------------------------------------
	
//to get the x,y coordinated for the ears
if (YesEar == 0) {
	selectWindow(image3);
	run("Duplicate...", "title=[Ear Holes]");
	//find the ear holes in the skin
	run("8-bit");
	//setAutoThreshold("Default dark");
	run("Threshold...");
	setThreshold(0, 1);
	run("Set Measurements...", "area mean standard min centroid center perimeter fit skewness kurtosis stack display redirect=None decimal=2");
	run("Analyze Particles...", "size="+earHole+"-Infinity show=Masks display exclude");
	image4 = getTitle();
	selectWindow(image4);
	run("Invert");
	setBatchMode(true);
	for (e = 0; e < 2; e += 1) { //This cuts up the image in to a right and left side with the help of the line to analyse ear holes separately
		selectWindow(image4);
		if (e == 0) { //Right ear on the left image side
			bp = "Right Ear";
			makeRectangle(0, 0, Xs1, h);
			run("Cut");
			newImage(""+bp+"-1", "8-bit white", Xs1, h, 1);
			run("Paste");
			newImage(""+bp+"-2", "8-bit white", w, h, 1);
			run("Add Image...", "image=["+bp+"-1] x=0 y=0 opacity=100");
		}	
		else if (e == 1) { //Left ear is on the right image side
			bp = "Left Ear";
			makeRectangle(Xs1, 0, w, h);
			run("Cut");
			newImage(""+bp+"-1", "8-bit white", w-Xs1, h, 1);
			run("Paste");
			newImage(""+bp+"-2", "8-bit white", w, h, 1);
			run("Add Image...", "image=["+bp+"-1] x="+Xs1+" y=0 opacity=100");
		}	
		run("Flatten");
		setBatchMode(false);
		run("8-bit");
		run("Create Selection");
		
		"List XY Coordinates" { //This will get all the xy coordinates for the outline
			requires("1.30k");
			getSelectionCoordinates(x, y);
			for (i=0; i<x.length; i++)
				print(i+" "+x[i]+" "+y[i]);
		}
		
		if (e == 0) { //chooses to either safe the right or left ear data
			selectWindow("Log");
			saveAs("Text", path+"_REarOutline.csv");
			selectWindow("Log");
			run("Close");
		} else if (e == 1) {
			selectWindow("Log");
			saveAs("Text", path+"_LEarOutline.csv");
			selectWindow("Log");
			run("Close");
		}	
	}
	
	//closes all open windoes
	close("Right Ear-1");
	close("Right Ear-2");
	close("Right Ear-3");
	close("Left Ear-1");
	close("Left Ear-2");
	close("Left Ear-3");
	close("Ear Holes");
	close("Mask of Ear Holes");
	close("Ear Holes");
	close("Red Green Merge (RGB)");
	selectWindow("Threshold");
	run("Close");
}

run("Collect Garbage");

//-------------------------------------------------------------------------------------------------------------------------------------

//saving patch and density images and the results before closing
dir2 = dir1+"/../../OutputDistance"+File.separator; //this creates the new path; "/../"selects parent folder
File.makeDirectory(dir2); //this creates the new path
path = dir2+imageName+"."+fileExt;//this creates the new path
dotIndex = lastIndexOf(path, "."); 
path = substring(path, 0, dotIndex); // remove ext.

//select image for distance measurement and close all else to recover RAM
selectWindow(image6);
close("\\Others");
run("Collect Garbage");

// open a sample image
mouseEar = nResults - 2;
selectWindow(image6);
rename("Distance");
run("Colors...", "foreground=white background=black selection=green"); //adjust marker line colour
run("Set Measurements...", "center display redirect=None decimal=2");

//run loop without opening new images
setBatchMode(true);

for (i = 1; i <= maxCount; i++) {
	//get coordinates from results table
	originX = getResult("XM", mouseEar);
	originY = getResult("YM", mouseEar);
	destinX = getResult("XM", 1 + i);
	destinY = getResult("YM", 1 + i);
	// make an example selection
	xpoints = newArray(originX/conv, destinX/conv);
	ypoints = newArray(originY/conv, destinY/conv);
	makeSelection("point", xpoints, ypoints);
	//measures the length of the drawn line	
	Roi.getCoordinates(xpoints, ypoints); // store point coordinates
	run("Line Width...", "line="+lineWidth); //adjusts line thickness
	makeLine(xpoints[0],ypoints[0],xpoints[1],ypoints[1]); // make line
	run("Measure"); // measure the line
	run("Flatten");
	close("\\Others");
	run("Collect Garbage");
}
setBatchMode(false);

//saving distance image and results
save(path+"_distance.tif");
close();
selectWindow("Results");
IJ.deleteRows(0, maxCount + 4);
saveAs("Results", path+"_distance.csv");
selectWindow("Results");
run("Close");

//-------------------------------------------------------------------------------------------------------------------------------------

//clear previously used RAM from previous runs
run("Collect Garbage");		
