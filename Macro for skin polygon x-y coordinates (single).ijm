//This macro returns all vertex points of the shape as approximated by polygons

//all variables in the code to be adjusted

//Set Image Scale in µm !!!!!!!!!!!!!!!! change appropriately !!!!!!!!!!!!!!!!
scale = 465.67; // this is for EXP2 and EXP4 Rag10,12,13,19,20,21,22,23,24,28,29,30,31,32,33
//scale = 927.02; // this is for EXP5 Rag14,15,16,17,18,25,26,27,34,35,36

//Marker line thickness
yW = 20; //white line; //full size = 20; half size = 10
lw1 = 1; // for measurement

//For Distance Measurement
//define earhole size to exclude abitrary skin tears / gaps
earHole = 30000000; //30,000,000 µm^2 = 30 mm^2 = 0.3 cm^2

//--------------------------------------------------------------------------------------------------------------------------------

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

//saving patch and density images and the results before closing
dir2 = dir1+"/../../OutputPolygon"+File.separator; //this creates the new path; "/../"selects parent folder
File.makeDirectory(dir2); //this creates the new path
path = dir2+imageName+"."+fileExt;//this creates the new path
dotIndex = lastIndexOf(path, "."); 
path = substring(path, 0, dotIndex); // remove ext.

//set scale
run("Set Scale...", "distance="+scale+" known=10000 pixel=1 unit=µm global");

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

//Getting raw image dimensions
w = getWidth;
h = getHeight;

//decides whether ear holes are present and line needs to be drawn
YesEar = getNumber("Are ear holes present? Yes = 0, No = 1", 0);

if (YesEar == 0) {
	//adjusts the line to measure
	run("Line Width...", "line="+lw1+"");
	
	//define spinal line
	waitForUser("Draw middle line along spin");
	getLine(x1, y1, x2, y2, lineWidth);
	Xs1 = x1;
}
//-----------------------------------------------	
//To get x,y coorcinates for the skin outline
selectWindow(image3);
run("Duplicate...", " ");
rename("Sink Outline");
run("RGB Color");
run("Colors...", "foreground=white background=black selection=yellow"); //adjust marker line colour
run("8-bit");
//setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(1, 255);
//setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(1, 255);
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
selectWindow("Results");
run("Close");
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
	run("Set Measurements...", "center display redirect=None decimal=2");
	run("Analyze Particles...", "size="+earHole+"-Infinity show=Masks display exclude");
	image4 = getTitle();
	selectWindow(image4);
	run("Invert");
	
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
	end
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
	selectWindow("Threshold");
	run("Close");
	selectWindow("Results");
	run("Close");
}

run("Collect Garbage");