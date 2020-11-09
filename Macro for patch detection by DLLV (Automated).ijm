//This macro measures patch distribution per area

//clear previously used RAM from previous runs
run("Collect Garbage");

//activates Bio-Formats Importer options
run("Bio-Formats Macro Extensions");

//get image directory and extract images
dir1 = getDirectory(" "); //change this according to file location
list = getFileList(dir1);

for (ii=0; ii<list.length; ii++) { //begin of outer loop
	path = dir1+list[ii];

//--------------------------------------------------------------------------------------------------------------------------------

//all variables in the code to be adjusted
//!!!!for some reason the macro works more relaibly with the variables within the loop to be redefined every loop around

//Variables that are called repeatedly in the loop

//Divisor that decides into how many parts the image is cut
Div = 10;

//Set Image Scale in µm !!!!!!!!!!!!!!!! change appropriately !!!!!!!!!!!!!!!!
//scale = 465.67; // this is for EXP2 and EXP4 Rag10,12,13,19,20,21,22,23,24,28,29,30,31,32,33
scale = 927.02; // this is for EXP5 Rag14,15,16,17,18,25,26,27,34,35,36
	
//Marker line thickness
yW = 20; //white outline line; //full size = 20; half size = 10
lw1 = 1; //line to measure lines
lw2 = 50; //white lines; //full size = 50; half size = 25 

//--------------------------------------------------------------------------------------------------------------------------------

    //opens file without opening Bio-Fomats importer window
    run("Bio-Formats Windowless Importer", "open=[path]");

	//saving patch and density images and the results before closing
	dir2 = dir1+"/../OutputDLLV"+File.separator; //this creates the new path; "/../"selects parent folder
	File.makeDirectory(dir2); //this creates the new path
	path = dir2+list[ii];//this creates the new path
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
	image1 = getTitle();
	run("Duplicate...", "title=[Cut up]");
	image2 = getTitle();
	run("Merge Channels...", "c1=["+Red+"] c2=["+Green+"] create");
	rename("Red Green Merge");
	
	//Managing Raw Image
	run("RGB Color");
	image3 = getTitle();
	run("Duplicate...", "title=[Patches by area]");
	image4 = getTitle();
	run("Duplicate...", "title=[Whole Skin]");
	image5 = getTitle();
		
	selectWindow(image3);
	
	//Getting raw image dimensions
	w = getWidth;
	h = getHeight;
	B = (w+1)/Div;
	
	//adjusts the line to measure
	run("Line Width...", "line="+lw1+"");
	
	//define spinal line
	waitForUser("Draw middle line along spin");
	getLine(x1, y1, x2, y2, lineWidth);
	Xs1 = x1;
	
	//print the variables
	print("w, h, B, Div, Xs1", "("+list[ii]+")");
	print(w, h, B, Div, Xs1);
	
	//Loop to measure back, flanks and abdomen		
	for (e = 0; e < 4; e += 1) { //begin of inner loop
	
	//--------------------------------------------------------------------------------------------------------------------------------
	
	//all variables in the code to be adjusted
	//!!!!for some reason the macro works more relaibly with the variables within the loop to be redefined every loop around
	
	//Variables that are called repeatedly in the loop
	
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
	
	//Patch size selection
	patchSize = 4148.65; //µm^2
	//according to Guertin & Sabatini (2015) [doi: 10.1038/npg.els.0003359], the average macrophage diameter is 20-30 µm;
	//if 25 µm is the mean diameter for macrophages, then one macrophage has an average 2D surface of ~490 µm^2
	//a pixel is in half size about 21.47 µm x 21.47 µm = 461 µm^2
	//Thus, 4 macrophages will require at least 3 x 3 = 9 pixels to be covered; Thus, 64.41 µm x 64.42 µm = 4148.65  µm^2
	//Since, the "scale" is adjusted for full sized images, the pixel number will dobule but the dimensions in µm will stay the same
		
	//Marker line thickness
	yY = 10; //yellow outline line; //full size = 10; half size = 5
	
	//--------------------------------------------------------------------------------------------------------------------------------
	
		//generate the count of rows in the results table for the if condition further down; this is needed in case no measurements are made
		if (isOpen("Results")) {
			numRow1 = nResults;		
		} else {
			numRow1 = 0;
		}
	
		selectWindow(image2);
	
		//select a region to work on
		if (e == 0) { //Back
			bp = "Back";
			makeRectangle(Xs1-B, 0, B*2, h);
			run("Cut");
			newImage(""+bp+"-1", "RGB black", B*2, h, 1);
			run("Paste");
			newImage(""+bp+"-2", "RGB black", w, h, 1);
			run("Add Image...", "image=["+bp+"-1] x="+Xs1-B+" y=0 opacity=100");
		}	
		else if (e == 1) { //Right Flanks (Sand Fly Exposed)
			bp = "Right Flank";
			makeRectangle(Xs1-B*3, 0, B*2, h);
			run("Cut");
			newImage(""+bp+"-1", "RGB black", B*2, h, 1);
			run("Paste");
			newImage(""+bp+"-2", "RGB black", w, h, 1);
			run("Add Image...", "image=["+bp+"-1] x="+Xs1-B*3+" y=0 opacity=100");
		}	
		else if (e == 2) { //Left Flanks (Unexposed)
			bp = "Left Flank";
			makeRectangle(Xs1+B, 0, B*2, h);
			run("Cut");
			newImage(""+bp+"-1", "RGB black", B*2, h, 1);
			run("Paste");
			newImage(""+bp+"-2", "RGB black", w, h, 1);
			run("Add Image...", "image=["+bp+"-1] x="+Xs1+B+" y=0 opacity=100");
		}	
		else { //Abdomen
			bp = "Abdomen";
			makeRectangle(Xs1-B*5, 0, B*Div, h);
			run("Cut");
			if (0 > Xs1-B*5) {
			newImage(""+bp+"-1", "RGB black", w-(w-(Xs1+B*5)), h, 1);
			}	else {
				newImage(""+bp+"-1", "RGB black", w-(Xs1-B*5), h, 1);
			}
			run("Paste");
			newImage(""+bp+"-2", "RGB black", w, h, 1);
			run("Add Image...", "image=["+bp+"-1] x="+Xs1-B*5+" y=0 opacity=100");
		}
	
		run("Flatten");
		rename(""+bp+" Area");
		image6 = getTitle();
	
		//Change name of image to be measured in for label adjustments
		selectWindow(image1);
		rename(""+bp+" Region");
		image1 = getTitle();
		
		if (e == 0) {
			run("Colors...", "foreground=white background=black selection=yellow"); //adjust marker line colour
		}	else if (e == 1) {
			run("Colors...", "foreground=white background=black selection=magenta"); //adjust marker line colour
		}	else if (e == 2) {
			run("Colors...", "foreground=white background=black selection=green"); //adjust marker line colour
		}	else {
			run("Colors...", "foreground=white background=black selection=cyan"); //adjust marker line colour
		}	
	
		selectWindow(image6);
		run("Duplicate...", "title=[Thresholding]");
		run("RGB Color");
		
		//identify patch areas
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
		for (i=0;i<3;i++){
		  selectWindow(""+i);
		  close();
		}
		selectWindow("Result of 0");
		close();
		selectWindow("Result of Result of 0");
		rename(a);
		// Colour Thresholding-------------end
		run("Set Measurements...", "area mean standard min centroid center perimeter fit skewness kurtosis display redirect=["+image1+"] decimal=2");
		run("Analyze Particles...", "size="+patchSize+"-Infinity show=Masks display");
		rename("Patch Mask");
		image7 = getTitle();
		
		//measure whole patch area per area
		if (isOpen("Results")) {//the following two if-conditions are supposed to help to get around the problem when no patches are detected and thus, no selection can be created
			numRow2 = nResults;
			
			if (0 < numRow2 - numRow1) {
				
				//measure whole patch area per area
				selectWindow(image7);
				run("Create Selection");
				run("Properties... ", "name=Patches width="+yY+"");
				run("ROI Manager...");
				roiManager("Add");
				roiManager("Select", 0);
					if (e == 0) {
							roiManager("Rename", "Whole Patch Area-1");
						}	else if (e == 1) {
							roiManager("Rename", "Whole Patch Area-2");
						}	else if (e == 2) {
							roiManager("Rename", "Whole Patch Area-3");		
						}	else {
							roiManager("Rename", "Whole Patch Area-4");
						}
				roiManager("Measure");
				selectWindow("ROI Manager");
				run("Close");
				selectWindow(image4);
				run("Restore Selection");
				run("Flatten");
				image4 = getTitle();
				selectWindow("Threshold Color");
				run("Close");
			}
		}
		close(image7);
		
		//measure whole skin per area
		selectWindow(image6);
		run("8-bit");
		//setAutoThreshold("Default dark");
		run("Threshold...");
		setThreshold(1, 255);
		run("Create Selection");
		run("Properties... ", "name=Patches width="+yY+"");
		run("ROI Manager...");
		roiManager("Add");
		roiManager("Select", 0);
			if (e == 0) {
					roiManager("Rename", "Whole Skin Area-1");
				}	else if (e == 1) {
					roiManager("Rename", "Whole Skin Area-2");
				}	else if (e == 2) {
					roiManager("Rename", "Whole Skin Area-3");
				}	else {
					roiManager("Rename", "Whole Skin Area-4");
				}
		roiManager("Measure");
		selectWindow("ROI Manager");
		run("Close");
		selectWindow("Threshold");
		run("Close");
		
		//Closing Windows
		close("Mask of "+bp+" Area-1");
		close("Mask of "+bp+" Area");
		close(""+bp+" Area-1");
		close(""+bp+" Area");
		close(""+bp+"-2");
		close(""+bp+"-1");
		run("Collect Garbage");
	} //end of inner loop
	
	rename("Patches by Area");
	image8 = getTitle();
	
	//Creat skin outline
	selectWindow(image5);
	run("Colors...", "foreground=white background=black selection=white"); //adjust marker line colour
	run("8-bit");
	//setAutoThreshold("Default dark");
	selectWindow(image5);
	run("Threshold...");
	setThreshold(1, 255);
	run("Create Selection");
	selectWindow(image8);
	run("Restore Selection");
	run("Properties... ", "  width="+yW);
	run("Flatten");
	image9 = getTitle();
	selectWindow("Threshold");
	run("Close");
	
	//make grid lines
	run("Line Width...", "line="+lw2+"");
	selectWindow(image9);
	makeLine(Xs1-B*3, 0, Xs1-B*3, h);
	run("Flatten");
	makeLine(Xs1-B, 0, Xs1-B, h);
	run("Flatten");
	makeLine(Xs1+B, 0, Xs1+B, h);
	run("Flatten");
	makeLine(Xs1+B*3, 0, Xs1+B*3, h);
	run("Flatten");
		
	//Close all unessential remaining images and clear RAM
	close("\\Others");
	run("Collect Garbage");
	
	//saving
	save(path+"_by_area.tif");
	close();
	saveAs("Results", path+"_Results.csv");
	selectWindow("Results");
	run("Close");
	
	//clear previously used RAM from previous runs
	run("Collect Garbage");
} //end of outer loop

//save the log sheet with the variables
selectWindow("Log");
saveAs("Text", path+"_DLLV_variables_All.txt");
selectWindow("Log");
run("Close");
