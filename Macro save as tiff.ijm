//This macro returns all vertex points of the shape as approximated by polygons

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

//Set Image Scale in µm !!!!!!!!!!!!!!!! change appropriately !!!!!!!!!!!!!!!!
scale = 465.67; // this is for EXP2 and EXP4 Rag10,12,13,19,20,21,22,23,24,28,29,30,31,32,33
//scale = 927.02; // this is for EXP5 Rag14,15,16,17,18,25,26,27,34,35,36
	
//--------------------------------------------------------------------------------------------------------------------------------

    //opens file without opening Bio-Fomats importer window
    run("Bio-Formats Windowless Importer", "open=[path]");

	//saving patch and density images and the results before closing
	dir2 = dir1+"/../../OutputTIFF"+File.separator; //this creates the new path; "/../"selects parent folder
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
	
	run("Merge Channels...", "c1=["+Red+"] c2=["+Green+"] create");
	rename("Red Green Merge");
	run("Flatten");
	image3 = getTitle();

	selectWindow(image3);
	save(path+".tif");
	close();
	run("Close All");
	run("Collect Garbage");
}
	