/* This first macro pulls the scar tissue in the .czi files and saves them in a folder to be analysed later.
 */
 
//CLEAR LOG
print("\\Clear");

// CLOSE ALL OPEN IMAGES
while (nImages>0) { 
	selectImage(nImages); 
    close(); 
}

//SET BATCH MODE 
setBatchMode(false);   //true or false, true if you don't want to see the images, which is faster

//START MESSAGE
print("**** STARTING THE MACRO ****");

//INPUT/OUPUT folders
inDir=getDirectory("Choose the input folder"); 
WoundOutputDir=getDirectory("Where should I save the wound picture?");
HeartOutputDir=getDirectory("And the whole heart pictures?");
MuscleOutputDir=getDirectory("Finally, in which folder should I save images of muscle only?");
myList=getFileList(inDir);  //an array

//Define your measurments and settings for 
run("Set Measurements...", "mean display redirect=None decimal=3");
run("Roi Defaults...", "color=green stroke=1 group=0");
HeartCounter=0;
ScarredHeartCounter=0;


for (j = 0 ; j < myList.length ; j++ ){
	path=inDir+myList[j];   //path to each .czi file
	run("Bio-Formats Macro Extensions");
	Ext.setId(path);
	Ext.getSeriesCount(seriesCount);
	/*the .czi files contains a number of pictures in pyramidal format (i.e. the same image at different levels of resolution).  
	 * These are unfortunately not regular: there are 3-5 images of the same region at different level of resolution.  
	 * I therefore need to find out how many regions have been photographed and how many layers each has so I only open a single image per region
	 */
	run("Bio-Formats Importer", "open=["+path+"] color_mode=Default view=Hyperstack stack_order=XYCZT series_"+(1));
	FileName=File.nameWithoutExtension;
	ImageID=File.name;
	ID=getImageID();
	print("Processing "+ImageID+" series 1");

	//finding the number of regions
	FileMetadata=getImageInfo();
	ScanRegion=lastIndexOf(FileMetadata, "ScanRegion");
	NumberOfRegions=substring(FileMetadata, ScanRegion+10, ScanRegion+11);
	NumberOfRegions=parseInt(NumberOfRegions);
	print("There are "+NumberOfRegions+1+" regions in "+myList[j]);

	//Finding the highest resolution layer for analysis
	LayerCount=1;
	ChosenLayersPerRegion=newArray(0);
	ChosenLayersPerRegion=Array.concat(ChosenLayersPerRegion,1);
	
	for (region = 0; region < NumberOfRegions+1; region++) {
		PyramidLayersCount=lastIndexOf(FileMetadata, "PyramidLayersCount #"+region+1);
		RegionLayers=substring(FileMetadata, PyramidLayersCount+24, PyramidLayersCount+25);
		RegionLayers=parseInt(RegionLayers)+1;
		LayerCount=LayerCount+RegionLayers;
		ChosenLayersPerRegion=Array.concat(ChosenLayersPerRegion,LayerCount);
		}
	ChosenLayersPerRegion=Array.deleteIndex(ChosenLayersPerRegion, NumberOfRegions+1);
	print("I will analyse the following series:");
	Array.print(ChosenLayersPerRegion);
	close();

	for (i=0; i<NumberOfRegions+1; i++) {
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default view=Hyperstack stack_order=XYCZT series_"+(ChosenLayersPerRegion[i]));
		FileName=File.nameWithoutExtension;
		ImageID=File.name;
		ID=getImageID();
		ScarExtraction(ID);
	}
	HeartCounter=HeartCounter+NumberOfRegions+1;
	
	
}

close("Results");

//saving log
print("There were "+HeartCounter+" hearts, of which "+ScarredHeartCounter+" had a scar.");
selectWindow("Log");
saveAs("Text", HeartOutputDir+FileName+"_Log.txt");
close("*");

print("***** Macro done *****");

function ScarExtraction (ID) {
	run("RGB Color");
	rename("RGB");
	selectWindow("RGB");
	Scar=getBoolean("Is there a scar in this heart?");

	//extracting the scar for Ilastik analysis
	if(Scar==1){
		ScarredHeartCounter=ScarredHeartCounter+1;
		setTool("zoom");
		waitForUser("Please zoom in on the wounding site.\nWhen done, click OK.");
		setTool("rectangle");
		run("Roi Defaults...", "color=green stroke=2 group=0");
		waitForUser("Please define the wounding site.\nEnsure you enclude the pericardium by the scar.\nWhen done, click OK.");
		run("Add Selection...");
		run("Duplicate...", "title=Wound");
		run("Measure");
		Area=getResult("Area");
		setMetadata("Info", Area);
		PaddedLayer=IJ.pad(ChosenLayersPerRegion[i], 2);
		saveAs("Tiff", WoundOutputDir+FileName+"_series_"+PaddedLayer+"_Wound");
	}
	else{
		Area=0;
	}

	//extracting muscle for general analysis
	selectWindow("RGB");
	run("Scale to Fit");
	setTool("rectangle");
	run("Roi Defaults...", "color=magenta stroke=2 group=0");
	makeRectangle(3084, 2580, 1500, 1500);
	roiManager("Set Color", "red");
	roiManager("Set Line Width", 2);
	waitForUser("Defining Muscle Region", "Please move the rectangle to a region of the heart that is pure muscle\nWhen done, click OK.");
	run("Add Selection...");
	run("Duplicate...", "title=muscle duplicate");
	PaddedLayer=IJ.pad(ChosenLayersPerRegion[i], 2);
	saveAs("Tiff", MuscleOutputDir+FileName+"_series_"+PaddedLayer+"_Muscle");

	//saving whole image with defined area highlighted (red for muscle, green for scar, if present)
	selectWindow("RGB");
	setMetadata("Scar_Area=", Area);
	saveAs("Tiff", HeartOutputDir+FileName+"_series_"+PaddedLayer+"_Heart");
	close("*");
}
