/* This macro extracts the 100um edge outside a scar
 *  It saves Channels 1, 3 and 4 and processes channel 1 through StarDist
 *  The files are saved in a foder for analysis in CellProfiler
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
outputDir=getDirectory("And the output folder");
myList=getFileList(inDir);  //an array

//Define your measurments and settings for 
run("Set Measurements...", "area mean display redirect=None decimal=0");
roiManager("Set Line Width", 1);

for (j = 0 ; j < myList.length ; j++ ){
	path=inDir+myList[j];   //path to each file
	run("Bio-Formats Importer", "open=["+path+"] autoscale color_mode=Composite view=Hyperstack stack_order=XYCZT series_"+(1));
	FileName=File.nameWithoutExtension;
	ImageID=File.name;
	print("Processing "+ImageID);

	//first defining the wound area
	roiManager("reset");
	setTool("freehand");
	waitForUser("Defining wound area", "Please define the wound area from the edge of the heart\nPress OK when done");
	roiManager("Add");
	//Enlarging area and substracting original area.  The selected region contains the cells of interest.
	run("Enlarge...", "enlarge=100");
	roiManager("Add");
	roiManager("Select", 0);
	setBackgroundColor(0, 0, 0);
	run("Clear", "stack");
	roiManager("Select", 1);
	run("Clear Outside", "stack");
	run("Duplicate...", "title=ROI duplicate");
	rename("ROI");
	close(ImageID);
	selectWindow("ROI");
	roiManager("reset");
	
	//The defined ROI contains a lot of blank space, wel'll get rid of it using C2
	Stack.setChannel(2);
	run("Duplicate...", "title=Heart channels=2");
	setThreshold(100, 65535);
	setOption("BlackBackground", false);
	run("Convert to Mask");
	run("Median...", "radius=2");
	run("Fill Holes");
	run("Set Measurements...", "area mean min display redirect=None decimal=0");
	run("Analyze Particles...", "size=100-Infinity add");
	
	//we have created one but possibly also a series of ROIs, the largest of which is the heart region
	NumberOfROIs=roiManager("count");
	BiggestROI=0;
	ROIIndex=0;
	for (ROI = 0; ROI < NumberOfROIs; ROI++) {
		ROIsize=getResult("Area",ROI);
		if(ROIsize>BiggestROI){
			BiggestROI=ROIsize;
			ROIIndex=ROI;
		}
	}
	close("Results");
	ROIsToDelete=Array.getSequence(NumberOfROIs);
	ROIsToDelete=Array.deleteIndex(ROIsToDelete, ROIIndex);
	roiManager("select", ROIsToDelete);

	//roiManager("delete");
	
	roiManager("Select", 0);
	selectWindow("ROI");
	roiManager("Select", 0);
	run("Crop");
	close("Heart");
	
	selectWindow("ROI");
	run("Select None");
	run("Split Channels");

	selectWindow("C4-ROI");
	saveAs("Tif", outputDir+"Mef2-"+FileName);
	close();
	selectWindow("C3-ROI");
	saveAs("Tif", outputDir+"PCNA-"+FileName);
	close();
	selectWindow("C2-ROI");
	close();
	selectWindow("C1-ROI");
	run("Duplicate...", "title=StarDist");
	getHistogram(values, counts, 256);
	values2=Array.deleteIndex(values, 0);
	counts2=Array.deleteIndex(counts, 0);
	Array.getStatistics(counts2, min, max, mean, stdDev);
	indexValue=0;
	for (i = 0; i < counts2.length-1; i++) {
		if(counts2[i]>indexValue){
			indexValue=counts2[i];
			MaxIndex=i;
		}
	}
	BackgroundCorrection=values2[MaxIndex];
	run("Subtract...", "value="+BackgroundCorrection);
	getDimensions(width, height, channels, slices, frames);
	WidthTile=round(width/750);
	HeightTile=round(height/750);
	Tiles=WidthTile*HeightTile;
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'StarDist', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.479071', 'nmsThresh':'0.3', 'outputType':'Label Image', 'nTiles':'"+Tiles+"', 'excludeBoundary':'3', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	selectWindow("C1-ROI");
	saveAs("Tif", outputDir+"DAPI-"+FileName);
	close();
	selectWindow("Label Image");
	saveAs("Tif", outputDir+"StarDist-"+FileName);
	close("*");
}

print("**** MACRO DONE ****\nAll Regions processed");
	