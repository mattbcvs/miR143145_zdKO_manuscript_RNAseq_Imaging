/* This macro quantifies the scar tissue composition from the Ilastik pixel segmentation.
 *  It measures the scar area and its components
 *  It projects outlines of the components on pictures of the scar tissue, for quality control
 */
 
//CLEAR LOG
print("\\Clear");

// CLOSE ALL OPEN IMAGES
while (nImages>0) { 
	selectImage(nImages); 
    close(); 
}

//SET BATCH MODE 
setBatchMode(true);   //true or false, true if you don't want to see the images, which is faster

//START MESSAGE
print("**** STARTING THE MACRO ****");

//INPUT/OUPUT folders
inDir=getDirectory("Choose the input folder"); 
OutputDir=getDirectory("Where should I save the results?");
myList=getFileList(inDir);  //an array

//Define your measurments and settings for 
run("Set Measurements...", "mean display redirect=None decimal=3");
run("Roi Defaults...", "color=green stroke=1 group=0");

FileNames=newArray();
Area=newArray();
rowIndex=0;

//opening the files and their corresponding simple segmentation maps
for (i = 0; i < myList.length; i++) {
	path=inDir+myList[i];
	if(endsWith(path, "Wound.tif")){
		open(path);
		CompleteFileName=File.name;
		FileName=File.nameWithoutExtension;
		print("Analysing "+FileName);
		FileNames=Array.concat(FileNames,FileName);
		run("Duplicate...", "title=Flat");
		open(inDir+FileName+"_Simple Segmentation.tif");
		rename("Segmentation");
		run("glasbey_on_dark");

		//selecting individual dyes
		Dyes=newArray("Collagen", "Muscle","Fibrin");
		DyeValues=newArray(4,5,6);
		DyesHighlight=newArray("cyan", "yellow", "magenta");

		//calculating area of individual dyes and applying QC
		for (d = 0; d < 3; d++) {
			selectWindow("Segmentation");
			run("Duplicate...", "title="+Dyes[d]);
			setThreshold(DyeValues[d], DyeValues[d]);
			setOption("BlackBackground", false);
			run("Convert to Mask");
			run("Analyze Particles...", "summarize add");

			//collecting data from summary
			TotalArea=Table.get("Total Area", rowIndex, "Summary");
			Area=Array.concat(Area,TotalArea);
			rowIndex=rowIndex+1;

			//applying ROIs as overlay to Wound image for QC
			selectWindow(CompleteFileName);
			roiManager("Show All without labels");
			roiManager("Set Color", DyesHighlight[d]);
			roiManager("Set Line Width", 2);
			Count=roiManager("count");
			if(Count>0){
				roiManager("select all");
				roiManager("Combine");
				run("Add Selection...");
				selectWindow("Flat");
				roiManager("Show All without labels");
				roiManager("Set Color", DyesHighlight[d]);
				roiManager("Set Line Width", 2);
				roiManager("select All");
				roiManager("Combine");
				run("Flatten");
				close("Flat");
				selectWindow("Flat-1");
				rename("Flat");
			}
			roiManager("reset");
		}
	selectWindow("Flat");
	run("Select None");
	saveAs("tif", OutputDir+FileName+"-QC");
	close("*");
	}

	//Collecting the numbers
	WoundArea=newArray();
	CollagenArea=newArray();
	MuscleArea=newArray();
	FibrinArea=newArray();
	PCCollagenArea=newArray();
	PCMuscleArea=newArray();
	PCFibrinArea=newArray();
	
	for (A = 0; A < Area.length; A=A+3) {
		ScarArea=Area[A]+Area[A+1]+Area[A+2];
		WoundArea=Array.concat(WoundArea,ScarArea);
		
		Collagen=Area[A];
		PCCollagen=(Collagen/ScarArea)*100;
		CollagenArea=Array.concat(CollagenArea,Collagen);
		PCCollagenArea=Array.concat(PCCollagenArea,PCCollagen);
		
		Muscle=Area[A+1];
		PCMuscle=(Muscle/ScarArea)*100;
		MuscleArea=Array.concat(MuscleArea,Muscle);
		PCMuscleArea=Array.concat(PCMuscleArea,PCMuscle);
		
		Fibrin=Area[A+2];
		PCFibrin=(Fibrin/ScarArea)*100;
		FibrinArea=Array.concat(FibrinArea,Fibrin);
		PCFibrinArea=Array.concat(PCFibrinArea,PCFibrin);
	}

}
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);

//Pooling the results in a table and saving
Array.show("Results", FileNames, WoundArea,CollagenArea, MuscleArea,FibrinArea,PCCollagenArea,PCMuscleArea,PCFibrinArea);
saveAs("Results", OutputDir+year+"_"+month+1+"_"+dayOfMonth+"_ScarArea_Results.csv");
close("Results");

print("Macro done!");


