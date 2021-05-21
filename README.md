# nrg-seminar-project
Seminar project for the advanced computer graphics course at UL FRI

The JAR file of the java part of the project is located in the build folder. 
Run the program with the following command:
java -jar procedural-generation-of-volumetric-fluid-environemt.jar
The program will output a generateVolume.out file into its root folder.

RUn the VPT framework with node/packer and node/server-node command (opens on localhost:3000).
In VPT select generateVolume.out as the input file.
Set percision to 8bit, hasGradient to true and Number of samples to a value from 1-10.
The last value specifies the amount of models loaded from the file.

After the file loads set the following parameters for best results:
Extinction: 3
Scattering albedo: Nearly 1
Midtones: 0.25
Load the saved transfer function named BlueWaterBrownGround from the Build folder.

Finally you can use the iteration slider to iterate through the different loaded volumes.

