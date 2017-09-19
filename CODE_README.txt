Note: C++ code was written using Microsoft Visual Studio 2015.
      Code requires Eigen library which is already included. If need to include again:
      1. Open project, right-click LemkeHowson in Solution Explorer. Go to Properties.
      2. Configuration Properties -> C/C++ -> General -> Additional Include Directories -> <Edit...>
      3. Add eigen folder to path.


Instructions for using code:

Run the code (Ctrl+F5) to display the output - by default this is the two matrices, the Nash equilibrium (x1,x2) and the no. of pivots.
Adjust initial label on line 25.
Game consists of two random matrices A and B by default. Lines 41 to 43 have options for B = I and B = A^T. Comment/uncomment where appropriate.
To initialise matrices with specific elements, use comma initialisation (commented example found on lines 51 to 58). Must specify dimensions on line 21 and uncomment line 32. Must also comment lines 40-43 and line 30.
To find the worst case (go through all labels), uncomment the for loop on line 66 and its corresponding closing bracket on line 198. Also uncomment lines 150 and 200. Comment out lines 25 and 151.