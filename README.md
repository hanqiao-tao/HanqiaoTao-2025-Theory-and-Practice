The source code for the article: Tao H Q, Song G P, Jiang J, et al. Research on flight test scheduling under complex constraints[J]. Systems Engineering — Theory & Practice (under review).
This code is intended solely for academic purposes. Unauthorized use for commercial or other non-academic purposes is strictly prohibited！
Here is a Visual Studio 2019 project provided.

To run this program, you need to:
1. Open the file "FTSP-P.h", and set the path of output file for the variable "output_file" defined in line 94;
2. Open the file "heuristic.cpp", and set the path of data file for the function "instance_reader_otto", see line 25;
3. Save the changes, and open the file "FTSP-P.sln", then run the program.

The code automatically generates test instances based on instances proposed by Otto et al. (2013). The function "generate_instance" defined in line 270 in 'heuristic.cpp' is responsible for generating instances.
Other parameters controlling the algorithm can be found in the header file "FTSP-P.h".
We also provide codes that use the CPLEX 20.10 to directly solve the problem.

If you discover any bugs in the code, please feel free to contact taohanqiao19@nudt.edu.cn.

References
Otto, A., Otto, C., & Scholl, A. (2013). Systematic data generation and test design for solution on the example of SALBPGen for assembly line balancing. European Journal of Operation Research, 240, 32-42.
