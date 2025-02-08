An executable file and source code are provided here for the article: Tao H Q, Song G P, Jiang J, et al. Research on flight test scheduling under complex constraints[J]. Systems Engineering — Theory & Practice (under review).

These materials are intended solely for academic purposes. Any unauthorized use for commercial or other non-academic purposes is strictly prohibited！

The executable file can automatically solve instances mentioned in this paper when a valid instance path is provided.

If you are interested in the source code or our work, please feel free contact us at taohanqiao19@nudt.edu.cn. We suggest readers to use Visual Studio 2019 (or newer revisions) for compiling the source code.

The code automatically generates test instances based on instances proposed by Otto et al. (2013). The function "generate_instance" defined in line 270 in 'heuristic.cpp' is responsible for generating instances.
Other parameters controlling the algorithm can be found in the header file "FTSP-P.h".
We also provide codes that use the CPLEX 20.10 to directly solve the problem.

If you encounter any issues or bugs in the code, please contact us at taohanqiao19@nudt.edu.cn.

**References**

Otto, A., Otto, C., & Scholl, A. (2013). Systematic data generation and test design for solution on the example of SALBPGen for assembly line balancing. European Journal of Operation Research, 240, 32-42.

**Additional Information**
1. If you use an Intel CPU with both P-cores and E-cores, we strongly recommend running the program on the same type of cores while using CPLEX. Due to the different speeds of P-cores and E-cores, CPLEX may report some errors when solving instances in parallel.
