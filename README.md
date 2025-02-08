An executable file and the source code for the article: Tao H Q, Song G P, Jiang J, et al. Research on flight test scheduling under complex constraints[J]. Systems Engineering — Theory & Practice (under review).

The executable file and source code is intended solely for academic purposes. Unauthorized use for commercial or other non-academic purposes is strictly prohibited！

The executable file can automatically solve instances mentioned in this paper once a feasible instance path is provided.

If you are interested in the source code or our works, please privately contact taohanqiao19@nudt.edu.cn. We suggest readers to use Visual Studio 2019 (or higher) to compile the source code.

The code automatically generates test instances based on instances proposed by Otto et al. (2013). The function "generate_instance" defined in line 270 in 'heuristic.cpp' is responsible for generating instances.
Other parameters controlling the algorithm can be found in the header file "FTSP-P.h".
We also provide codes that use the CPLEX 20.10 to directly solve the problem.

If you find any problem or bug in the code, please feel free to contact taohanqiao19@nudt.edu.cn.

**References**

Otto, A., Otto, C., & Scholl, A. (2013). Systematic data generation and test design for solution on the example of SALBPGen for assembly line balancing. European Journal of Operation Research, 240, 32-42.
