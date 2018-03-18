# How to Cite our Work
	
D. Griebler, J. Loff, G. Mencagli, M. Danelutto and L. G. Fernandes. **Efficient NAS Benchmark Kernels with C++ Parallel Programming**. *In proceedings of the 26th Euromicro International Conference on Parallel, Distributed and Network-Based Processing (PDP)*. Cambridge, United Kingdom, 2018.

# The NPB-CPP Benchmark

These codes were converted to **C++** from the original [NPB3.3.1](https://www.nas.nasa.gov/publications/npb.html). We achieved similar performance in **C++** compared to the **Fortran** version.

	==================================================================
		NAS Parallel Benchmarks in C++, OpenMP, FastFlow, and TBB
	 												
			Code contributors: 
					Dalvan Griebler    		
					Júnior Löff
													
		Warning: in case of problems send an email to us:					
			dalvan.griebler@acad.pucrs.br			
			junior.loff@acad.pucrs.br				
	==================================================================


This folder contains:

	- NPB-FF - Directory with the parallel version implemented in FastFlow
	- NPB-OMP - Directory with the parallel version translated from the original NPB version
	- NPB-SER - Directory with the serial version of the NPB ported to C++
	- NPB-TBB - Directory with the parallel version implemented in Thread Building Blocks

Each directory is independent and contains its own implemented version of the kernels:

	IS - Integer Sort, random memory access
	EP - Embarrassingly Parallel
	CG - Conjugate Gradient, irregular memory access and communication
	MG - Multi-Grid on a sequence of meshes, long- and short-distance communication, memory intensive
	FT - discrete 3D fast Fourier Transform, all-to-all communication

# Software Requiriments

*Warning: our tests were made with GCC-5*

**TBB**

*Installation*

	apt-get install libtbb-dev

**FastFlow** 

*Installation*

	svn co https://svn.code.sf.net/p/mc-fastflow/code/ $HOME/fastflow


# How to Compile 

Enter the directory from the version desired and execute:

	make _BENCHMARK CLASS=_VERSION


_BENCHMARKs are: 
		
	EP, CG, MG, IS and FT 
																										
_VERSIONs are: 
	
	Class S: small for quick test purposes
	Class W: workstation size (a 90's workstation; now likely too small)	
	Classes A, B, C: standard test problems; ~4X size increase going from one class to the next	
	Classes D, E, F: large test problems; ~16X size increase from each of the previous Classes  


Command:

	make ep CLASS=B