//=================================
// include guard
#ifndef __GETEIGENINCLUDED__
#define __GETEIGENINCLUDED__


//--------------------
// Linux:
#ifdef __linux__ 
   #include <eigen/Eigen/Dense>
   #include  <eigen/Eigen/Core>
//--------------------
// Windows:
#elif _WIN32
   #include <Eigen\Dense>
   #include  <Eigen\Core>
//--------------------
// OSX (not correct yet)
#elif __APPLE__ 
   #include <eigen/Eigen/Dense>
   #include  <eigen/Eigen/Core>
#else
#endif



#endif