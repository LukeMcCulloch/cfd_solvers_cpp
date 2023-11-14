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
   //#include <Eigen\Dense> 
   //#include  <Eigen\Core>
   // //quick hack to get Eigen into my particular build (you will need to link to Eigen your own way)
   #include "../../Eigen/Dense"
   #include "../../Eigen/Core"
#elif __CYGWIN__
   #include "../../Eigen/Dense"
   #include "../../Eigen/Core"
#elif __MINGW32__
   #include "../../Eigen/Dense"
   #include "../../Eigen/Core"

//--------------------
// OSX (not correct yet)
#elif __APPLE__ 
   #include <eigen/Eigen/Dense>
   #include  <eigen/Eigen/Core>
#else
#endif



#endif