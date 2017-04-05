#include <dmsim.h>
#include <sstream>
#include <iostream>

/****************************************/
/****************************************/

template <typename T>
T from_string(const char* str) {
   T val;
   std::istringstream iss(str);
   iss >> val;
   return val;
}

/****************************************/
/****************************************/

int main(int argc, char** argv) {
   /*
    * Parameters parsing.
    */
   if(argc < 10) {
      std::cerr << "[FATAL] Expecting 10 parameters, got "
                << argc << std::endl;
      abort();
   }
   uint32_t n  = from_string<uint32_t>(argv[1]); // number of robots
   double k    = from_string<double>  (argv[2]); // constant
   double r1   = from_string<double>  (argv[3]); // switch rate for 1 neighbor
   double r2   = from_string<double>  (argv[4]); // switch rate for 2 neighbors
   uint32_t ce = from_string<uint32_t>(argv[5]); // choose every
   double d    = from_string<double>  (argv[6]); // density
   double m    = from_string<double>  (argv[7]); // max speed
   int s       = from_string<uint32_t>(argv[8]); // random seed
   int g       = from_string<uint32_t>(argv[9]); // gnuplot
   /*
    * Simulation loop
    */
   sim_t sim(n, k, r1, r2, ce, d, m, s, g);
   sim.exec();
}

/****************************************/
/****************************************/
