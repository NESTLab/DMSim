#ifndef DMSIM_H
#define DMSIM_H

#include <cstdint>
#include <fstream>
#include <list>
#include <vector>

struct robot_t;

/****************************************/
/****************************************/

/*
 * A grid cell
 */
struct cell_t {
   uint64_t t;            // time of last update
   std::list<robot_t*> r; // list of robots in this cell
   /* Constructor */
   cell_t() : t(0) {}
};

/*
 * A grid. Used for sensing and communication.
 */
class grid_t {

public:

   /**
    * Constructor.
    */
   explicit grid_t(uint32_t s, // cells on a side
                   double w);  // world side size

   /**
    * Destructor.
    */
   ~grid_t();

   /**
    * Returns cell at given position.
    */
   cell_t& operator()(uint32_t i,
                      uint32_t j);

   /**
    * Returns cell at given position.
    */
   cell_t& operator[](robot_t& r);

   /** Places a robot in the grid */
   void place(uint64_t t,  // current time
              robot_t& r); // robot

   /** Returns the number of cells per side */
   uint32_t cells_per_side() const {
      return _s;
   }

private:

   /** Updates a specific cell */
   void update(uint64_t t,  // current time
               robot_t& r,  // robot
               uint32_t i,  // grid position i
               uint32_t j); // grid position j

private:

   uint32_t _s; // number of cells on a grid side
   double _cs;  // cell side length
   cell_t*  _c; // the cells
};

/****************************************/
/****************************************/

/*
 * A grid of robots.
 */
class sim_t {

public:

   /**
    * Constructor.
    */
   sim_t(uint32_t n,  // number of robots
         double k,    // k = n * epsilon
         double r1,   // rec rate 1
         double r2,   // rec rate 2
         uint32_t ce, // choose every
         double d,    // robot density
         double m,    // max speed
         int s,       // random seed
         bool g       // gnuplot
      );

   /**
    * Destructor.
    */
   ~sim_t();

   /**
    * Executes the simulation.
    */
   int exec();

private:

   /**
    * Simulation step.
    */
   void step();

   /**
    * Returns true when the simulation is done.
    */
   bool done();

   /**
    * Sets the seed of a Mersenne-Twister random number generator.
    */
   void mt_setseed(uint32_t seed);

   /**
    * Returns a random value between 0 and 0xFFFFFFFFUL.
    */
   uint32_t mt_uniform32();

   /**
    * Returns a random value between 0. and 1.
    */
   double uniform();

   /**
    * Returns the result of a coin toss with probability p.
    */
   bool bernoulli(double p);

   void gnuplot();

private:

   /** GNUPlot output */
   bool _g;
   /** Epsilon */
   double _e;
   /** Rec rates */
   double _r1, _r2;
   /** Choose every */
   uint32_t _ce;
   /** Robot density */
   double _d;
   /** Max speed */
   double _m;
   /** Current step */
   uint64_t _t;
   /** Decision counter */
   uint32_t _dc;
   /** Choice counters (true, false) */
   uint32_t _ct, _cf;
   /** The world side length */
   double _ws;
   /** The communication grid of the robots */
   grid_t _cmg;
   /** The list of robots */
   std::vector<robot_t*> _r;
   /* Base file name */
   std::string _bfn;
   /** The decision data file */
   std::ofstream _df;
   /** The experiment data file */
   std::ofstream _ef;
   /** The grid data file */
   std::ofstream _gf;
   /** The communication data file */
   std::ofstream _mf;
   /** The neighbor data file */
   std::ofstream _nf;
   /** The gnuplot script file */
   std::ofstream _sf;
   /** Random number generator state */
   int32_t* _rngstate;
   /** Random number generator index */
   uint32_t _rngidx;
};

#endif
