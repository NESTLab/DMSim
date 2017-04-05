#include "dmsim.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cstdio>

size_t DEBUG_INDENTATION = 0;

#define DEBUG(MSG, ...) { fprintf(stderr, "[DEBUG] "); for(size_t DEBUG_I = 0; DEBUG_I < DEBUG_INDENTATION; ++DEBUG_I) fprintf(stderr, "  "); fprintf(stderr, MSG, ##__VA_ARGS__); }

#define DEBUG_FUNCTION_ENTER { ++DEBUG_INDENTATION; DEBUG("%s - START\n", __PRETTY_FUNCTION__ ); }

#define DEBUG_FUNCTION_EXIT { DEBUG("%s - END\n", __PRETTY_FUNCTION__ ); --DEBUG_INDENTATION; }

#define TRACE(LINE) LINE; DEBUG(#LINE "\n");

/****************************************/
/****************************************/

static const double   PI            = 3.14159265358979323846;
static const uint32_t MAX_DECISIONS = 200;

/****************************************/
/****************************************/

/* RNG period parameters */
#define MT_N          624
#define MT_M          397
#define MT_MATRIX_A   0x9908b0dfUL /* constant vector a */
#define MT_UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define MT_LOWER_MASK 0x7fffffffUL /* least significant r bits */
#define MT_INT_MAX    0xFFFFFFFFUL

/*
 * A robot and its data.
 */
struct robot_t {
   uint32_t id;             // robot id
   double px, py;           // position
   double vx, vy;           // velocity
   bool cd, nd;             // current decision, next decision
   std::vector<int> gi, gj; // grid cells occupied
   std::vector<int> com;    // robots communicated with
};

std::ostream& operator<<(std::ostream& o, const robot_t& r) {
   o << r.px << "\t"
     << r.py << "\t"
     << r.vx << "\t"
     << r.vy << "\t"
     << (r.cd ? 2 : 3);
   return o;
}

/****************************************/
/****************************************/

grid_t::grid_t(uint32_t n,
               double w) :
   _s(log2(n)+1),
   _cs(_s / w),
   _c(new cell_t[_s * _s]) {}

/****************************************/
/****************************************/

grid_t::~grid_t() {
   delete[] _c;
}

/****************************************/
/****************************************/

cell_t& grid_t::operator()(uint32_t i,
                           uint32_t j) {
   if(i >= _s || j >= _s) {
      std::cerr << "[FATAL] Cell index ("
                << i
                << ","
                << j
                << ") out of bounds ("
                << _s
                << ","
                << _s
                << ")"
                << std::endl;
      abort();
   }
   return _c[j * _s + i];
}

/****************************************/
/****************************************/

#define W_TO_I(X) ((int32_t)::floor((X) * _cs))

cell_t& grid_t::operator[](robot_t& r) {
   return (*this)(W_TO_I(r.px), W_TO_I(r.py));
}

/****************************************/
/****************************************/

#define UPDATE_CELL(A,B) update(t,r,(A),(B));

void grid_t::place(uint64_t t,
                   robot_t& r) {
   /* Erase grid data */
   r.gi.clear();
   r.gj.clear();
   /* Take cells within communication perimeter */
   int32_t mini = std::max<int32_t>(0,    W_TO_I(r.px-1.));
   int32_t maxi = std::min<int32_t>(_s-1, W_TO_I(r.px+1.));
   int32_t minj = std::max<int32_t>(0,    W_TO_I(r.py-1.));
   int32_t maxj = std::min<int32_t>(_s-1, W_TO_I(r.py+1.));
   /* Update cell and its 8-radius */
   for(int32_t i = mini; i <= maxi; ++i)
      for(int32_t j = minj; j <= maxj; ++j)
         UPDATE_CELL(i, j);
}

/****************************************/
/****************************************/

void grid_t::update(uint64_t t,
                    robot_t& r,
                    uint32_t i,
                    uint32_t j) {
   /* Get cell */
   cell_t& c = _c[j * _s + i];
   /* Clear cell if its contents are out of date */
   if(c.t < t) {
      c.t = t;
      c.r.clear();
   }
   /* Add robot to cell */
   if(std::find(c.r.begin(), c.r.end(), &r) == c.r.end())
      c.r.push_back(&r);
   /* Update robot grid info */
   r.gi.push_back(i);
   r.gj.push_back(j);
}

/****************************************/
/****************************************/

sim_t::sim_t(uint32_t n,
             double k,
             double r1,
             double r2,
             uint32_t ce,
             double d,
             double m,
             int s,
             bool g) :
   _g(g),
   _e(k / n),
   _r1(r1), _r2(r2),
   _ce(ce),
   _d(d),
   _m(m),
   _t(0),
   _dc(0),
   _ct(n/2), _cf(n/2),
   _ws(::sqrt(n * PI / d)),
   _cmg(n, _ws) {
   /* Initialize random seed */
   _rngstate = new int32_t[MT_N];
   mt_setseed(s);
   /* Scatter robots */
   for(uint32_t i = 0; i < n; ++i) {
      robot_t* r = new robot_t;
      r->id = i;
      r->px = _ws * uniform();
      r->py = _ws * uniform();
      r->vx = _m * _ws * (2.0 * uniform() - 1.0);
      r->vy = _m * _ws * (2.0 * uniform() - 1.0);
      r->cd = (i < n/2);
      r->nd = r->cd;
      _r.push_back(r);
      _cmg.place(0, *r);
   }
   /*
    * Open data files
    */
   std::ostringstream bfn;
   bfn << n << "_"
       << _r1 << "_"
       << _r2 << "_"
       << _e << "_"
       << d << "_"
       << _m << "_"
       << _ce << "_"
       << s;
   _bfn = bfn.str();
   std::string fn;
   /* Open decision data file */
   fn = _bfn + ".dec";
   _df.open(fn.c_str(),
            std::ofstream::out | std::ofstream::trunc);
   if(_g) {
      /* Open experiment data file */
      fn = _bfn + ".exp";
      _ef.open(fn.c_str(),
               std::ofstream::out | std::ofstream::trunc);
      /* Open grid data file */
      fn = _bfn + ".grd";
      _gf.open(fn.c_str(),
               std::ofstream::out | std::ofstream::trunc);
      /* Open communication data file */
      fn = _bfn + ".com";
      _mf.open(fn.c_str(),
               std::ofstream::out | std::ofstream::trunc);
      /* Open neighbor data file */
      fn = _bfn + ".nbr";
      _nf.open(fn.c_str(),
               std::ofstream::out | std::ofstream::trunc);
      /* Open gnuplot script file */
      fn = _bfn + ".gpl";
      _sf.open(fn.c_str(),
               std::ofstream::out | std::ofstream::trunc);
      /* Output first data dump */
      gnuplot();
   }
}
   
/****************************************/
/****************************************/

sim_t::~sim_t() {
   if(_g) {
      /* Show last robot state */
      ++_t;
      gnuplot();
      /* Output script file */
      _sf << "ROBOTS=" << _r.size() << std::endl
          << "STEPS=" << _t << std::endl
          << "EVERY=" << _ce << std::endl
          << "DF='" << _bfn << ".dec'" << std::endl
          << "EF='" << _bfn << ".exp'" << std::endl
          << "GF='" << _bfn << ".grd'" << std::endl
          << "CF='" << _bfn << ".com'" << std::endl
          << "NF='" << _bfn << ".nbr'" << std::endl
          << "OF='" << _bfn << ".gif'" << std::endl
          << "WS=" << _ws << std::endl
          << "CS=" << _cmg.cells_per_side() << std::endl
          << "set terminal gif size 800,400 animate loop 0 delay 10" << std::endl
          << "set output OF" << std::endl
          << "set linetype 1 lw 1 pt 7 ps 1" << std::endl
          << "set linetype 2 lc rgb 'red'" << std::endl
          << "set linetype 3 lc rgb 'blue'" << std::endl
          << "set linetype 4 lc rgb 'green' lw 1 pt 4 ps 1" << std::endl
          << "set linetype 5 lc rgb 'black' lw 1" << std::endl
          << "do for [II=0:STEPS] {" << std::endl
          << "   set multiplot title sprintf('N=%d, e=%f, r=(%f,%f), choose every=%d, d=%f, S=%f, step %d',ROBOTS," << _e << "," << _r1 << "," << _r2 << "," << _ce << "," << _d << "," << _m << ",II)" << std::endl
          << "   set margins 3,0,2,2.5" << std::endl
          << "   unset grid" << std::endl
          << "   unset xlabel" << std::endl
          << "   unset ylabel" << std::endl
          << "   unset key" << std::endl
          << "   unset ytics" << std::endl
          << "   set origin 0,0" << std::endl
          << "   set size 0.5,1" << std::endl
          << "   set xrange [0:WS]" << std::endl
          << "   set yrange [0:WS]" << std::endl
          << "   set xtics (0,WS) format '%h'" << std::endl
          << "   plot EF every :::II::II using 1:2:5 with points lt 1 lc var notitle, \\" << std::endl
          << "        EF every :::II::II with vectors lt 1 lc var notitle" << std::endl
          << "   if(II % EVERY == 0) {" << std::endl
          << "      set origin 0,0" << std::endl
          << "      set size 0.5,1" << std::endl
          << "      unset xtics" << std::endl
          << "      plot CF every :::II::II with vectors lt 5 notitle" << std::endl
          << "   }" << std::endl
          << "   set origin 0,0" << std::endl
          << "   set size 0.5,1" << std::endl
          << "   set xrange [0:CS]" << std::endl
          << "   set yrange [0:CS]" << std::endl
          << "   set xtics 1 format ''" << std::endl
          << "   set ytics 1 format ''" << std::endl
          << "   set grid xtics ytics back linestyle 0" << std::endl
          << "   plot GF every :::(II*ROBOTS)::(II*ROBOTS+ROBOTS) using ($1+0.5):($2+0.5) with points lt 4 notitle" << std::endl
          << "   set origin 0.5,0.5" << std::endl
          << "   set size 0.5,0.5" << std::endl
          << "   set margins 7,3,2,2.5" << std::endl
          << "   set xrange [0:(STEPS / EVERY)]" << std::endl
          << "   set yrange [0:ROBOTS]" << std::endl
          << "   set key default" << std::endl
          << "   set xtics auto format '%h'" << std::endl
          << "   set ytics auto format '%h'" << std::endl
          << "   set xlabel 'Decision #' offset 0,0.5" << std::endl
          << "   set ylabel 'Robots' offset 2,0" << std::endl
          << "   plot DF every ::::(floor(II / EVERY)) using 1:2 with lines lt 2 title 'A', \\" << std::endl
          << "        DF every ::::(floor(II / EVERY)) using 1:3 with lines lt 3 title 'B'" << std::endl
          << "   set origin 0.5,0" << std::endl
          << "   set size 0.5,0.5" << std::endl
          << "   set margins 7,3,3,1.5" << std::endl
          << "   plot NF every ::::(floor(II / EVERY)) using 1:2 with lines lt 2 title 'no neigh', \\" << std::endl
          << "        NF every ::::(floor(II / EVERY)) using 1:3 with lines lt 3 title 'one neigh', \\" << std::endl
          << "        NF every ::::(floor(II / EVERY)) using 1:4 with lines lt 4 title 'two neigh'" << std::endl
          << "   unset multiplot" << std::endl
          << "}" << std::endl;
   }
   /* Dispose of RNG state */
   delete[] _rngstate;
   /* Close data files */
   if(_df.is_open()) _df.close();
   if(_g) {
      if(_ef.is_open()) _ef.close();
      if(_gf.is_open()) _gf.close();
      if(_mf.is_open()) _mf.close();
      if(_nf.is_open()) _nf.close();
      if(_sf.is_open()) _sf.close();
   }
   /* Remove robots */
   while(!_r.empty()) {
      delete _r.back();
      _r.pop_back();
   }
}
   
/****************************************/
/****************************************/

int sim_t::exec() {
   while(!done()) step();
   return 0;
}

/****************************************/
/****************************************/

void sim_t::step() {
   /* Increase time step */
   ++_t;
   /* Sense and control only when necessary */
   if(_t % _ce == 0) {
      /* Decision counter */
      ++_dc;
      /* Update grid */
      for(int i = 0; i < _r.size(); ++i)
         _cmg.place(_t, *_r[i]);
      /* Neighbor counters */
      uint32_t nc0 = 0, nc1 = 0, nc2 = 0;
      /* Sense + control */
      for(int i = 0; i < _r.size(); ++i) {
         /* Erase list of robots in communication */
         _r[i]->com.clear();
         /* Pick new random speed */
         _r[i]->vx = _m * _ws * (2.0 * uniform() - 1.0);
         _r[i]->vy = _m * _ws * (2.0 * uniform() - 1.0);
         /* Get cell */
         cell_t& c = _cmg[*_r[i]];
         /* Make a list of robots in range */
         std::list<robot_t*> rr;
         for(std::list<robot_t*>::iterator it = c.r.begin();
             it != c.r.end(); ++it) {
            if(((*it)->px - _r[i]->px) * ((*it)->px - _r[i]->px) +
               ((*it)->py - _r[i]->py) * ((*it)->py - _r[i]->py) < 1.0) {
               rr.push_back(*it);
            }
         }
         /* Get the number of neighbors */
         uint32_t nn = rr.size() - 1;
         if(nn == 0) {
            /* No neighbors, spontaneous switching */
            if(bernoulli(_e)) {
               _r[i]->nd = !_r[i]->nd;
            }
            if(_g) {
               /* Increase no-neighbor counter */
               ++nc0;
            }
         }
         else {
            /* Make a list of neighbors */
            std::list<robot_t*> nl = rr;
            nl.erase(std::find(nl.begin(), nl.end(), _r[i]));
            /* Make a decision depending on how many neighbors we have */
            if(nn == 1) {
               /* Single neighbor */
               if(_r[i]->cd != (*nl.begin())->cd && // robots disagree
                  bernoulli(_r1)) {                 // switch rate
                  _r[i]->nd = !_r[i]->cd;
               }
               if(_g) {
                  /* Add to communication list */
                  _r[i]->com.push_back((*nl.begin())->id);
                  /* Increase one-neighbor counter */
                  ++nc1;
               }
            }
            else {
               /* At least two neighbors */
               robot_t *n1, *n2;
               if(nn == 2) {
                  /* Exactly two neighbors */
                  n1 = *nl.begin();
                  n2 = *std::next(nl.begin(), 1);
               }
               else {
                  /* More than two neighbors, pick two at random */
                  /* Get neighbor 1 */
                  uint32_t ni = nn * uniform();
                  std::list<robot_t*>::iterator it = nl.begin();
                  for(uint32_t j = 0; j < ni; ++j) ++it;
                  n1 = *it;
                  nl.erase(it);
                  /* Get neighbor 2 */
                  ni = nl.size() * uniform();
                  it = nl.begin();
                  for(uint32_t j = 0; j < ni; ++j) ++it;
                  n2 = *it;
               }
               /* Make decision */
               if(n1->cd == n2->cd    && // neighbors agree
                  _r[i]->cd != n1->cd && // robot disagrees with neighbors
                  bernoulli(_r2)) {      // switch rate
                  _r[i]->nd = !_r[i]->cd;
               }
               if(_g) {
                  /* Add to communication list */
                  _r[i]->com.push_back(n1->id);
                  _r[i]->com.push_back(n2->id);
                  /* Increase two-neighbor counter */
                  ++nc2;
               }
            }
         }
         /* Update choice counts */
         if(_r[i]->nd != _r[i]->cd) {
            /* Robot switched decision */
            if(_r[i]->nd) {
               /* Robot switched from false to true */
               --_cf;
               ++_ct;
            }
            else {
               /* Robot switched from true to false */
               ++_cf;
               --_ct;
            }
         }
      }
      /* Write decision data line */
      _df << (_t / _ce) << "\t"
          << _ct << "\t"
          << _cf << std::endl;
      /* Write neighbor data line */
      _nf << (_t / _ce) << "\t"
          << nc0 << "\t"
          << nc1 << "\t"
          << nc2 << std::endl;
   }
   /* Show robots state */
   if(_g) gnuplot();
   /* Act */
   double opx, opy;
   for(int i = 0; i < _r.size(); ++i) {
      /* Update decision */
      _r[i]->cd = _r[i]->nd;
      /* Save old position */
      opx = _r[i]->px;
      opy = _r[i]->py;
      /* Update position */
      _r[i]->px += _r[i]->vx;
      _r[i]->py += _r[i]->vy;
      /* Manage collisions with walls */
      if(_r[i]->px < 0. ||
         _r[i]->px >= _ws) {
         _r[i]->vx = -_r[i]->vx;
         _r[i]->px = opx + _r[i]->vx;
      }
      if(_r[i]->py < 0. ||
         _r[i]->py >= _ws) {
         _r[i]->vy = -_r[i]->vy;
         _r[i]->py = opy + _r[i]->vy;
      }
      if(_r[i]->px < 0. ||
         _r[i]->px >= _ws) {
         std::cerr << "x="
                   << _r[i]->px
                   << " is out of bounds (0, "
                   << _ws
                   << ")"
                   << std::endl;
         abort();
      }
      if(_r[i]->py < 0. ||
         _r[i]->py >= _ws) {
         std::cerr << "y="
                   << _r[i]->py
                   << " is out of bounds (0, "
                   << _ws
                   << ")"
                   << std::endl;
         abort();
      }
   }
}

/****************************************/
/****************************************/

bool sim_t::done() {
   return (_dc == MAX_DECISIONS) || (_ct == 0) || (_cf == 0);
}

/****************************************/
/****************************************/

void sim_t::mt_setseed(uint32_t seed) {
   _rngstate[0] = seed & 0xffffffffUL;
   for(_rngidx = 1; _rngidx < MT_N; ++_rngidx) {
      _rngstate[_rngidx] =
         (1812433253UL * (_rngstate[_rngidx-1] ^ (_rngstate[_rngidx-1] >> 30)) + _rngidx);
      _rngstate[_rngidx] &= 0xffffffffUL;
   }
}

/****************************************/
/****************************************/

uint32_t sim_t::mt_uniform32() {
   uint32_t y;
   static uint32_t mag01[2] = { 0x0UL, MT_MATRIX_A };
   /* mag01[x] = x * MATRIX_A  for x=0,1 */
   if (_rngidx >= MT_N) { /* generate N words at one time */
      int32_t kk;
      for (kk = 0; kk < MT_N - MT_M; ++kk) {
         y = (_rngstate[kk] & MT_UPPER_MASK) | (_rngstate[kk+1] & MT_LOWER_MASK);
         _rngstate[kk] = _rngstate[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      for (; kk < MT_N - 1; ++kk) {
         y = (_rngstate[kk] & MT_UPPER_MASK) | (_rngstate[kk+1] & MT_LOWER_MASK);
         _rngstate[kk] = _rngstate[kk+(MT_M-MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      y = (_rngstate[MT_N-1] & MT_UPPER_MASK) | (_rngstate[0] & MT_LOWER_MASK);
      _rngstate[MT_N-1] = _rngstate[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
      _rngidx = 0;
   }
   y = _rngstate[_rngidx++];
   /* Tempering */
   y ^= (y >> 11);
   y ^= (y << 7) & 0x9d2c5680UL;
   y ^= (y << 15) & 0xefc60000UL;
   y ^= (y >> 18);
   return y;
}

/****************************************/
/****************************************/

double sim_t::uniform() {
   return (double)mt_uniform32() / (double)MT_INT_MAX;
}

/****************************************/
/****************************************/

bool sim_t::bernoulli(double p) {
   return uniform() < p;
}

/****************************************/
/****************************************/

void sim_t::gnuplot() {
   /* Experiment data block */
   _ef << "# Time step " << _t << std::endl;
   /* Communication data block */
   _mf << "# Time step " << _t << std::endl;
   /* Create spare communication data block for when no communication occurs */
   _mf << "0 0 0 0" << std::endl;      
   /* Go through all the robots */
   for(uint32_t i = 0; i < _r.size(); ++i) {
      /* Experiment data record */
      _ef << *_r[i] << std::endl;
      /* Grid data block */
      _gf << "# Time step " << _t << ", robot " << i << std::endl;
      for(uint32_t j = 0; j < _r[i]->gi.size(); ++j) {
         /* Grid data record */
         _gf << _r[i]->gi[j]
             << "\t"
             << _r[i]->gj[j]
             << std::endl;
      }
      /* End of grid data block */
      _gf << std::endl;
      /* Communication data records */
      for(uint32_t j = 0; j < _r[i]->com.size(); ++j) {
         /* Communication data record */
         _mf << _r[i]->px << "\t"
             << _r[i]->py << "\t"
             << (_r[_r[i]->com[j]]->px - _r[i]->px) << "\t"
             << (_r[_r[i]->com[j]]->py - _r[i]->py)
             << std::endl;
      }
   }
   /* End of experiment data block */
   _ef << std::endl;
   /* End of communication experiment data block */
   _mf << std::endl;
}

/****************************************/
/****************************************/
