# DMSim #

DMSim is a simulator for decision-making with movable agents.

## Compilation ##

To compile DMSim, execute these commands:

    $ mkdir build
    $ cd build
    $ cmake -DCMAKE_BUILD_TYPE=Release ../src
    $ make
    
This will create the executable file `build/dmsim`.

## Running DMSim ##

During execution, DMSim creates several data files. To have them all
neatly in one place, create a directory `data` and run DMSim for
there:

    $ mkdir data
    $ cd data

Then, launch DMSim with the following syntax:

    $ ../build/dmsim <n> <k> <r1> <r2> <ce> <d> <m> <s> <g>

All the options are mandatory and have the following meaning:

| Option | Meaning                             |
|--------|-------------------------------------|
| `n`    | number of agents                    |
| `k`    | `k` = `n`*`epsilon`                 |
| `r1`   | switch rate for 1 neighbor          |
| `r2`   | switch rate for 2 neighbors         |
| `ce`   | choose every `ce` simulation steps  |
| `d`    | robot density                       |
| `m`    | maximum agent speed                 |
| `s`    | random seed                         |
| `g`    | gnuplot output, `0` = off, `1` = on |

As said above, when you run dmsim you'll get several data files:

| Extension | Data type                                                           |
|-----------|---------------------------------------------------------------------|
| `com`     | Communication among robots for every decision step                  |
| `dec`     | Decision counts at every decision step                              |
| `exp`     | Position, speed, and decision of each robot at each simulation step |
| `nbr`     | How many robots have 0, 1, or >= 2 neighbors at every decision step |
| `grd`     | Used to debug the grid-based algorithm that optimizes communication |
| `gpl`     | GNUPlot script that generates GIF animations (created if `g`=`1`)   |

To generate the GIF animation, run GNUPlot as follows:

    $ gnuplot <file>.gpl

Be patient - GNUPlot will take quite some time to generate the animation!
