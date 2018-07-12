# sims\_altSched
Scheduler code optimized for speed that observes blocks near the meridian. 

# Installation
## Required Dependencies
1. sims\_ocs
2. sims\_speedObservatory

Note that sims\_ocs requrires salpy to be installed. If you do not wish to
install salpy, a workaround is to comment out the salpy import statement
in the sims\_ocs code that throws an error when you try to run sims\_altSched.

## Optional Dependencies
1. pygame: pygame is not in the LSST stack and is only needed for visualizing
the scheduler as it runs. The visualization is very useful, but you can still
generate schedules using this code without pygame installed. You can install
pygame even on a read-only stack installation by running
```
pip install --user pygame
```
after loading the LSST environment.

# Running
First, read through the configuration flags/variables located at the top of
the `Simulator.py` file. These control the output of the simulator, and you
may want to change before running for the first time.

To run the simulator, run
```
python -m Simulator
```
in the terminal.

# Motivation Behind This Repository
This section describes why I wrote a new scheduler and a new simulator instead
of using opsim's. If you're just interested in how the code works and don't care
about my petty soap box speech, you can safely skip this section.

## Simulator

I wrote a new simulator to simulate my scheduler for a few reasons:
1. The opsim simulator and scheduler are integrated together well enough that
plugging in a new scheduler to the opsim simulator is difficult. Even with the
newly refactored opsim v4 (as of early September 2017), it would still be
necessary to pretend that the scheduler is actually a Proposal (which is not
an abstraction for an entire Scheduler). And mixing abstractions usually ends
poorly in my experience.
2. The interface that any scheduler running under opsim's simulator must
conform to is complicated to install, platform-dependent, and not suitable
for the request/reply pattern of inter-component communication demanded by a
simulator/scheduler combination.
3. Writing a new simulator is quite simple if the scheduler interface is
simple. The implementation in this repo is under 500 lines, but the vast
majority of that code controls the visualization and simulator output. 

## Scheduler

I wrote a new scheduler because opsim's scheduler is based around a type of
algorithm called a "greedy algorithm" that chooses whatever field appears to
be most meritorious at every time step based on an intuitive merit function.
I believe it is impossible to achieve optimal schedules with such an algorithm.

For those who don't already know about "greedy algorithms": a greedy algorithm
is an algorithm that optimizes a merit function at each time step. One important
limitation of these algorithms is that they are only able to optimize 
the merit function globally in a very limited set of cases.
Scheduling a telescope is clearly not one of these cases, since it is possible
that, by following the path of greatest local merit, the telescope is led away
from a region of very high merit that could have been reachable had the
scheduler previously taken a few lower-merit visits. Once the telescope has
slewed away from this region, the slew to reach this region is prohibitively
large. Using lookahead can only partly address this problem, as lookahead
(even with pruning) is exponentially computationally expensive in the amount
of lookahead time.

The problem is that we are giving the scheduler information about scheduling
constraints in a form that is hard for it to parse (i.e. as a merit
function instead of as an algorithm that satisfies the scheduling constraints
directly). We already know that the telescope should scan the meridian
as much as possible to minimize airmass. We know that it must achieve revisits
after ~30 minutes, and so to minimize slew times, we must arrange each 30
minute observing block to start and end in the same place (in alt/az space).
Combining just these two constraints yields a scheduler that scans up and down
on the meridian, executing N up scans and N down scans every ~30 minutes.
Without a rolling cadence, N must be 1 to avoid excessive intra-scan slews.

Opsim chose to encode these considerations into a merit function that it
attempts to optimize greedily. However, modulo filter changes, there are not
very many schedules that actually satisfy these requirements. In other words,
there aren't that many ways to scan up and down the meridian, and there aren't
any reasons obvious to me why the telescope should leave the meridian with any
regularity.

I implemented a scheduler that chooses one such meridian scanning strategy.
It was hastily constructed and has hardly been optimized, and yet
it largely reproduces or exceeds opsim's performance on the static survey.
For example, it achieves a 10% higher single visit depth (expressed in
effective observing seconds) than opsim v3.

This "marching-meridian" scheduler, however, was not intended to outperform
opsim on the static survey except for by the improvement in
seeing and skybrightness that result from the average normalized airmass
being reduced to 1.03. And in fact, I hear that opsim v4 stays closer to
the meridian, keeping the normalized airmass low as well.
Instead, this scheduler excels in time domain metrics, where,
from a theoretical standpoint, we would expect a greedy scheduling algorithm
to have the most trouble. According to the SNLots metric, which measures the
number of well-sampled transients that last on the order of a month, this
scheduler gets 5 times as many well sampled transients as opsim v3. It achieves
nearly every asteroid linkage revisit, while minion\_1020 misses ~20% (though
I haven't measured the baseline minion\_1016, which may do better). The
median single-band intranight gaps go down by a factor of several.

(These comparisons are all between my code and the opsim v3 results databases
that I have access to. I hear that opsim v4 has seen improvements, but there
are no opsim v4 results databases available to me as of this writing.)

The takeaway message is that there are alternatives to using a greedy algorithm
to choose every exposure. Although it's possible to make the static survey
relatively optimized simply by staying near the meridian (which I believe opsim
v4 does), it is much harder to optimize time domain metrics. Even knowing what
time domain metric performance is achievable is difficult, and if nothing else,
this scheduler demonstrates that a significant (5x) improvement in ~1-month
transient sampling is possible while maintaining static survey performance.

Feel free to email me at berkeley.edu, username drothchild, for a more
complete evaluation of this scheduler's performance.

# Overview of Public Interfaces

## Simulator

The `Simulator` class has one important public method: `run()`. This method
takes a `lsst.sims.speedObservatory.Telescope` as an argument, and it simulates
a survey with parameters sepcified in `config.py`. The output of the `Simulator`
is controlled by configuration variables at the top of `Simulator.py`.

Rough pseudocode for the `Simulator` is as follows:
```
initialize scheduler
for nightNum in range(numberOfNightsToSimulate):
    time = start time of nightNum'th night
    for Visit in scheduler.scheduleNight(nightNum):
        time += time to slew to and execute Visit 
        save Visit to disk
        if the night is over:
            break
```

The `Simulator`'s output is in the form of a csv with columns `mjd` (the
MJD of the beginning of the exposure), `ra`, `dec`, `prop` (which is either
wfd or dd, just for bookkeeping), and `filter`.

## Scheduler

The `Scheduler`, located at `minis/Scheduler.py`, has four public methods,
described in the sections below. The most important thing to note is that
the `Scheduler` preplans every night and has no ability to dynamically
reschedule if conditions change during the night or if the telescope falls
behind (or gets ahead of) the precomputed schedule. Rudimentary support for
this should probably be added, but on cloudy nights, it would probably be
better to switch over to a merit-based scheduler entirely.

### `scheduleNight(nightNum)`

This generator yields the visits that the scheduler decides to schedule on
the night that is passed in as the argument. It yields `None` if the scheduler
runs out of visits to schedule for the night. (In operations, this of course
should never happen, but I did not have time to dynamically generate new
fields when the preplanned ones ran out.)

The `Scheduler` keeps a reference to a `context` object that it can query
for the current time. `scheduleNight(...)` should therefore be used as a
generator, with the `context` object (i.e. the simulator) updating the simulated
time after every yielded visit. Calling `list(scheduleNight(nightNum))` will not
work, since the scheduler will have an incorrect time reference (and because
`notifyVisitComplete(...)` must also be called between each yielded visit)

### `notifyVisitComplete(visit, time)`

This method must be called every time a visit that `schedulNight(...)` yielded
was actually executed by the telescope (simulated or real).

### `notifyNightEnd()`

This method must be called every time the night is over. It allows the
scheduler to learn which fields it had planned to visit were never executed.

### `notifyDomeClosed(durationClosed)`

This method must be called every time the dome is closed so that the scheduler
knows to skip visits it had planned to execute during the closure time.

# Overview of the Scheduling Innards

I won't describe every public method of every class. Descriptions of those
can be found in comments within the code. However, I'll briefly mention the
structure of the code to orient newcomers to this repo.

## `NightScheduler`

The main work of planning out a night is accomplished in the `NightScheduler`
class, located under `minis/`. The main method is `schedule()`, which is a
generator that yields the visits for the night. There are also notification
methods analogous to those of `Scheduler`.

## Choosing Filters

The `NightScheduler` calls upon `minis/filtersequence.py` in order to assign
filters to visits. There are currently two strategies of filter scheduling
supported by `filtersequence.py`: one that changes filters between every
visit and revisit, and another that changes filters only twice a night at
the twilight boundaries. There are two methods for each of these
filter-changing strategies, one that returns the filters that should be
used, and the other to notify the `filtersequence` that a night has been
scheduled.

## Visits

There are two visit abstractions: a `Visit` and a `VisitPair`. The `VisitPair`
is just a collection of two `Visit`s, and each `Visit` specifies the parameters
for the exposure.

## Rotations

There is currently not support for varying the telescope rotator angle.

## Dithering

This scheduler does not use fixed fields. Instead, it uses a fixed tiling,
speficied in `minis/MiniSurvey.py`, that is rotated randomly to generate
a set of pointings that differ every night. Thus, positional dithering
is not required, and the resulting skymaps are much more homogeneous on
the 1 degree angular scale than those that can be obtained using fixed
fields and dithering.

## Deep Drilling

The scheduler has rudimentary support for deep drilling fields. If enabled
in `config.py`, the scheduler will schedule some hour-long blocks at the
field centers listed in `ddFields.csv`. These blocks will be Visits in the
`u` filter (for no particular reason) marked in the csv output as `dd` visits.
Right now (September 2017), the total time spent on these visits is about 4\%
of the survey duration.

Note that there is no control to make sure that these exposures don't hit the
zenith avoidance zone before they're completed. I think the code tries to avoid
that but doesn't enforce it. I added this support mostly just to show that
taking an hour out of the WFD observing sequence doesn't adversely affect the
WFD cadence by more than the reduction in total exposure time.

## Peripheral Code

As part of my development process, I wrote code that visualizes the scheduler
and code that keeps track of how many times each pixel is visited. This code
can be found in the `graphics/` directory, in `SkyMap.py`, and in
`SummaryPlots.py`. The graphics code relies on the `pygame` module, which is
not in the LSST stack. None of this code is necesassry for generating schedules.
