Fixed the automatic placement of ``WCSAxes`` tick labels so that a coordinate
whose ticks and tick labels are both hidden is no longer assigned a spine.
Previously such a coordinate could be assigned a spine because of its tick
count, pushing a visible coordinate onto a spine where it had no ticks. This
affected, for example, WCSes where one pixel axis maps to several world
coordinates.
