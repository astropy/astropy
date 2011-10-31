/* SKY CALCULATOR PROGRAM 
   John Thorstensen, Dartmouth College.  
   This program computes many quantities frequently needed by the 
   observational astronomer.  It is written as a completely 
   self-contained program in standard c, so it should be
   very transportable; the only issue I know of that really affects
   portability is the adequacy of the double-precision floating
   point accuracy on the machine.  Experience shows that c compilers
   on various systems have idiosyncracies, though, so be sure
   to check carefully.

   This is intended as an observatory utility program; I assume the
   user is familiar with astronomical coordinates and nomenclature.
   While the code should be very transportable, I also 
   assume it will be installed by a conscientious person who
   will run critical tests before it is released at a new site.
   Experience shows that some c compilers generate unforseen errors
   when the code is ported, so the output should be checked meticulously
   against data from other sites.
 

   After prompting the user for site information, the program enters
   its main command loop,   

	- enter RA, dec, proper motion, epoch, date, time, 
	     new site parameters, etc. ... (when the code starts,
	     it initializes these things to the current date and time
             and the celestial coordinates of the zenith at the
	     specified site).

	- print an almanac of rise/set and twilight for the sun and
             moon for the specified night

	- compute and display circumstances of observation for the
	   current parameters, including precessed coordinates,
	   airmass, interference from moon or twilight, parallactic
	   angle, etc; the program also gives calendar date in 
	   both UT and local, Julian date, and barycentric corrections. 

	- compute and display a table of airmasses (etc) at
	   hourly intervals through the night.  This is very useful
	   at the telescope.  Also, if one has a modest number of 
	   objects, it may be convenient (using system utilities)
	   to redirect the output and print a hard copy of these 
	   tables for ready reference.

	- compute and display galactic and ecliptic coordinates.
   
	- compute and display rough (of order 0.1 degree, but often
	  much better) positions of the major planets.

    The program is self-contained.  It has been developed primarily on 
   UNIX and Linux machines, and should adapt easily to any system with 
   a c compiler.  
   
	** BUT CAUTION ... **
   Because many of the routines take a double-precision floating point
   Julian Date as their time argument, one must be sure that the machine
   and compiler carry sufficient mantissa to reach the desired accuracy.
   On most architectures the double-precision floating point julian date
   has an accuracy of order 0.01 seconds of time, which is just adequate.

LEGALITIES: 

   I make no guarantee as to the accuracy, reliability, or
   appropriateness of this program, though I have found it to be 
   reasonably accurate and quite useful to the working astronomer.
   
   The program is COPYRIGHT 2000 BY JOHN THORSTENSEN.  
   Permission is hereby granted for non-profit scientific or educational use.
   For-profit use (e. g., by astrologers!) must be through negotiated
   license.  The author requests that observatories and astronomy 
   departments which install this as a utility notify the author
   by paper mail, just so I know how widely it is used.     

   Credits:  
    * The julian date and sidereal time routines were 
    originally coded in PL/I by  Steve Maker of Dartmouth College.  
    They were based on routines in the old American Ephemeris.
    Many of the routines were coded from Jean Meeus' "Astronomical
    Formulae for Calculators", published by Willman-Bell.  This is
    an extraordinarily helpful little book!

    BUGS:

    The program can enter an infinite loop of prompting in response
    to incorrect input.  This is fairly rare in practice, but
    may be bothersome on multi-user systems.  If anyone can document
    this properly, I may be able to fix it.  In the meantime you
    may want to run this only on systems on which you can easily 
    stop runaway processes (e.g., by CRTL-C on Linux or UNIX machines.)

    APOLOGIES/DISCLAIMER:
    I am aware that the code here does not always conform to
    the best programming practices.  Not every possible error condition
    is anticipated, and no guarantee is given that this is bug-free.
    Nonetheless, most of this code has been shaken down at several
    hundred sites for several years, and I have never received any
    actual bug reports.  Many users have found this program
    to be useful.

    CHANGES SINCE THE ORIGINAL DISTRIBUTION ....

	The program as listed here is for the most part similar to that
	posted on the IRAF bulletin board in 1990.  Some changes
	include:

	01 In the original code, many functions returned structures, which
	   some c implementations do not like.  These have been eliminated.

	02 The original main() was extremely cumbersome; much of it has
	   been broken into smaller (but still large) functions.

	03 The hourly airmass includes a column for the altitude of the
	   sun, which is printed if it is greater than -18 degrees.

	04 The planets are included (see above).  As part of this, the
	   circumstances calculator issues a warning when one is within
	   three degrees of a major planet.  This warning is now also
	   included in the hourly-airmass table.
    
	05 The changeover from standard to daylight time has been rationalized.
	   Input times between 2 and 3 AM on the night when DST starts (which
	   are skipped over and  hence don't exist) are now disallowed; input 
	   times between 1 and 2 AM on the night when DST ends (which are
	   ambiguous) are interpreted as standard times.  Warnings are printed
	   in both the almanac and calculator mode when one is near to the
	   changeover.

	06 a much more accurate moon calculation has been added; it's used
	   when the moon's coordinates are given explicitly, but not for
	   the rise/set times, which iterate and for which a lower precision 
	   is adequate.

	07 It's possible now to set the observatory elevation; in a second 
	   revision there are now two separate elevation parameters specified.
	   The elevation above the horizon used only in rise/set calculations 
	   and adjusts rise/set times assuming the parameter is the elevation 
	   above flat surroundings (e. g., an ocean).  The true elevation above
	   sea level is used (together with an ellipsoidal earth figure) in
	   determining the observatory's geocentric coordinates for use in 
	   the topocentric correction of the moon's position and in the 
	   calculation of the diurnal rotation part of the barycentric velocity
	   correction.  These refinements are quite small.

	08 The moon's altitude above the horizon is now printed in the 
	   hourly airmass calculation; in the header line, its illuminated
	   fraction and angular separation from the object are included,
	   as computed for local midnight.

	09 The helio/barycentric corrections have been revised and improved.
	   The same routines used for planetary positions are used to
	   compute the offset from heliocentric to solar-system
	   barycentric positions and velocities.  The earth's position
	   (and the sun's position as well) have also been improved somewhat.

	10 The printed day and date are always based on the same truncation
	   of the julian date argument, so they should now always agree
	   arbitrarily close to midnight.

	11 A new convention has been adopted by which the default is that the
	   date specified is the evening date for the whole night.  This way,
	   calculating an almanac for the night of July 3/4 and then specifying
	   a time after midnight gives the circumstances for the *morning of
	   July 4*.  Typing 'n' toggles between this interpretation and a
	   literal interpretation of the date.

	12 The planetary proximity warning is now included in the hourly airmass
	   table.

	13 A routine has been added which calculates how far the moon is from
	   the nearest cardinal phase (to coin a phrase) and prints a 
	   description.  This information is now included in both the almanac
	   and the calculator mode.

	14 The output formats have been changed slightly; it's hoped this 
	   will enhance comprehensibility.

	15 A substantial revision affecting the user interface took place 
	   in September of 1993.  A command 'a' has been added to the 
	   'calculator' menu, which simply prints the almanac (rise, set,
	   and so on) for the current night.  I'd always found that it was
	   easy to get disoriented using the '=' command -- too much 
	   information about the moment, not enough about the time
	   context.  Making the almanac info *conveniently* available
	   in the calculator mode helps your get oriented.

	   When the 'a' almanac is printed, space is saved over the 
	   almanac printed on entry, because there does not need
	   to be a banner introducing the calculator mode.  Therefore some 
	   extra information is included with the 'a' almanac; this includes
	   the length of the night from sunset to sunrise, the number of 
	   hours the sun is below -18 degrees altitude, and the number of hours 
	   moon is down after twilight.  In addition, moonrise and moonset
	   are printed in the order in which they occur, and the occasional
	   non-convergence of the rise/set algorithms at high latitude are
	   signalled more forcefully to the user.

	16 I found this 'a' command to be convenient in practice, and never
	   liked the previous structure of having to 'quit' the calculator
	   mode to see almanac information for a different night.  D'Anne
	   Thompson of NOAO also pointed out how hokey this was, especially the
	   use of a negative date to exit. So, I simply removed the outer 
	   'almanac' loop and added a 'Q' to the main menu for 'quit'.  The
	   use of upper case -- for this one command only --  should guard
	   against accidental exit.

	17 The menu has been revised to be a little more readable.  

	18 More error checking was added in Nov. 1993, especially for numbers.  
	   If the user gives an illegal entry (such as a number which isn't
	   legal), the rest of the command line is thrown away (to avoid
	   having scanf simply chomp through it) and the user is prompted
	   as to what to do next.  This seems to have stopped all situations
	   in which execution could run away.  Also, typing repeated carriage
	   returns with nothing on the line -- which a frustrated novice
	   user may do because of the lack of any prompts -- causes a
	   little notice to be printed to draw attention to the help and menu
	   texts.
	 
	19 I found in practice that, although the input parameters and
	   conditions are deducible *in principle* from such things as the
	   'a' and '=' output, it takes too much digging to find them.  So
	   I instituted an 'l' command to 'look' at the current parameter
	   values.  To make room for this I put the 'Cautions and legalities'
	   into the 'w' (inner workings) help text.  This looks as though
	   it will be be very helpful to the user. 

	20 The more accurate moon calculation is used for moonrise and
	   moonset; the execution time penalty appears to be trivial. 
	   Low precision moon is still used for the summary moon information
	   printed along with the hourly airmass table.

	21 A calculation of the expected portion of the night-sky 
	   brightness due to moonlight has been added.  This is based on
	   Krisciunas and Schaefer's analytic fits (PASP, 1991).  Obviously,
	   it's only an estimate which will vary considerably depending on
	   atmospheric conditions.

	22 A very crude estimator of the zenith sky brightness in twilight
	   has been added.

	23 A topocentric correction has been added for the sun, in anticipation
	   of adding eclipse prediction.

	24 The code now checks for eclipses of the sun and moon, by making
	   very direct use of the predicted positions.  If an eclipse is
	   predicted, a notice is printed in print_circumstances; also, a
	   disclaimer is printed for the lunar sky brightness if a lunar
	   eclipse is predicted to be under way.

	25 In the driver of the main calculator loop, a provision has been
	   added for getting characters out of a buffer rather than reading
	   them directly off the command line.  This allows one to type any
	   valid command character (except Q for quit) directly after a number 
	   in an argument without generating a complaint from the program 
	   (see note 18).  This had been an annoying flaw.

	26 In 1993 December/1994 January, the code was transplanted
	   to a PC and compiled under Borland Turbo C++, with strict
	   ANSI rules.  The code was cut into 9 parts -- 8 subroutine
	   files, the main program, and an incude file containing
	   global variables and function prototypes.

	27 An "extra goodies" feature has been added -- at present it
	   computes geocentric times of a repeating phenomenon as
	   viewed from a site.  This can be used for relatively little-used
           commands to save space on the main menu.  

	28 The sun and moon are now included in the "major planets"
	   printout.  This allows one to find their celestial positions
	   even when they are not visible from the selected site.
	
	29 A MAJOR new feature was added in February 1994, which computes
           the observability of an object at new and full moon over a 
           range of dates.  The galactic/ecliptic coordinate converter
           was moved to the extra goodies menu to make room for this.

	30 Inclusion of a season-long timescale means that it's not
           always necessary to specify a date on entry to the program.
           Accordingly, the program immediately starts up in what used
           to be called "calculator" mode -- only the site is prompted
           for.  It is thought that the site will be relevant to nearly
           all users.

	31 Because the user is not led by the hand as much as before, the
           startup messages were completely revised to direct new users
           toward a short `guided tour' designed to show the program's 
	   command structure and capabilities very quickly.  Tests on 
	   volunteers showed that users instinctively go for anything 
	   called the `menu', despite the fact that that's a slow way to 
	   learn, so all mention of the menu option is removed from the 
	   startup sequence; they'll catch on soon enough.

	32 Code has been added to automatically set the time and
           date to the present on startup.  A menu option 'T' has been
           added to set the time and date to the present plus a settable
           offset.  This should be very useful while observing.
	
	33 Because Sun machines apparently do not understand ANSI-standard
           function declarations, the code has been revised back to K&R
           style.  It's also been put back in a monolithic block for
           simplicity in distribution.
	
	34 The startup has been simplified still further, in that the
           coordinates are set automatically to the zenith on startup.
	   An 'xZ' menu item sets to the zenith for the currently specified
           time and date (not necessarily the real time and date.)

	35 Another MAJOR new capability was added in early 1994 --
           the ability to read in a list of objects and set the current
	   coordinates to an object on the list.  The list can be sorted
           in a number of ways using information about the site, date
           and time.
	
	35 Calculator-like commands were added to the extra goodies menu
           to do precessions and to convert Julian dates to calendar
           dates.  An option to set the date and time to correspond to
           a given julian date was also added.
	
	36 Another substantial new capability was added Aug 94 -- one can
           toggle open a log file (always named "skyclg") and keep
           a record of the output.  This is done simply by replacing
           most occurrences of "printf" with "oprintf", which mimics
           printf but writes to a log file as well if it is open.
	   This appears to slow down execution somewhat.

	37 12-degree twilight has been added to the almanac.  While the
	   awkward "goto" statements have been retained, the statement
           labels have been revised to make them a little clearer.	   

	38 The precession routine was generalized to include nutation and
	   aberration, and routines to calculate these effects are now
           included.  Nearly all calls to the new routine leave the
           aberration and nutation out, but 'xa' extra goodies gives
           apparent place up to aberration and nutation.  Tests against the
	   FK5 and "Apparent Places of Fundamental Stars" shows the
           agreement to 0.1 arcsec or smaller.  

	39 Atmospheric refraction is computed and reported in the "xa"
           apparent place calculation.  Barometric pressure is guesstimated
           from the observatory elevation above sea level.

	40 An airmass based on a series expansion is printed in the "=" 
           command, in place of the secant of the zenith distance, provided
           sec z is < 12.   The approximation breaks down larger than this.

	41 The parallactic angle calculation has been simplified and folded
           inside the routine which calculates altitude and azimuth. 

	42 A few utility routines (atan_circ, xyz_cel) were cleaned up.

	43 I added a new option for reading the system clock every time
	   'timely' output is asked for, and computing for the updated 
	   time & date.  This is toggled on and off with 'xU', for 
	   updating.  The date and time can optionally be offset from that 
	   read from the system clock, as in the 'T' option.

	44 The sexigesimal coordinate-output routine "put_coords" was
           revised to be more general and a little less ugly.  It's not
           necessarily shorter!  This required gussying up the outputs
           because the new routine doesn't space just the same. 
           It's now possible to include a "+" sign explicitly in the
           output when desired.

        45 The calendrical conversion routines have been replaced with
           routines coded from Meeus' "Astronomical Formulae for 
           Calculators", to avoid the use of proprietary Numerical
           Recipes code.

        46 The routine used to sort objects by airmass, etc. has been
           replaced with an original coding of the heapsort algorithm
           to avoid the use of proprietary "Numerical Recipes" code.
           All the code *should* now be freely distributable.

        47 The routine which reads strings from input has been rewritten
           so that it can take delimited strings like 18:02:12.33 ...
           delimiter can be any character which is not a numeral,
           a + or - sign, or a decimal point (i.e., anything which can't
           be part of a number.)

	48 Another 'xtra goodies' option has been added to print out the
	   parallax factors (i.e., parallax displacement in ra and dec
           a star would have at exactly 1 pc distance at the nominal
           time and date).  Also prints XYZ for earth and annual aberration
           (displacement of star due to finite speed of light and earth's
           motion about the sun.)

        49 In 2000 October I fixed some infelicities with the 
           xv and xf commands, and updated the xd (TDT - UT) command.
           The TeX documentation and the on-line documentation were 
	   updated.  In 2002 Jan. finally fixed a bug which displayed the 
	   wrong ephemeris in the header in the xv command 
	   (computations were OK).
	
	50 Added the capability to read and write site-parameter files.
	
	51 replaced "round" with "roundx" to avoid name conflict arising
	   with a math library function.

	52 2003 -- updated the etcorr routine again to reflect Delta T 
	   values in 2003 almanac and adjust forward extrapolation to 
	   match last tabulated value.

        53 Added code adapted from F. Ochsenbein (CDS Strasbourg) to give the
	   constellation location of the point in question, and report the
	   location in the "=" and "xc" commands.
	
	54 Added a "parallactic penalty" factor, i.e. the magnitude of 
	   tan z times the magnitude of the sine of the parallactic angle.
	   This severity of the parallactic differential refraction effect
	   will be proportional to this, for a slit which defaults to 
	   north-south. 
	
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>
#include <string.h>

/* a couple of the system-dependent magic numbers are defined here */

#define SYS_CLOCK_OK 1    /* 1 means ANSI-standard time libraries do work,
   2 means they don't.  This is used by compiler switches in file 5 and
   the main program.  */

#define LOG_FILES_OK 1  /* 1 means that log files are enabled.
			Any other value means they're not.  */

#define MAX_OBJECTS 2000 
#define MINSHORT -32767   /* min, max short integers and double precision */
#define MAXSHORT 32767
#define MAXDOUBLE 1.0e38
#define MINDOUBLE -1.0e38
#define BUFSIZE 150 

#define XFORM_FROMSTD  1  /* defined quantities for apparent place transforms .. */
#define XFORM_TOSTDEP  -1
#define XFORM_JUSTPRE  1
#define XFORM_DOAPPAR  0

/* some (not all) physical, mathematical, and astronomical constants
   used are defined here. */

#define  PI                3.14159265358979
#define  TWOPI             6.28318530717959
#define  PI_OVER_2         1.57079632679490  /* From Abramowitz & Stegun */   
#define  ARCSEC_IN_RADIAN  206264.8062471
#define  DEG_IN_RADIAN     57.2957795130823
#define  HRS_IN_RADIAN     3.819718634205
#define  KMS_AUDAY         1731.45683633   /* km per sec in 1 AU/day */
#define  SPEED_OF_LIGHT    299792.458      /* in km per sec ... exact. */
#define  SS_MASS           1.00134198      /* solar system mass in solar units */
#define  J2000             2451545.        /* Julian date at standard epoch */
#define  SEC_IN_DAY        86400.
#define  FLATTEN           0.003352813   /* flattening of earth, 1/298.257 */
#define  EQUAT_RAD         6378137.    /* equatorial radius of earth, meters */
#define  ASTRO_UNIT        1.4959787066e11 /* 1 AU in meters */
#define  RSUN              6.96000e8  /* IAU 1976 recom. solar radius, meters */
#define  RMOON             1.738e6    /* IAU 1976 recom. lunar radius, meters */
#define  PLANET_TOL        3.          /* flag if nearer than 3 degrees 
						to a major planet ... */
#define  KZEN              0.172       /* zenith extinction, mag, for use
				     in lunar sky brightness calculations. */
#define FIRSTJD            2415387.  /* 1901 Jan 1 -- calendrical limit */
#define LASTJD             2488070.  /* 2099 Dec 31 */

/* MAGIC NUMBERS which might depend on how accurately double-
   precision floating point is handled on your machine ... */

#define JDRESOLUTION 3.5e-8   /* 3 milliseconds */

#define  EARTH_DIFF        0.05            /* used in numerical
   differentiation to find earth velocity -- this value gives
   about 8 digits of numerical accuracy on the VAX, but is 
   about 3 orders of magnitude larger than the value where roundoff
   errors become apparent. */

#define  MIDN_TOL          0.00001         /* this is no longer
   used -- it was formerly
   how close (in days) a julian date has to be to midnight
   before a warning flag is printed for the reader.  VAX
   double precision renders a Julian date considerably
   more accurately than this.  The day and date are now based
   on the same rounding of the julian date, so they should
   always agree. */

/*  FUNCTION PROTOTYPES and type definitions ....
    These are used in breaking the code into function libraries.
    They work properly on a strictly ANSI compiler, so they
    apparently comply with the ANSI standard format.  */

/* these global variables determine whether the program updates the
time and date whenever asked for output. */

#if SYS_CLOCK_OK == 1 
int update_on = 0; 
double update_delta = 0.;  
#endif

struct coord
   {
     short sign;  /* carry sign explicitly since -0 not neg. */
     double hh;
     double mm;
     double ss;
   };

struct date_time
   {
	short y;
	short mo;
	short d;
	short h;
	short mn;
	float s;
   };

FILE *sclogfl = NULL;

double star_tzero, star_terr, 
	star_period, star_perr;  /* for ephemeris calculations ... global */

void oprntf(char *fmt, ...)

/* This routine should look almost exactly like printf in terms of its
   arguments (format list, then a variable number of arguments
   to be formatted and printed).  It is designed to behave just
   like printf (though perhaps not all format types are supported yet)
   EXCEPT that IF the globally-defined file pointer "sclogfl" is
   defined, IT ALSO WRITES TO THAT FILE using fprintf.  The skeleton
   for this came from Kernighan and Ritchie, 2nd edition, page 156 --
   their "minprintf" example.  I modified it to include the 
   entire format string (e.g., %8.2f, %7d) and to write to the 
   file as well as standard output.  */

{
	va_list ap;        /* see K&R for explanation of these macros */
	char *p, *sval;
	char outform[10];  /* an item's output format, e.g. %8.2f */
	char strout[150]; 
	int ival, i;
	short shval;
	char cval;
	double dval;

	va_start(ap,fmt);
	for (p = fmt; *p; p++) {
		if (*p != '%') {
			putchar(*p);
/* overkill to put in these preprocessor flags, perhaps. */
#if LOG_FILES_OK == 1 
			if(sclogfl != NULL) fputc(*p,sclogfl);
#endif
			continue;
		}
		i = 0;
		outform[i] = '%';
		p++;
		while(*p != 'd' && *p != 'f' && *p != 's' && *p != 'c'
		   && *p != 'h') {
			outform[++i] = *p++;
		}				
		switch (*p) {
		case 'd':			
			ival = va_arg(ap, int);
			outform[++i] = *p;
			outform[++i] = '\0';
			printf(outform, ival);
#if LOG_FILES_OK == 1 
			if(sclogfl != NULL)
				fprintf(sclogfl,outform,ival);
#endif
			break;
		case 'h':    /* signals short argument ... */			
			shval = va_arg(ap, int);
			outform[++i] = 'd';
			outform[++i] = '\0';
			++p;  /* skip the 'd' in '%hd' */
			printf(outform, shval);
#if LOG_FILES_OK == 1 
			if(sclogfl != NULL)
				fprintf(sclogfl,outform,shval);
#endif
			break;
		case 'c':			
			/* cval = va_arg(ap, char); 
			  ... compiler problem on Sun machines */
			cval = va_arg(ap, int);
			outform[++i] = *p;
			outform[++i] = '\0';
			printf(outform, cval);
#if LOG_FILES_OK == 1 
			if(sclogfl != NULL)
				fprintf(sclogfl,outform,cval);
#endif
			break;
		case 'f':
			dval = va_arg(ap, double);
			outform[++i] = *p;
			outform[++i] = '\0';
			printf(outform, dval);
#if LOG_FILES_OK == 1 
			if(sclogfl != NULL) 
				fprintf(sclogfl,outform,dval);
#endif
			break;
		case 's':
			outform[++i] = *p;
			outform[++i] = '\0';
			i = 0;
			for (sval = va_arg(ap, char *); *sval; sval++) {
				strout[i++] = *sval;
			}
			strout[i] = '\0';
			printf(outform,strout);
#if LOG_FILES_OK == 1 
			if(sclogfl != NULL) fprintf(sclogfl,outform,strout);
#endif
			break;
		default:
	                ;
		}
	}
	va_end(ap);
}

/* elements of K&R hp calculator, basis of commands */

char buf[BUFSIZE];
int bufp=0;


char getch() /* get a (possibly pushed back) character */
{
	return((bufp > 0) ? buf[--bufp] : getchar());
}

void ungetch(c) /* push character back on input */
	int c;
{
	if(bufp > BUFSIZE)
		printf("Ungetch -- too many characters.\n");
	else 
		buf[bufp++] = c;
}

/* some functions for getting well-tested input. */

int legal_num_part(c)
	char c; 

{
	if((c != '.') && (c != '-') && (c != '+') && (c < '0' || c > '9')) 
		return(-1);  /* not a legal number part */
	else return(0);
}
 
int legal_int_part(c)
	char c; 
   
{
	if((c != '-') && (c < '0' || c > '9')) 
		return(-1);  /* not a legal number part */
	else return(0);
} 

int legal_command_char(c)
	char c;
{
	/* Allows more sophisticated argument checking by seeing if
	      a character appended to an argument is actually a 
	      legal commmand. */
	switch(c) {
	  case '?': return(1);
		break;
	  case 'i': return(1);
		break;
	  case 'f': return(1);
		break;
	  case 'w': return(1);
		break;  
	  case 'r': return(1);
		break;
	  case 'd': return(1);
		break;
	  case 'C': return(1);
		break;
	  case 'y': return(1);
		break;
	  case 't': return(1);
		break;
	  case 'T': return(1);
		break;
	  case 'n': return(1);
		break;
	  case 'g': return(1);
		break;
	  case 'e': return(1);
		break;
	  case 'p': return(1);
		break;
	  case 's': return(1);
		break;
	  case 'l': return(1);
		break;
	  case '=': return(1);
		break;
	  case 'a': return(1);
		break;
	  case 'h': return(1);
		break;
	  case 'o': return(1);
		break;
	  case 'm': return(1);
		break;
	  case 'c': return(1);
		break;
	  case 'x': return(1);
		break;
	  /* let's not allow 'Q' here! */
	  default: return(0);
	}       
}

int parsedouble(s,d)
               
	char *s;
	double *d;
                  
	/* return values 0 = ok, with number, 1 = found a valid command,
	   but no number, and -1 = an error of some sort (unexpected char)*/

{
   short i=0, legal = 0;

   while((*(s+i) != '\0') && (legal == 0)) {
	if(legal_num_part(*(s+i)) == 0) i++;       
	else if(legal_command_char(*(s+i)) == 1) {
		/* to allow command to follow argument without blanks */
		ungetch(s[i]);
		*(s+i) = '\0';  /* will terminate on next pass. */
	}
	else legal = -1;
   } 
   
   if(legal == 0) {
	if(i > 0) {
		sscanf(s,"%lf",d);
		return(0);
	} 
	else if (i == 0) { /* ran into a command character -- no input */
		*d = 0.;
		return(1);  /* ok, actually */
	}
   }
   else {
	printf("%s is not a legal number!! Try again!\n",s);
	return(-1);
   }
}


int getdouble(d,least,most,errprompt) 
             
	double *d,least,most;
	char *errprompt; 
             

{
    char s[30], buf[200], c;
    short success = -1, ndiscard = 0;

    scanf("%s",s);
    while(success < 0) {
	success = parsedouble(s,d);
	if((success == 0) && ((*d < least) || (*d > most))) {
	   printf("%g is out of range; allowed %g to %g -- \n",
			*d,least,most);
	   success = -1;
	} 
	if(success < 0) {
	   /* if there's error on input, clean out the rest of the line */
	   ndiscard = 0;
	   while((c = getchar()) != '\n')  {
		buf[ndiscard] = c;
		ndiscard++;
	   }
	   if(ndiscard > 0) {
		buf[ndiscard] = '\0';  /* terminate the string */         
		printf("Rest of input (%s) has been discarded.\n",buf);
	   }
	   printf("%s",errprompt);
	   printf("\nTry again:");
	   scanf("%s",s);
	}
    }
    return((int) success);
}    

int parseshort(s,d)
              
	char *s;
	short *d;
              
{
   short i=0, legal = 0;

   while((*(s+i) != '\0') && (legal == 0)) {
	if(legal_int_part(*(s+i)) == 0) i++;
	else if(legal_command_char(*(s+i)) == 1) {
		/* to allow command to follow argument without blanks */
		ungetch(s[i]);
		*(s+i) = '\0';  /* will terminate on next pass. */
	}
	else legal = -1;
   } 

   if(legal == 0) {
	if(i > 0) {
		sscanf(s,"%hd",d);
		return(0);
	} 
	else if (i == 0) { /* ran into a command character -- no input */
		*d = 0.;
		return(1);  /* didn't get a number, but something else legal */
	}
   } 
   else {
	printf("%s is not a legal integer number!! Try again!\n",s);
	return(-1);
   }
}

int getshort(d,least,most,errprompt) 

	short *d,least,most;
	char *errprompt;

{
    char s[30];
    short success = -1, ndiscard = 0;
    char c, buf[200];

    scanf("%s",s);
    while(success < 0) {
	success = parseshort(s,d);
	if((success == 0) && ((*d < least) || (*d > most))) {
	   printf("%d is out of range; allowed %d to %d -- try again.\n",
			*d,least,most);
	   success = -1;
	} 
	if(success < 0) {
	   /* if there's error on input, clean out the rest of the line */
	   ndiscard = 0;
	   while((c = getchar()) != '\n')  {
		buf[ndiscard] = c;
		ndiscard++;
	   }
	   if(ndiscard > 0) {
		buf[ndiscard] = '\0';  /* cap the string */       
		printf("Rest of input (%s) has been discarded.\n",buf);
	   }
	   printf("%s",errprompt);
	   printf("Try again:");
	   scanf("%s",s);
	}
    }
    return( (int) success);
}    


double bab_to_dec(bab)
                 
	struct coord bab;
                 
   /* converts a "babylonian" (sexigesimal) structure into 
      double-precision floating point ("decimal") number. */
   {
   double x;
   x = bab.sign * (bab.hh + bab.mm / 60. + bab.ss / 3600.);
   return(x);
   }

void dec_to_bab (deci,bab)
                
	double deci;
	struct coord *bab;
                
   /* function for converting decimal to babylonian hh mm ss.ss */

{
   int hr_int, min_int;

   if (deci >= 0.) bab->sign = 1; 
   else {
      bab->sign = -1;
      deci = -1. * deci;
   }
   hr_int = deci;   /* use conversion conventions to truncate */
   bab->hh = hr_int;
   min_int = 60. * (deci - bab->hh);
   bab->mm = min_int;
   bab->ss = 3600. * (deci - bab->hh - bab->mm / 60.);
}

short get_line(s)
              
	char *s;   	

/* gets a line terminated by end-of-line and returns number of characters. */
{       
	char c;
	short i = 0;

	c = getchar(); /* get the first character */
	/* chew through until you hit non white space */
	while((c == '\n') || (c == ' ') || (c == '\t')) c = getchar();

	s[i]=c;
	i++;

	/* keep going til the next newline */
	while((c=getchar()) != '\n') {
		s[i]=c;
		i++;
	}
	s[i]='\0';  /* terminate with null */
	return(i);
}

int is_delimited(char *instrng) 

/* Checks to see if the instrng has non-numeric
   characters, which are then assumed to be delimiters
   (like colons). */

{
	int val = 0;

	while(*instrng != '\0') {
		if((*instrng >= '0' && *instrng <= '9') || *instrng == '-'
			|| *instrng == '+') {
			instrng++;
		}
		else return(1);
	}
	return(0);
}


double conv_delimited_coord(char *instrng) 

/* Handed a string instr, converts to a double-precision
   coordinate under the assumption that the string is a
   continguous, delimiter-separated sexagesimal triplet, e.g.
   -0 degr 12m 02.33s would be -0:12:2.33 or -0:12:02.33.
   The whole thing has to be continguous (no whitespace!).
*/

{

   char substr[20], *place, *substrplace;
   int hrs, min, sign = 1;
   double secs, val; 
   
   place = instrng;
   substrplace = substr;
   while(*place == ' ' || *place == '\t') place++;
   if(*place == '+') {
	sign = 1;
	place++;
   }
   if(*place == '-') {
	sign = -1;
	place++;
   }

   while(*place >= '0' && *place <= '9') {
	*substrplace = *place;
	place++; substrplace++;
   }
   *substrplace = '\0';
   sscanf(substr,"%d",&hrs);
   substrplace=substr;
   place++;  /* skip the colon */

   while(*place >= '0' && *place <= '9') {
	*substrplace = *place;
	place++; substrplace++;
   }
   *substrplace = '\0';
   sscanf(substr,"%d",&min);
   substrplace=substr;
   place++;  /* skip the colon */

   while((*place >= '0' && *place <= '9') || *place == '.' ) {
	*substrplace = *place;
	place++; substrplace++;
   }
   *substrplace = '\0';
   sscanf(substr,"%lf",&secs);

   val = sign * (hrs + min / 60. + secs / 3600.);
   return(val);
}

double get_coord()

/* Reads a string from the terminal and converts it into
   a double-precision coordinate.  This is trickier than 
   it appeared at first, since a -00 tests as non-negative; 
   the sign has to be picked out and handled explicitly. */
/* Prompt for input in the calling routine.*/
{
   short sign;
   double hrs, mins, secs;
   char hh_string[6];  /* string with the first coord (hh) */
   char hh1;   
   char errprompt[80];
   short i = 0;
   int end_in = 0;

   /* read and handle the hour (or degree) part with sign */

   scanf("%s",hh_string);

   if(is_delimited(hh_string) == 1) {  /* It's colon-delimited, or something */
	return(conv_delimited_coord(hh_string));
   }
	
   hh1 = hh_string[i];

   while(hh1 == ' ') {
       /* discard leading blanks */
       i++;
       hh1 = hh_string[i];
   }

   if(hh1 == '-') sign = -1;   /* get sign explicitly */
   else if(hh1 == '+') sign = 1;
   else sign = 1;

   if((end_in = parsedouble(hh_string,&hrs)) < 0) {
	printf("Didn't parse correctly -- set parameter to zero!!\n");
	return(0.);
   }

   if(sign == -1) hrs = -1. * hrs;

   /* read in the minutes and seconds normally */           
   if(end_in == 0)      
       	end_in = getdouble(&mins,0.,60.,
	  "Give minutes again, then seconds; no further prompts.\n");
   else return(sign * hrs);
   if(end_in == 0) end_in = getdouble(&secs,0.,60.,
     "Give seconds again, no further prompts.\n");
   else if(end_in == 1) secs = 0.;
   return(sign * (hrs + mins / 60. + secs / 3600.));
}

double roundx(x,places) 

double x;
int places;

/* rounds argument x to places places, e.g. 2.32839,1 -> 2.3. */

{
	double tmp, base = 1.;
	int i, ip;	

	for(i = 1; i <= places; i++) {  /* bet this is faster than pow ... */
		base *= 10.;
	}
	tmp = x * base;
	if(tmp >= 0.) 
        	tmp += 0.5;
        else tmp -= 0.5;	
	ip = (int) tmp;
	tmp = ((double) ip) / base;
	return(tmp);
}

round_coord(incoord, outcoord, prec) 

struct coord *incoord, *outcoord;
int prec;  

/* Rounds the seconds of a struct coord to a specified precision;
   if they turn out to be sixty, does the carry to the other fields.

	precision 0 -- whole minutes  (seconds set to zero)
	          1 -- tenths of a minute (seconds set to zero)
		  2 -- whole seconds
		  3 -- tenths of a second
	  	  4 -- hundredths ... 
				etc.
*/

{
	outcoord->sign = incoord->sign;

/* initialize */
	
	outcoord->ss = incoord->ss,prec;
	outcoord->mm = incoord->mm;
	outcoord->hh = incoord->hh;

	if(prec <= 1) {
		outcoord->mm = roundx((outcoord->mm + outcoord->ss / 60.),prec);
		outcoord->ss = 0.;
		if(outcoord->mm >= 59.99) {  /* permissible because of limit 
				on prec */
			outcoord->mm -= 60.;
			outcoord->hh += 1.;
		}
	}
	else {
		outcoord->ss = roundx(outcoord->ss,(prec-2));
		if( outcoord->ss >= 59.999999999) {  /* as many digits as
				one would ever want ... */
			outcoord->ss -= 60.;
			outcoord->mm += 1.;
			if(outcoord->mm >= 59.999999999) {
				outcoord->mm -= 60.;
				outcoord->hh += 1.;
			}
		}
	}
}

void put_hrs(hrs,sign,width,showpos,alignsign) 

double hrs;   
short sign;
int width, showpos, alignsign;

/* Puts out the hours (or decimal degrees) with the 
   following format information:

	-- Allows "width" digits of space for hours;
		e.g. -20 would be width 2. 
	-- if showpos == 1, prints a + sign if result is
                positive.
        -- if alignsign == 1, prints sign before the field;
		otherwise places sign flush with digit.
*/

{
	int i, digitsout, leadblanks;
	char outform[20];
	double tmp;

	if(alignsign == 1) {
		if(sign < 0) oprntf("-");
		else if(showpos == 1) oprntf("+");
		sprintf(outform,"%%%d.0f",width);
		oprntf(outform,hrs);
	}
	else {
		tmp = fabs(hrs);
		digitsout = 1;
		while(tmp >= 10.) {
			digitsout++;
			tmp /= 10.;
		}
		if(digitsout >= width) {
			if(sign < 0) oprntf("-");
			else if (showpos == 1) oprntf("+");
			oprntf("%.0f",hrs);
		}
		else {
			for(i = 1; i < width - digitsout; i++) 
				oprntf(" ");
			if(sign < 0) oprntf("-");	
			else if (showpos == 1) oprntf("+");
			else oprntf(" ");
			sprintf(outform,"%%%d.0f",digitsout);
			oprntf(outform,hrs);
		}
	}
}	


void put_coords(deci, prec, showsign) 

double deci;  /* decimal version */
int prec, showsign;   

{
	struct coord bab, babout;
	char formstr[20];
	int outstringlen;

	dec_to_bab(deci,&bab);
	
	round_coord(&bab,&babout,prec);
	
	if(prec == 0) {
	        put_hrs(babout.hh, babout.sign, 3, showsign, 0);
		oprntf(" %02.0f",babout.mm);
	}
	else if (prec == 1) {
	        put_hrs(babout.hh, babout.sign, 3, showsign, 0);
		oprntf(" %04.1f",babout.mm);
	}	
	else {
		if(prec == 2) {
	        	put_hrs(babout.hh, babout.sign, 3, showsign, 0);
		    	oprntf(" %02.0f %02.0f",
				babout.mm,babout.ss);
		}
		else {
	            put_hrs(babout.hh, babout.sign, 3, showsign, 0);
		    sprintf(formstr," %%02.0f %%0%d.%df",
			prec+1,prec-2);
		    oprntf(formstr,babout.mm,babout.ss);	
		}
	}
}	

void chomp(char *s) {

/* remplaces first newline character in a string with end-of-string. */

	char *p;
	
	p = s;
	while(*p != '\n' && *p != '\0') p++;
	if(*p == '\n') *p = '\0';
	
}
	
void load_site(longit,lat,stdz,use_dst,
	zone_name,zabr,elevsea,elev,horiz,site_name)
              
	double *longit,*lat;
   	double *stdz;
	short *use_dst;
   	char *zone_name, *zabr;
	double *elevsea;
   	double *elev,*horiz;
	char *site_name;
               
/* sets the site-specific quantities; these are
		longit     = W longitude in decimal hours
		lat        = N latitude in decimal degrees
		stdz       = standard time zone offset, hours
		elevsea    = elevation above sea level (for absolute location)
		elev       = observatory elevation above horizon, meters
		horiz      = (derived) added zenith distance for rise/set due
				to elevation
		use_dst    = 0 don't use it
			     1 use USA convention
			     2 use Spanish convention
			     < 0 Southern hemisphere (reserved, unimplimented)
		zone_name  = name of time zone, e. g. Eastern
		zabr       = single-character abbreviation of time zone
		site_name  = name of site.  */

{
	short nch;
	char obs_code[3];  /* need only one char, but why not? */
	char errprompt[50], buf[100];
	FILE *sitef = NULL;
	char sitefname[200];
	int ok;

	printf("*SELECT SITE* - Enter single-character code:\n");
	printf("   n .. NEW SITE, prompts for all parameters.\n");
	printf("   W .. Write site parameters to a file.\n");
	printf("   R .. Read new site parameters from a file.\n");
	printf("   x .. exit without change (current: %s)\n",site_name);
	printf("   k .. Kitt Peak [MDM Obs.]\n");
	printf("   s .. Shattuck Obs., Dartmouth College, Hanover NH\n");
	printf("   z .. SALT/SAAO, Sutherland, SA.\n");
	printf("   e .. European Southern Obs, La Silla\n");
	printf("   v .. Very Large Telescope, Cerro Paranal\n");
	printf("   a .. Anglo-Australian Telelescope, Siding Spring\n");
	printf("   h .. Mt. Hopkins, AZ (MMT, FLWO)\n");
	printf("   p .. Palomar Observatory\n");
	printf("   c .. Las Campanas Observatory \n");
	printf("   t .. Cerro Tololo \n");
	printf("   r .. Roque de los Muchachos, La Palma, Canary Is.\n");
	printf("   d .. Dominion Astrophysical Obs., Victoria, BC\n");
	printf("   o .. McDonald Observatory, Mt. Locke, Texas\n");
	printf("   m .. Mauna Kea, Hawaii\n");
	printf("   l .. Lick Observatory\n");
	printf("Your answer --> ");
	scanf("%s",obs_code);
	if(obs_code[0] == 'x') {
		printf("No action taken.  Current site = %s.\n",site_name);
		return;
	}
	else if(obs_code[0] == 'W') {
		printf("Filename for saving site parameters:");
		scanf("%s",sitefname);
		if((sitef = fopen(sitefname,"w")) != NULL) {
			fprintf(sitef,"%s\n",site_name);
			fprintf(sitef,"%s\n",zone_name);
			fprintf(sitef,"%c %hd %lf %lf %lf %lf %lf",
					*zabr,*use_dst,*longit,*lat,*stdz,
					*elevsea,*elev);
			printf("Site param file %s written.\n",
				sitefname);
		}
		else printf("ERROR -- Couldn't open file %s\n",
			sitefname);
		fclose(sitef);
	}
	else if(obs_code[0] == 'R') {
		printf("Filename for site parameters:");
		scanf("%s",sitefname);
		sitef = fopen(sitefname,"r");
		ok = 1;
		if(sitef != NULL) {
			if(fgets(buf,100,sitef) != NULL) {
				/* careful -- allow blanks in site name */
				chomp(buf);
				strcpy(site_name,buf);
			}
			else ok = 0;
			if(fgets(buf,100,sitef) != NULL) {
				/* also allow blanks in zone name */
				chomp(buf);
				strcpy(zone_name,buf);
			}
			else ok = 0;
			if(fgets(buf,100,sitef) != NULL) {
				/* the rest is numeric or single-char */
				if(sscanf(buf,"%c %hd %lf %lf %lf %lf %lf",
					zabr,use_dst,longit,lat,stdz,elevsea,elev)
					!= 7) ok = 0; 
			}
			else ok = 0;
			fclose(sitef);
		}
		else {
			printf("FAILED: Couldn't open file %s ... \n",sitefname);
			ok = 0;
		}
		if(ok == 1) printf("Loaded parameters for %s\n",site_name);
		else printf("CAREFUL -- Something wrong with the file!\n");
	}
	else if(obs_code[0] == 'k') {
		strcpy(site_name,"Kitt Peak [MDM Obs.]");
		strcpy(zone_name, "Mountain");
		*zabr = 'M';
		*use_dst = 0;
		*longit = 7.44111; /* decimal hours */
		*lat = 31.9533;    /* decimal degrees */
		*stdz = 7.;
		*elevsea = 1925.;  /* for MDM observatory, strictly */ 
		*elev = 700.;  /* approximate -- to match KPNO tables */
	}
	else if (obs_code[0] == 's') {
		strcpy(site_name, "Shattuck Observatory");
		strcpy(zone_name,"Eastern");
		*zabr = 'E';
		*use_dst = 1;
                *longit = 4.81907;  /* from GPS */
                *lat = 43.705;
		*stdz = 5.;
		*elevsea = 183.;
		*elev = 0.;  /* below surrouding horizon */
	}	
	else if (obs_code[0] == 'z') {
		strcpy(site_name, "SALT/SAAO, Sutherland");
		strcpy(zone_name,"South African");
		*zabr = 'S';
		*use_dst = 0;
                *longit = -1.38744;  
                *lat = -32.3783;
		*stdz = -2;
		*elevsea = 1771.;
		*elev = 00.;  /* guess! */
	}

	else if (obs_code[0] == 'e') {
		strcpy(site_name, "ESO, Cerro La Silla");
		strcpy(zone_name, "Chilean");
		*zabr = 'C';
		*use_dst = -1;
		*longit = 4.7153;
		*lat = -29.257;
		*stdz = 4.;
		*elevsea = 2347.; 
		*elev = 2347.; /* for ocean horizon, not Andes! */
		printf("\n\n** Will use daylght time, Chilean date conventions. \n\n");
	}
	else if (obs_code[0] == 'v') {
		strcpy(site_name, "VLT, Cerro Paranal");
		strcpy(zone_name, "Chilean");
		*zabr = 'C';
		*use_dst = -1;
		*longit = 4.69356;
		*lat = -24.625;
		*stdz = 4.;
		*elevsea = 2635.; 
		*elev = 2635.; /* for ocean horizon, not Andes! */
		printf("\n\n** Will use daylght time, Chilean date conventions. \n\n");
	}
	else if (obs_code[0] == 'p') {
		strcpy(site_name, "Palomar Observatory");
		strcpy(zone_name, "Pacific");
		*zabr = 'P';
		*use_dst = 1;
		*longit = 7.79089;
		*lat = 33.35667;
		*elevsea = 1706.;
		*elev = 1706.;  /* not clear if it's appropriate ... */
		*stdz = 8.;
	}
	else if (obs_code[0] == 't') {
		strcpy(site_name, "Cerro Tololo");
		strcpy(zone_name, "Chilean");
		*zabr = 'C';
		*use_dst = -1;
		*longit = 4.721;
		*lat = -30.165;
		*stdz = 4.;
		*elevsea = 2215.;
		*elev = 2215.; /* for ocean horizon, not Andes! */
		printf("\n\n** Will use daylght time, Chilean date conventions. \n\n");
	}
	else if (obs_code[0] == 'c') {
		strcpy(site_name, "Las Campanas Observatory");
		strcpy(zone_name, "Chilean");
		*zabr = 'C';
		*use_dst = -1;
		*longit = 4.71333;
		*lat = -29.00833;
		*stdz = 4.;
		*elevsea = 2282.;
		*elev = 2282.; /* for ocean horizon, not Andes! */
		printf("\n\n** Will use daylght time, Chilean date conventions. \n\n");
	}
	else if (obs_code[0] == 'h') {
		strcpy(site_name, "Mount Hopkins, Arizona");
		strcpy(zone_name, "Mountain");
		*zabr = 'M';
		*use_dst = 0;
		*longit = 7.39233;
		*lat = 31.6883;
		*elevsea = 2608.;
		*elev = 500.;  /* approximate elevation above horizon mtns */
		*stdz = 7.;
	}
/*      else if (obs_code[0] == 'c') {
		strcpy(site_name,"Harvard College Observatory");
		strcpy(zone_name,"Eastern");
		*zabr = 'E';
		*use_dst = 1;
		*longit = 4.742;
		*lat = 42.38;
		*elevsea = 0.;  /* small, anyhow *?
		*elev = 0.;
		*stdz = 5.;
	}                                      --- COMMENTED OUT */
        else if (obs_code[0] == 'o') {
		strcpy(site_name,"McDonald Observatory");
		strcpy(zone_name,"Central");
		*zabr = 'C';
		*use_dst = 1;
		*longit = 6.93478;
                *lat = 30.6717;
                *elevsea = 2075;
                *elev = 1000.;  /* who knows? */
                *stdz = 6.;
        }
	else if (obs_code[0] == 'a') {
		strcpy(site_name, "Anglo-Australian Tel., Siding Spring");
		strcpy(zone_name, "Australian");
		*zabr = 'A';
		*use_dst = -2;
		*longit = -9.937739;
		*lat = -31.277039;
		*elevsea = 1149.;
		*elev = 670.;
		*stdz = -10.;
	}
/*	else if (obs_code[0] == 'b') {
		strcpy(site_name, "Black Moshannon Observatory");
		strcpy(zone_name, "Eastern");
		*zabr = 'E';
		*use_dst = 1;
		*longit = 5.20033;
		*lat = 40.92167;
		*elevsea = 738.;
		*elev = 0.;  
		*stdz = 5.;
	}
*/
	else if (obs_code[0] == 'd') {
		strcpy(site_name, "DAO, Victoria, BC");
		strcpy(zone_name, "Pacific");
		*zabr = 'P';
		*use_dst = 1;
		printf("\n\nWARNING: United States conventions for DST assumed.\n\n");
		*longit = 8.22778;
		*lat = 48.52;
		*elevsea = 74.;
		*elev = 74.;  /* not that it makes much difference */
		*stdz = 8.;             
	}
	else if (obs_code[0] == 'm') {
		strcpy(site_name, "Mauna Kea, Hawaii");
		strcpy(zone_name, "Hawaiian");
		*zabr = 'H';
		*use_dst = 0;
		*longit = 10.36478;
		*lat = 19.8267;
		*elevsea = 4215.;
		*elev = 4215.;  /* yow! */
		*stdz = 10.;
	}
	else if (obs_code[0] == 'l') {
		strcpy(site_name, "Lick Observatory");
		strcpy(zone_name, "Pacific");
		*zabr = 'P';
		*use_dst = 1;
		*longit = 8.10911;
		*lat = 37.3433;
		*elevsea = 1290.;
		*elev = 1290.; /* for those nice Pacific sunsets */
		*stdz = 8.;
	}
	else if (obs_code[0] == 'r') {
		strcpy(site_name, "Roque de los Muchachos");
		strcpy(zone_name, "pseudo-Greenwich");
		*zabr = 'G';
		*use_dst = 2;
		*longit = 1.192;
		*lat = 28.75833;
		*elevsea = 2326.;
		*elev = 2326.;
		*stdz = 0.;
	}
	else if (obs_code[0] == 'n') {
		printf("Enter new site parameters; the prompts give current values.\n");
		printf("(Even if current value is correct you must re-enter explicitly.)\n");
		printf("WEST longitude, (HOURS min sec); current value ");
		put_coords(*longit,3,1);
		printf(": ");
		*longit = get_coord();
		printf("Latitude, (d m s); current value ");
		put_coords(*lat,2,1);
		printf(": ");
		*lat = get_coord();
		printf("Actual elevation (meters) above sea level,");
		printf(" currently %5.0f:",*elevsea);
		strcpy(errprompt,"Because of error,");
		    /* situation is uncomplicated, so simple errprompt */
		getdouble(elevsea,-1000.,100000.,errprompt);
		printf("Effective elevation, meters (for rise/set),");
		printf(" currently %5.0f:",*elev); 
		getdouble(elev,-1000.,20000.,errprompt); /* limits of approx. ... */
		printf("Site name (< 30 char): ");
		nch=get_line(site_name);
		printf("Std time zone, hours W; currently %3.0f :",*stdz);
		getdouble(stdz,-13.,13.,errprompt);
		printf("Time zone name, e. g., Central: ");
		nch = get_line(zone_name);
		printf("Single-character zone abbreviation, currently %c : ",*zabr);
		scanf("%c",zabr);
		printf("Type daylight savings time option --- \n");
		printf("   0  ... don't use it, \n");
		printf("   1  ... use United States convention for clock change.\n");
		printf("   2  ... use Spanish (Continental?) convention.\n");
		printf("  -1  ... use Chilean convention.\n");
		printf("  -2  ... use Australian convention.\n");
		printf("(Other options would require new code). Answer: --> ");
		getshort(use_dst,-100,100,errprompt);
	}
		else {
		printf("UNKNOWN SITE '%c' -- left as %s. Note input is case-sensitive.\n",
			      obs_code[0],site_name);
	} 
	/* now compute derived quantity "horiz" = depression of horizon.*/
	*horiz = sqrt(2. * *elev / EQUAT_RAD) * DEG_IN_RADIAN;   
}

double atan_circ(x,y)
                
	double x,y;
                
{
	/* returns radian angle 0 to 2pi for coords x, y --
	   get that quadrant right !! */

	double theta;

	if((x == 0.) && (y == 0.)) return(0.);  /* guard ... */

	theta = atan2(y,x);  /* turns out there is such a thing in math.h */
	while(theta < 0.) theta += TWOPI;
	return(theta);
}

void min_max_alt(lat,dec,min,max)

	double lat,dec,*min,*max;

{
	/* computes minimum and maximum altitude for a given dec and
	    latitude. */

	double x;
	lat = lat / DEG_IN_RADIAN; /* pass by value! */
	dec = dec / DEG_IN_RADIAN;
	x = cos(dec)*cos(lat) + sin(dec)*sin(lat);
	if(fabs(x) <= 1.) {
		*max = asin(x) * DEG_IN_RADIAN;
	}
	else oprntf("Error in min_max_alt -- arcsin(>1)\n");
	x = sin(dec)*sin(lat) - cos(dec)*cos(lat);
	if(fabs(x) <= 1.) {
		*min = asin(x) * DEG_IN_RADIAN;
	}
	else oprntf("Error in min_max_alt -- arcsin(>1)\n");
}

double altit(dec,ha,lat,az,parang)
            
	double dec,ha,lat,*az,*parang;
            
/* returns altitude(degr) for dec, ha, lat (decimal degr, hr, degr); 
    also computes and returns azimuth through pointer argument,
    and as an extra added bonus returns parallactic angle (decimal degr)
    through another pointer argument. */
{
	double x,y,z;
	double sinp, cosp;  /* sin and cos of parallactic angle */
	double cosdec, sindec, cosha, sinha, coslat, sinlat;
			/* time-savers ... */

	dec = dec / DEG_IN_RADIAN;
	ha = ha / HRS_IN_RADIAN;
	lat = lat / DEG_IN_RADIAN;  /* thank heavens for pass-by-value */
	cosdec = cos(dec); sindec = sin(dec);
	cosha = cos(ha); sinha = sin(ha);
	coslat = cos(lat); sinlat = sin(lat);
	x = DEG_IN_RADIAN * asin(cosdec*cosha*coslat + sindec*sinlat);
	y =  sindec*coslat - cosdec*cosha*sinlat; /* due N comp. */
	z =  -1. * cosdec*sinha; /* due east comp. */
	*az = atan2(z,y);   

	/* as it turns out, having knowledge of the altitude and 
           azimuth makes the spherical trig of the parallactic angle
           less ambiguous ... so do it here!  Method uses the 
	   "astronomical triangle" connecting celestial pole, object,
           and zenith ... now know all the other sides and angles,
           so we can crush it ... */
	
	if(cosdec != 0.) { /* protect divide by zero ... */ 
	   sinp = -1. * sin(*az) * coslat / cosdec;
		/* spherical law of sines .. note cosdec = sin of codec,
			coslat = sin of colat .... */
	   cosp = -1. * cos(*az) * cosha - sin(*az) * sinha * sinlat;
		/* spherical law of cosines ... also transformed to local
                      available variables. */
	   *parang = atan2(sinp,cosp) * DEG_IN_RADIAN;
		/* let the library function find the quadrant ... */
	}
	else { /* you're on the pole */
	   if(lat >= 0.) *parang = 180.;
	   else *parang = 0.;
	}

	*az *= DEG_IN_RADIAN;  /* done with taking trig functions of it ... */ 
	while(*az < 0.) *az += 360.;  /* force 0 -> 360 */
	while(*az >= 360.) *az -= 360.;
	
	return(x);
}

double secant_z(alt)
	double alt;
{
	/* Computes the secant of z, assuming the object is not
           too low to the horizon; returns 100. if the object is
           low but above the horizon, -100. if the object is just
           below the horizon. */

	double secz;
	if(alt != 0) secz = 1. / sin(alt / DEG_IN_RADIAN);
	else secz = 100.;
	if(secz > 100.) secz = 100.;
	if(secz < -100.) secz = -100.;
	return(secz);
} 
 
double true_airmass(secz) 
	double secz;

{
	/* returns the true airmass for a given secant z. */
	/* The expression used is based on a tabulation of the mean KPNO 
           atmosphere given by C. M. Snell & A. M. Heiser, 1968,
	   PASP, 80, 336.  They tabulated the airmass at 5 degr 
           intervals from z = 60 to 85 degrees; I fit the data with 
           a fourth order poly for (secz - airmass) as a function of
           (secz - 1) using the IRAF curfit routine, then adjusted the
           zeroth order term to force (secz - airmass) to zero at
           z = 0.  The poly fit is very close to the tabulated points
	   (largest difference is 3.2e-4) and appears smooth.  
           This 85-degree point is at secz = 11.47, so for secz > 12
           I just return secz - 1.5 ... about the largest offset 
           properly determined. */

	double seczmin1;
	int i, ord = 4;
	double coef[5];
	double result = 0;
	
	coef[1] = 2.879465E-3;  /* sun compilers do not allow automatic
			initializations of arrays. */
        coef[2] = 3.033104E-3; 	
	coef[3] = 1.351167E-3; 
	coef[4] = -4.716679E-5;
	if(secz < 0.) return(-1.);  /* out of range. */
	if(secz > 12) return (secz - 1.5);  /* shouldn't happen .... */
	seczmin1 = secz - 1.;
	/* evaluate polynomial ... */
	for(i = ord; i > 0; i--) 
		result = (result + coef[i]) * seczmin1;
	/* no zeroth order term. */
	result = secz - result;
	return(result);			

}

double ha_alt(dec,lat,alt)
       
	double dec,lat,alt;

{
	/* returns hour angle at which object at dec is at altitude alt.
	   If object is never at this altitude, signals with special 
	   return values 1000 (always higher) and -1000 (always lower). */
	
	double x,coalt,min,max;
	
	min_max_alt(lat,dec,&min,&max);
	if(alt < min) 
		return(1000.);  /* flag value - always higher than asked */
	if(alt > max)
		return(-1000.); /* flag for object always lower than asked */
	dec = PI_OVER_2 - dec / DEG_IN_RADIAN;
	lat = PI_OVER_2 - lat / DEG_IN_RADIAN;
	coalt = PI_OVER_2 - alt / DEG_IN_RADIAN;
	x = (cos(coalt) - cos(dec)*cos(lat)) / (sin(dec)*sin(lat));
	if(fabs(x) <= 1.) return(acos(x) * HRS_IN_RADIAN);
	else {
		oprntf("Error in ha_alt ... acos(>1).\n");
		return(1000.);
	}
}

double subtend(ra1,dec1,ra2,dec2)
              
	double ra1,dec1,ra2,dec2;
              
  /*args in dec hrs and dec degrees */

{
	/* angle subtended by two positions in the sky --
	   return value is in radians.  Hybrid algorithm works down
	   to zero separation except very near the poles. */

	double x1, y1, z1, x2, y2, z2;
	double theta;

	ra1 = ra1 / HRS_IN_RADIAN;
	dec1 = dec1 / DEG_IN_RADIAN;
	ra2 = ra2 / HRS_IN_RADIAN;
	dec2 = dec2 / DEG_IN_RADIAN;
	x1 = cos(ra1)*cos(dec1);
	y1 = sin(ra1)*cos(dec1);
	z1 = sin(dec1);
	x2 = cos(ra2)*cos(dec2);
	y2 = sin(ra2)*cos(dec2);
	z2 = sin(dec2);
	theta = acos(x1*x2+y1*y2+z1*z2);
     /* use flat Pythagorean approximation if the angle is very small
	*and* you're not close to the pole; avoids roundoff in arccos. */
	if(theta < 1.0e-5) {  /* seldom the case, so don't combine test */
		if(fabs(dec1) < (PI_OVER_2 - 0.001) && 
		    fabs(dec2) < (PI_OVER_2 - 0.001))    {
			/* recycled variables here... */
			x1 = (ra2 - ra1) * cos((dec1+dec2)/2.);
			x2 = dec2 - dec1;
			theta = sqrt(x1*x1 + x2*x2);
		}
	}
	return(theta);
}

int get_pm(dec, mura, mudec) 

double dec, *mura, *mudec;

{
	/* This gets the proper motions.  New routine
	(2/94) assumes that primary PM convention will
	be rate of change of right ascension in seconds
	of time per year.  Either one can be entered here,
	but the value passed out is seconds of time per year at the
        equator (i.e., rate of change of RA itself). */

	char pmtype[3];
        int status;

	printf("Note -- two different conventions for RA proper motion.\n");
	printf("Enter RA p.m., either as delta RA(sec) or arcsec per yr.:");
	scanf("%lf",mura);
        if(*mura != 0.) {
	  printf("Type s if this is RA change in time sec per yr,\n");
	  printf("or a if this is motion in arcsec per yr:");
	  scanf("%s",pmtype);
        }
	else pmtype[0] = 's';  /* if pm is zero, it doesn't matter. */

	if(pmtype[0] == 's') status = 0;
	else if (pmtype[0] == 'a') {
		*mura = *mura /(15. * cos(dec / DEG_IN_RADIAN));
		printf("Equivalent to %8.5f sec of time per yr. at current dec.\n",
		     *mura);
                               printf("(Will only be correct at this dec.)\n");
                status = 0;
	}
	else {
		printf("UNRECOGNIZED PM TYPE ... defaults to sec/yr\n");
		status = -1;
        }
	printf("Give declination PM in arcsec/yr:");
	scanf("%lf",mudec);
	printf("(Note: Proper motion correction only ");
	printf("by mu * delta t; inaccurate near pole.)\n");
	printf("Don't forget to reset for new object.\n");

        return(status);  /* 0 = success */
}


int get_date(date) 

	struct date_time *date;
{	
	
	int valid_date = 0;
	int niter = 0;

	while(valid_date == 0 && niter++ < 3) {        
	/* have to put in a trap or there's a possible runaway ... */
	 	getshort(&(date->y),-10,2100,      
           	   "Give year again, then month and day.\n");
		if(date->y <= 0) return(-1);

	/* scan for mo and day here, *then* error check. */
		getshort(&(date->mo),1,12,
		   "Give month again (as number 1-12), then day.\n");
		getshort(&(date->d),0,32,"Give day again.\n");
		/* a lot of this error checking is redundant with the
		    checks in the new getshort routine.... */
		if(date->y < 100)  {
			date->y = date->y + 1900;
			printf("(Your answer assumed to mean %d)\n",date->y);
		}
		else if((date->y <= 1900 ) | (date->y >= 2100)){
			printf("Date out of range: only 1901 -> 2099 allowed.\n");
			printf("Try it again!\n");
		}
		/* might be nice to allow weird input dates, but calendrical
		   routines will not necessarily handle these right ... */
		else if((date->d < 0) || (date->d > 32)) 
			printf("day-of-month %d not allowed -- try again!\n",
				date->d);
		else if((date->mo < 1) || (date->mo > 12)) 
			printf("month %d doesn't exist -- try again!\n",
				date->mo);
		else {
			valid_date = 1;
			return(0);  /* success */
		}
	}
	if(valid_date == 0) return(-1);
}

int get_time(date) 

	struct date_time *date;

{
	
	struct coord ttime;
	
	dec_to_bab(get_coord(),&ttime);
	date->h = (short) ttime.hh;
		/* awkward, because h and m of date are short. */
	date->mn = (short) ttime.mm;
	date->s = ttime.ss;
	return(0);
}

double date_to_jd(struct date_time date) 

/* From Meeus' Astronomical Formulae for Calculators.  The two JD
   conversion routines routines were replaced 1998 November 29 to
   avoid inclusion of copyrighted "Numerical Recipes" code.  A test
   of 1 million random JDs between 1585 and 3200 AD gave the same
   conversions as the NR routines. */

{
	double jd;
 	int y, m;
	long A, B;

	if(date.mo <= 2) {
		y = date.y - 1;
		m = date.mo + 12;
	}
	else {
		y = date.y;
		m = date.mo;
	}

	A = (long) (y / 100.);
	B = 2 - A + (long) (A / 4.);

	jd = (long) (365.25 * y) + (long) (30.6001 * (m + 1)) + date.d + 
		1720994.5;

	jd += date.h / 24. + date.mn / 1440. + date.s / 86400.;

	if(date.y > 1583) return(jd + B);  
		else return(jd);
		/* Not quite right, since Gregorian calendar first
		adopted around Oct 1582.  But fine for modern. */
}

void caldat(jdin, date, dow)

double jdin;
struct date_time *date;
short *dow;  /* day of week */

{
	/* from Jean Meeus, Astronomical Formulae for Calculators,
	   published by Willman-Bell Inc.
           Avoids a copyrighted routine from Numerical Recipes.
           Tested and works properly from the beginning of the 
           Gregorian calendar era (1583) to beyond 3000 AD. */

	double jdtmp;
	long alpha;
	long Z;
	long A, B, C, D, E;
	double F; 
	double x;   /* for day-of-week calculation */

	int rounded_ok = 0;

	while(rounded_ok == 0) {  /* needed to handle round-ups to 24h */
		jdtmp = jdin + 0.5;
	
		Z = (long) jdtmp;
	
	        x = Z/7.+0.01;
	        *dow = 7.*(x - (long) x);   /* truncate for day of week */
	
		F = jdtmp - Z;
	
		if(Z < 2299161) A = Z;
		else {
			alpha = (long) ((Z - 1867216.25) / 36524.25);
			A = Z + 1 + alpha - (long) (alpha / 4);
		}
	
		B = A + 1524;
		C = ((B - 122.1) / 365.25);
		D =  (365.25 * C);\
		E =  ((B - D) / 30.6001);
	
		date->d = B - D - (long)(30.6001 * E);
		if(E < 13.5) date->mo = E - 1;
			else date->mo = E - 13;
		if(date->mo  > 2.5)  date->y = C - 4716;
			else date->y = C - 4715;
		
		date->h = F * 24.;  /* truncate */
	        date->mn = (F - ((float) date->h)/24.) * 1440.;
	        date->s = (F - ((float) date->h)/24. -
	                        ((float) date->mn)/1440.) * 86400;
		if (date->h < 24.) rounded_ok = 1;
		else {
			printf("rounding up from 24 h\n");
			jdin = jdin + JDRESOLUTION;
		}
	}
}

short day_of_week(jd)
	double jd;
{ 
	/* returns day of week for a jd, 0 = Mon, 6 = Sun. */

	double x,y;
	long i;
	short d;
	
	jd = jd+0.5;
	i = jd; /* truncate */
	x = i/7.+0.01; 
	d = 7.*(x - (long) x);   /* truncate */
	return(d);
}

double day_of_year(jd) 
	double jd;

{
	double jdjan0;
	struct date_time date;
	short dow;
	
	caldat(jd,&date,&dow);
	/* find jd of "jan 0" = Dec 31 of previous year */
	date.y = date.y - 1;
	date.mo = 12;
	date.d = 31;
	date.h = 0;
	date.mn = 0;
	date.s = 0.;
	jdjan0 = date_to_jd(date);
	return(jd - jdjan0);
}

void print_day(d)
	short d;

{
	/* prints day of week given number 0=Mon,6=Sun */
	char *days = "MonTueWedThuFriSatSun";
	char day_out[4];

	day_out[0] = *(days+3*d);
	day_out[1] = *(days+3*d+1);
	day_out[2] = *(days+3*d+2);
	day_out[3] = '\0';  /* terminate with null char */
	
	oprntf("%s",day_out);
}


void print_all(jdin)

	double jdin;
{
	/* given a julian date,
	prints a year, month, day, hour, minute, second */

	struct date_time date;
        int ytemp, dtemp; /* compiler bug workaround ... SUN
         and silicon graphics */
	char *months = "JanFebMarAprMayJunJulAugSepOctNovDec";
	char mo_out[4];
	double out_time;
	short dow;

	caldat(jdin,&date,&dow);

	print_day(dow);
	oprntf(", ");

	mo_out[0] = *(months + 3*(date.mo - 1));
	mo_out[1] = *(months + 3*(date.mo - 1) + 1);
	mo_out[2] = *(months + 3*(date.mo - 1) + 2);
	mo_out[3] = '\0';

	/* going through the rigamarole to avoid 60's */

	out_time = date.h + date.mn / 60. + date.s / 3600.;

	ytemp = (int) date.y;
	dtemp = (int) date.d;
	oprntf("%d %s %2d, time ",
		ytemp,mo_out,dtemp);
	put_coords(out_time,3,0);
}

void print_current(date,night_date,enter_ut) 
        struct date_time date;
	short night_date, enter_ut;
{
	/* just prints out the date & time and a little statement
           of whether time is "local" or "ut".  Functionalized to 
           compactify some later code. */
	
	double jd;

	jd = date_to_jd(date);
        if((night_date == 1) && (date.h < 12)) jd = jd + 1.0;
	print_all(jd);
	if(enter_ut == 0) oprntf(" local time.");
	else oprntf(" Universal time.");
}

void print_calendar(jdin,dow)
                   
	double jdin;
	short *dow;
                   
{
	/* given a julian date prints a year, month, day.
	   Returns day of week (0 - 6) for optional printing ---
	   idea is to do one roundoff only to eliminate possibility
	   of day/date disagreement. */
	
	struct date_time date;
	char *months = "JanFebMarAprMayJunJulAugSepOctNovDec";
	char mo_out[4]; 
	int ytemp, dtemp;  /* compiler bug workaround -- SUN and
           Silicon Graphics machines */

	caldat(jdin,&date,dow); 
	mo_out[0] = *(months + 3*(date.mo - 1));
	mo_out[1] = *(months + 3*(date.mo - 1) + 1);
	mo_out[2] = *(months + 3*(date.mo - 1) + 2);
	mo_out[3] = '\0';
        ytemp = (int) date.y;
        dtemp = (int) date.d;
	oprntf("%d %s %d",ytemp,mo_out,dtemp);
}

void print_time(jdin,prec)

	double jdin;
	short prec;

{
	/* given a julian day, prints time only; 
	   special precision of "-1" prints only hours!  */
	struct date_time date;
	double temptime;
	short dow;
	int temp; /* to get around a Sun compiler bug ... */
	
	caldat(jdin,&date,&dow);
	temp = date.h;  /* cast from short to int */

	temptime = date.h + date.mn/60. + date.s/3600.;

	if(prec >= 0) put_coords(temptime,prec,0);
	else if(date.mn < 30) oprntf("%2.0d hr",temp);
	else oprntf("%2.0d hr",(temp+1)); /* round it up */
}       

double frac_part(x) 

	double x;
{
	long i;
	i = x;
	x = x - i;
	return(x);
}
	

double lst(jd,longit) 
           
	double jd,longit; 
           
{
	/* returns the local MEAN sidereal time (dec hrs) at julian date jd
	   at west longitude long (decimal hours).  Follows
	   definitions in 1992 Astronomical Almanac, pp. B7 and L2. 
	   Expression for GMST at 0h ut referenced to Aoki et al, A&A 105,
	   p.359, 1982.  On workstations, accuracy (numerical only!)
	   is about a millisecond in the 1990s. */

	double t, ut, jdmid, jdint, jdfrac, sid_g, sid;
	long jdin, sid_int;

	jdin = jd;         /* fossil code from earlier package which 
			split jd into integer and fractional parts ... */
	jdint = jdin;
	jdfrac = jd - jdint;
	if(jdfrac < 0.5) {
		jdmid = jdint - 0.5;
		ut = jdfrac + 0.5;
	}
	else {
		jdmid = jdint + 0.5;
		ut = jdfrac - 0.5;
	}
	t = (jdmid - J2000)/36525;
	sid_g = (24110.54841+8640184.812866*t+0.093104*t*t-6.2e-6*t*t*t)/SEC_IN_DAY;
	sid_int = sid_g;
	sid_g = sid_g - (double) sid_int;
	sid_g = sid_g + 1.0027379093 * ut - longit/24.;
	sid_int = sid_g;
	sid_g = (sid_g - (double) sid_int) * 24.;
	if(sid_g < 0.) sid_g = sid_g + 24.;
	return(sid_g);
}

double adj_time(x)
	double x;

{
	/* adjusts a time (decimal hours) to be between -12 and 12, 
	   generally used for hour angles.  */

	if(fabs(x) < 100000.) {  /* too inefficient for this! */
		while(x > 12.) {
			x = x - 24.;
		}
		while(x < -12.) {
			x = x + 24.;
		}
	}
	else oprntf("Out of bounds in adj_time!\n");
	return(x);
}

void lpmoon(jd,lat,sid,ra,dec,dist)  

	double jd,lat,sid,*ra,*dec,*dist;  

/* implements "low precision" moon algorithms from
   Astronomical Almanac (p. D46 in 1992 version).  Does
   apply the topocentric correction. 
Units are as follows
jd,lat, sid;   decimal hours 
*ra, *dec,   decimal hours, degrees 
	*dist;      earth radii */
{

	double T, lambda, beta, pie, l, m, n, x, y, z, alpha, delta,
		rad_lat, rad_lst, distance, topo_dist;
	char dummy[40];  /* to fix compiler bug on IBM system */

	T = (jd - J2000) / 36525.;  /* jul cent. since J2000.0 */

	lambda = 218.32 + 481267.883 * T 
	   + 6.29 * sin((134.9 + 477198.85 * T) / DEG_IN_RADIAN)
	   - 1.27 * sin((259.2 - 413335.38 * T) / DEG_IN_RADIAN)
	   + 0.66 * sin((235.7 + 890534.23 * T) / DEG_IN_RADIAN)
	   + 0.21 * sin((269.9 + 954397.70 * T) / DEG_IN_RADIAN)
	   - 0.19 * sin((357.5 + 35999.05 * T) / DEG_IN_RADIAN)
	   - 0.11 * sin((186.6 + 966404.05 * T) / DEG_IN_RADIAN);
	lambda = lambda / DEG_IN_RADIAN;
	beta = 5.13 * sin((93.3 + 483202.03 * T) / DEG_IN_RADIAN)
	   + 0.28 * sin((228.2 + 960400.87 * T) / DEG_IN_RADIAN)
	   - 0.28 * sin((318.3 + 6003.18 * T) / DEG_IN_RADIAN)
	   - 0.17 * sin((217.6 - 407332.20 * T) / DEG_IN_RADIAN);
	beta = beta / DEG_IN_RADIAN;
	pie = 0.9508 
	   + 0.0518 * cos((134.9 + 477198.85 * T) / DEG_IN_RADIAN)
	   + 0.0095 * cos((259.2 - 413335.38 * T) / DEG_IN_RADIAN)
	   + 0.0078 * cos((235.7 + 890534.23 * T) / DEG_IN_RADIAN)
	   + 0.0028 * cos((269.9 + 954397.70 * T) / DEG_IN_RADIAN);
	pie = pie / DEG_IN_RADIAN;
	distance = 1 / sin(pie);

	l = cos(beta) * cos(lambda);
	m = 0.9175 * cos(beta) * sin(lambda) - 0.3978 * sin(beta);
	n = 0.3978 * cos(beta) * sin(lambda) + 0.9175 * sin(beta);

	x = l * distance; 
	y = m * distance; 
	z = n * distance;  /* for topocentric correction */
	/* lat isn't passed right on some IBM systems unless you do this
	   or something like it! */
	sprintf(dummy,"%f",lat);
	rad_lat = lat / DEG_IN_RADIAN;
	rad_lst = sid / HRS_IN_RADIAN;
	x = x - cos(rad_lat) * cos(rad_lst);
	y = y - cos(rad_lat) * sin(rad_lst);
	z = z - sin(rad_lat);


	topo_dist = sqrt(x * x + y * y + z * z);

	l = x / topo_dist; 
	m = y / topo_dist; 
	n = z / topo_dist;

	alpha = atan_circ(l,m);
	delta = asin(n);
	*ra = alpha * HRS_IN_RADIAN;
	*dec = delta * DEG_IN_RADIAN;
	*dist = topo_dist;
}


void lpsun(jd,ra,dec)

	double jd, *ra, *dec;

/* Low precision formulae for the sun, from Almanac p. C24 (1990) */
/* ra and dec are returned as decimal hours and decimal degrees. */

{
	double n, L, g, lambda,epsilon,alpha,delta,x,y,z;

	n = jd - J2000;
	L = 280.460 + 0.9856474 * n;
	g = (357.528 + 0.9856003 * n)/DEG_IN_RADIAN;
	lambda = (L + 1.915 * sin(g) + 0.020 * sin(2. * g))/DEG_IN_RADIAN;
	epsilon = (23.439 - 0.0000004 * n)/DEG_IN_RADIAN;

	x = cos(lambda); 
	y = cos(epsilon) * sin(lambda); 
	z = sin(epsilon)*sin(lambda);

	*ra = (atan_circ(x,y))*HRS_IN_RADIAN;
	*dec = (asin(z))*DEG_IN_RADIAN;
}

void eclrot(jd, x, y, z)
           
	double jd, *x, *y, *z;
           
/* rotates ecliptic rectangular coords x, y, z to
   equatorial (all assumed of date.) */

{
	double incl;
	double xpr,ypr,zpr;
	double T;

	T = (jd - J2000) / 36525;  /* centuries since J2000 */
	
	incl = (23.439291 + T * (-0.0130042 - 0.00000016 * T))/DEG_IN_RADIAN; 
		/* 1992 Astron Almanac, p. B18, dropping the 
		   cubic term, which is 2 milli-arcsec! */
	ypr = cos(incl) * *y - sin(incl) * *z;
	zpr = sin(incl) * *y + cos(incl) * *z;
	*y = ypr;
	*z = zpr;
	/* x remains the same. */       
}

double circulo(x)
	double x;
{
	/* assuming x is an angle in degrees, returns 
	   modulo 360 degrees. */

	int n;

	n = (int)(x / 360.);
	return(x - 360. * n);
}       


void geocent(geolong,geolat,height,x_geo,y_geo,z_geo)
            
	double geolong, geolat, height, *x_geo, *y_geo, *z_geo;
            
/* computes the geocentric coordinates from the geodetic 
(standard map-type) longitude, latitude, and height. 
These are assumed to be in decimal hours, decimal degrees, and
meters respectively.  Notation generally follows 1992 Astr Almanac, 
p. K11 */


{
	
	double denom, C_geo, S_geo;

	geolat = geolat / DEG_IN_RADIAN;
	geolong = geolong / HRS_IN_RADIAN;      
	denom = (1. - FLATTEN) * sin(geolat);
	denom = cos(geolat) * cos(geolat) + denom*denom;
	C_geo = 1. / sqrt(denom);
	S_geo = (1. - FLATTEN) * (1. - FLATTEN) * C_geo;
	C_geo = C_geo + height / EQUAT_RAD;  /* deviation from almanac
		       notation -- include height here. */
	S_geo = S_geo + height / EQUAT_RAD;
	*x_geo = C_geo * cos(geolat) * cos(geolong);
	*y_geo = C_geo * cos(geolat) * sin(geolong);
	*z_geo = S_geo * sin(geolat);
}
 
double etcorr(jd)

double jd;

{

	/* Given a julian date in 1900-2100, returns the correction
           delta t which is:
		TDT - UT (after 1983 and before 1998)
		ET - UT (before 1983)
		an extrapolated guess  (after 2001). 

	For dates in the past (<= 2001 and after 1900) the value is linearly
        interpolated on 5-year intervals; for dates after the present,
        an extrapolation is used, because the true value of delta t
	cannot be predicted precisely.  Note that TDT is essentially the
	modern version of ephemeris time with a slightly cleaner 
	definition.  

	Where the algorithm shifts there will be a small (< 0.1 sec)
        discontinuity.  Also, the 5-year linear interpolation scheme can 
        lead to errors as large as 0.5 seconds in some cases, though
 	usually rather smaller.   One seldom has actual UT to work with anyway,
	since the commonly-used UTC is tied to TAI within an integer number
	of seconds.  */

	double jd1900 = 2415019.5;
	double dates[22];
	double delts[22];  /* can't initialize this look-up table
            with stupid old sun compiler .... */
	double year, delt;
	int i;
	
	/* this stupid patch for primitive sun C compilers .... 
		do not allow automatic initialization of arrays! */

	for(i = 0; i <= 20; i++) dates[i] = 1900 + (double) i * 5.;
	dates[21] = 2001.;  /* the last accurately tabulated one in the
		2003 Almanac ... */

	delts[0] = -2.72;  delts[1] = 3.86; delts[2] = 10.46;
	delts[3] = 17.20;  delts[4] = 21.16; delts[5] = 23.62;
	delts[6] = 24.02;  delts[7] = 23.93; delts[8] = 24.33;
	delts[9] = 26.77;  delts[10] = 29.15; delts[11] = 31.07;
	delts[12] = 33.15;  delts[13] = 35.73; delts[14] = 40.18;
	delts[15] = 45.48;  delts[16] = 50.54; delts[17] = 54.34;
	delts[18] = 56.86;  delts[19] = 60.78; delts[20] = 63.83;
	delts[21] = 64.09;

	year = 1900. + (jd - jd1900) / 365.25;

	if(year < 2001. && year >= 1900.) {
		i = (year - 1900) / 5;
		delt = delts[i] + 
		 ((delts[i+1] - delts[i])/(dates[i+1] - dates[i])) * (year - dates[i]);
	}

	else if (year >= 2001. && year < 2100.)
		delt = 31.69 + (2.164e-3) * (jd - 2436935.4);  /* rough extrapolation */
                /* the 31.69 is adjusted to give 64.09 sec at the start of 2001. */
	else if (year < 1900) {
		printf("etcorr ... no ephemeris time data for < 1900.\n");
       		delt = 0.;
	}

	else if (year >= 2100.) {
		printf("etcorr .. very long extrapolation in delta T - inaccurate.\n");
		delt = 180.; /* who knows? */
	} 

	return(delt);
}

void accumoon(jd,geolat,lst,elevsea,geora,geodec,geodist,
     topora,topodec,topodist)
             
	double jd,geolat,lst,elevsea;
     	double *geora,*geodec,*geodist,*topora,*topodec,*topodist;
             
  /* jd, dec. degr., dec. hrs., meters */
/* More accurate (but more elaborate and slower) lunar 
   ephemeris, from Jean Meeus' *Astronomical Formulae For Calculators*,
   pub. Willman-Bell.  Includes all the terms given there. */

{       
/*      double *eclatit,*eclongit, *pie,*ra,*dec,*dist; geocent quantities,
		formerly handed out but not in this version */
	double pie, dist;  /* horiz parallax */
	double Lpr,M,Mpr,D,F,Om,T,Tsq,Tcb;
	double e,lambda,B,beta,om1,om2;
	double sinx, x, y, z, l, m, n;
	double x_geo, y_geo, z_geo;  /* geocentric position of *observer* */    

	jd = jd + etcorr(jd)/SEC_IN_DAY;   /* approximate correction to ephemeris time */
	T = (jd - 2415020.) / 36525.;   /* this based around 1900 ... */
	Tsq = T * T;
	Tcb = Tsq * T;

	Lpr = 270.434164 + 481267.8831 * T - 0.001133 * Tsq 
			+ 0.0000019 * Tcb;
	M = 358.475833 + 35999.0498*T - 0.000150*Tsq
			- 0.0000033*Tcb;
	Mpr = 296.104608 + 477198.8491*T + 0.009192*Tsq 
			+ 0.0000144*Tcb;
	D = 350.737486 + 445267.1142*T - 0.001436 * Tsq
			+ 0.0000019*Tcb;
	F = 11.250889 + 483202.0251*T -0.003211 * Tsq 
			- 0.0000003*Tcb;
	Om = 259.183275 - 1934.1420*T + 0.002078*Tsq 
			+ 0.0000022*Tcb;

	Lpr = circulo(Lpr);
	Mpr = circulo(Mpr);     
	M = circulo(M);
	D = circulo(D);
	F = circulo(F);
	Om = circulo(Om);

	
	sinx =  sin((51.2 + 20.2 * T)/DEG_IN_RADIAN);
	Lpr = Lpr + 0.000233 * sinx;
	M = M - 0.001778 * sinx;
	Mpr = Mpr + 0.000817 * sinx;
	D = D + 0.002011 * sinx;
	
	sinx = 0.003964 * sin((346.560+132.870*T -0.0091731*Tsq)/DEG_IN_RADIAN);

	Lpr = Lpr + sinx;
	Mpr = Mpr + sinx;
	D = D + sinx;
	F = F + sinx;


	sinx = sin(Om/DEG_IN_RADIAN);
	Lpr = Lpr + 0.001964 * sinx;
	Mpr = Mpr + 0.002541 * sinx;
	D = D + 0.001964 * sinx;
	F = F - 0.024691 * sinx;
	F = F - 0.004328 * sin((Om + 275.05 -2.30*T)/DEG_IN_RADIAN);

	e = 1 - 0.002495 * T - 0.00000752 * Tsq;

	M = M / DEG_IN_RADIAN;   /* these will all be arguments ... */
	Mpr = Mpr / DEG_IN_RADIAN;
	D = D / DEG_IN_RADIAN;
	F = F / DEG_IN_RADIAN;

	lambda = Lpr + 6.288750 * sin(Mpr)
		+ 1.274018 * sin(2*D - Mpr)
		+ 0.658309 * sin(2*D)
		+ 0.213616 * sin(2*Mpr)
		- e * 0.185596 * sin(M) 
		- 0.114336 * sin(2*F)
		+ 0.058793 * sin(2*D - 2*Mpr)
		+ e * 0.057212 * sin(2*D - M - Mpr)
		+ 0.053320 * sin(2*D + Mpr)
		+ e * 0.045874 * sin(2*D - M)
		+ e * 0.041024 * sin(Mpr - M)
		- 0.034718 * sin(D)
		- e * 0.030465 * sin(M+Mpr)
		+ 0.015326 * sin(2*D - 2*F)
		- 0.012528 * sin(2*F + Mpr)
		- 0.010980 * sin(2*F - Mpr)
		+ 0.010674 * sin(4*D - Mpr)
		+ 0.010034 * sin(3*Mpr)
		+ 0.008548 * sin(4*D - 2*Mpr)
		- e * 0.007910 * sin(M - Mpr + 2*D)
		- e * 0.006783 * sin(2*D + M)
		+ 0.005162 * sin(Mpr - D);

		/* And furthermore.....*/

	lambda = lambda + e * 0.005000 * sin(M + D)
		+ e * 0.004049 * sin(Mpr - M + 2*D)
		+ 0.003996 * sin(2*Mpr + 2*D)
		+ 0.003862 * sin(4*D)
		+ 0.003665 * sin(2*D - 3*Mpr)
		+ e * 0.002695 * sin(2*Mpr - M)
		+ 0.002602 * sin(Mpr - 2*F - 2*D)
		+ e * 0.002396 * sin(2*D - M - 2*Mpr)
		- 0.002349 * sin(Mpr + D)
		+ e * e * 0.002249 * sin(2*D - 2*M)
		- e * 0.002125 * sin(2*Mpr + M)
		- e * e * 0.002079 * sin(2*M)
		+ e * e * 0.002059 * sin(2*D - Mpr - 2*M)
		- 0.001773 * sin(Mpr + 2*D - 2*F)
		- 0.001595 * sin(2*F + 2*D)
		+ e * 0.001220 * sin(4*D - M - Mpr)
		- 0.001110 * sin(2*Mpr + 2*F)
		+ 0.000892 * sin(Mpr - 3*D)
		- e * 0.000811 * sin(M + Mpr + 2*D)
		+ e * 0.000761 * sin(4*D - M - 2*Mpr)
		+ e * e * 0.000717 * sin(Mpr - 2*M)
		+ e * e * 0.000704 * sin(Mpr - 2 * M - 2*D)
		+ e * 0.000693 * sin(M - 2*Mpr + 2*D)
		+ e * 0.000598 * sin(2*D - M - 2*F)
		+ 0.000550 * sin(Mpr + 4*D)
		+ 0.000538 * sin(4*Mpr)
		+ e * 0.000521 * sin(4*D - M)
		+ 0.000486 * sin(2*Mpr - D);
	
/*              *eclongit = lambda;  */

	B = 5.128189 * sin(F)
		+ 0.280606 * sin(Mpr + F)
		+ 0.277693 * sin(Mpr - F)
		+ 0.173238 * sin(2*D - F)
		+ 0.055413 * sin(2*D + F - Mpr)
		+ 0.046272 * sin(2*D - F - Mpr)
		+ 0.032573 * sin(2*D + F)
		+ 0.017198 * sin(2*Mpr + F)
		+ 0.009267 * sin(2*D + Mpr - F)
		+ 0.008823 * sin(2*Mpr - F)
		+ e * 0.008247 * sin(2*D - M - F) 
		+ 0.004323 * sin(2*D - F - 2*Mpr)
		+ 0.004200 * sin(2*D + F + Mpr)
		+ e * 0.003372 * sin(F - M - 2*D)
		+ 0.002472 * sin(2*D + F - M - Mpr)
		+ e * 0.002222 * sin(2*D + F - M)
		+ e * 0.002072 * sin(2*D - F - M - Mpr)
		+ e * 0.001877 * sin(F - M + Mpr)
		+ 0.001828 * sin(4*D - F - Mpr)
		- e * 0.001803 * sin(F + M)
		- 0.001750 * sin(3*F)
		+ e * 0.001570 * sin(Mpr - M - F)
		- 0.001487 * sin(F + D)
		- e * 0.001481 * sin(F + M + Mpr)
		+ e * 0.001417 * sin(F - M - Mpr)
		+ e * 0.001350 * sin(F - M)
		+ 0.001330 * sin(F - D)
		+ 0.001106 * sin(F + 3*Mpr)
		+ 0.001020 * sin(4*D - F)
		+ 0.000833 * sin(F + 4*D - Mpr);
     /* not only that, but */
	B = B + 0.000781 * sin(Mpr - 3*F)
		+ 0.000670 * sin(F + 4*D - 2*Mpr)
		+ 0.000606 * sin(2*D - 3*F)
		+ 0.000597 * sin(2*D + 2*Mpr - F)
		+ e * 0.000492 * sin(2*D + Mpr - M - F)
		+ 0.000450 * sin(2*Mpr - F - 2*D)
		+ 0.000439 * sin(3*Mpr - F)
		+ 0.000423 * sin(F + 2*D + 2*Mpr)
		+ 0.000422 * sin(2*D - F - 3*Mpr)
		- e * 0.000367 * sin(M + F + 2*D - Mpr)
		- e * 0.000353 * sin(M + F + 2*D)
		+ 0.000331 * sin(F + 4*D)
		+ e * 0.000317 * sin(2*D + F - M + Mpr)
		+ e * e * 0.000306 * sin(2*D - 2*M - F)
		- 0.000283 * sin(Mpr + 3*F);
	
	
	om1 = 0.0004664 * cos(Om/DEG_IN_RADIAN);        
	om2 = 0.0000754 * cos((Om + 275.05 - 2.30*T)/DEG_IN_RADIAN);
	
	beta = B * (1. - om1 - om2);
/*      *eclatit = beta; */
	
	pie = 0.950724 
		+ 0.051818 * cos(Mpr)
		+ 0.009531 * cos(2*D - Mpr)
		+ 0.007843 * cos(2*D)
		+ 0.002824 * cos(2*Mpr)
		+ 0.000857 * cos(2*D + Mpr)
		+ e * 0.000533 * cos(2*D - M)
		+ e * 0.000401 * cos(2*D - M - Mpr)
		+ e * 0.000320 * cos(Mpr - M)
		- 0.000271 * cos(D)
		- e * 0.000264 * cos(M + Mpr)
		- 0.000198 * cos(2*F - Mpr)
		+ 0.000173 * cos(3*Mpr)
		+ 0.000167 * cos(4*D - Mpr)
		- e * 0.000111 * cos(M)
		+ 0.000103 * cos(4*D - 2*Mpr)
		- 0.000084 * cos(2*Mpr - 2*D)
		- e * 0.000083 * cos(2*D + M)
		+ 0.000079 * cos(2*D + 2*Mpr)
		+ 0.000072 * cos(4*D)
		+ e * 0.000064 * cos(2*D - M + Mpr)
		- e * 0.000063 * cos(2*D + M - Mpr)
		+ e * 0.000041 * cos(M + D)
		+ e * 0.000035 * cos(2*Mpr - M)
		- 0.000033 * cos(3*Mpr - 2*D)
		- 0.000030 * cos(Mpr + D)
		- 0.000029 * cos(2*F - 2*D)
		- e * 0.000029 * cos(2*Mpr + M)
		+ e * e * 0.000026 * cos(2*D - 2*M)
		- 0.000023 * cos(2*F - 2*D + Mpr)
		+ e * 0.000019 * cos(4*D - M - Mpr);

	beta = beta/DEG_IN_RADIAN;
	lambda = lambda/DEG_IN_RADIAN;
	l = cos(lambda) * cos(beta);    
	m = sin(lambda) * cos(beta);
	n = sin(beta);
	eclrot(jd,&l,&m,&n);
	
	dist = 1/sin((pie)/DEG_IN_RADIAN);
	x = l * dist;
	y = m * dist;
	z = n * dist;

	*geora = atan_circ(l,m) * HRS_IN_RADIAN;
	*geodec = asin(n) * DEG_IN_RADIAN;        
	*geodist = dist;

	geocent(lst,geolat,elevsea,&x_geo,&y_geo,&z_geo);
	
	x = x - x_geo;  /* topocentric correction using elliptical earth fig. */
	y = y - y_geo;
	z = z - z_geo;

	*topodist = sqrt(x*x + y*y + z*z);
	
	l = x / (*topodist);
	m = y / (*topodist);
	n = z / (*topodist);

	*topora = atan_circ(l,m) * HRS_IN_RADIAN;
	*topodec = asin(n) * DEG_IN_RADIAN; 

}

void flmoon(n,nph,jdout) 
           
	int n,nph;
	double *jdout;
           
/* Gives jd (+- 2 min) of phase nph on lunation n; replaces
less accurate Numerical Recipes routine.  This routine 
implements formulae found in Jean Meeus' *Astronomical Formulae
for Calculators*, 2nd edition, Willman-Bell.  A very useful
book!! */

/* n, nph lunation and phase; nph = 0 new, 1 1st, 2 full, 3 last 
 *jdout   jd of requested phase */

{
	double jd, cor;
	double M, Mpr, F;
	double T;
	double lun;

	lun = (double) n + (double) nph / 4.;
	T = lun / 1236.85;
	jd = 2415020.75933 + 29.53058868 * lun  
		+ 0.0001178 * T * T 
		- 0.000000155 * T * T * T
		+ 0.00033 * sin((166.56 + 132.87 * T - 0.009173 * T * T)/DEG_IN_RADIAN);
	M = 359.2242 + 29.10535608 * lun - 0.0000333 * T * T - 0.00000347 * T * T * T;
	M = M / DEG_IN_RADIAN;
	Mpr = 306.0253 + 385.81691806 * lun + 0.0107306 * T * T + 0.00001236 * T * T * T;
	Mpr = Mpr / DEG_IN_RADIAN;
	F = 21.2964 + 390.67050646 * lun - 0.0016528 * T * T - 0.00000239 * T * T * T;
	F = F / DEG_IN_RADIAN;
	if((nph == 0) || (nph == 2)) {/* new or full */
		cor =   (0.1734 - 0.000393*T) * sin(M)
			+ 0.0021 * sin(2*M)
			- 0.4068 * sin(Mpr)
			+ 0.0161 * sin(2*Mpr)
			- 0.0004 * sin(3*Mpr)
			+ 0.0104 * sin(2*F)
			- 0.0051 * sin(M + Mpr)
			- 0.0074 * sin(M - Mpr)
			+ 0.0004 * sin(2*F+M)
			- 0.0004 * sin(2*F-M)
			- 0.0006 * sin(2*F+Mpr)
			+ 0.0010 * sin(2*F-Mpr)
			+ 0.0005 * sin(M+2*Mpr);
		jd = jd + cor;
	}
	else {
		cor = (0.1721 - 0.0004*T) * sin(M)
			+ 0.0021 * sin(2 * M)
			- 0.6280 * sin(Mpr)
			+ 0.0089 * sin(2 * Mpr)
			- 0.0004 * sin(3 * Mpr)
			+ 0.0079 * sin(2*F)
			- 0.0119 * sin(M + Mpr)
			- 0.0047 * sin(M - Mpr)
			+ 0.0003 * sin(2 * F + M)
			- 0.0004 * sin(2 * F - M)
			- 0.0006 * sin(2 * F + Mpr)
			+ 0.0021 * sin(2 * F - Mpr)
			+ 0.0003 * sin(M + 2 * Mpr)
			+ 0.0004 * sin(M - 2 * Mpr)
			- 0.0003 * sin(2*M + Mpr);
		if(nph == 1) cor = cor + 0.0028 - 
				0.0004 * cos(M) + 0.0003 * cos(Mpr);
		if(nph == 3) cor = cor - 0.0028 +
				0.0004 * cos(M) - 0.0003 * cos(Mpr);
		jd = jd + cor;

	}
	*jdout = jd;
}

float lun_age(jd, nlun) 
             
	double jd; 
	int *nlun; 
             
{
	/* compute age in days of moon since last new,
	   and lunation of last new moon. */

	int n; /* appropriate lunation */
	int nlast;
	double newjd, lastnewjd;
	short kount=0;
	float x;

	nlast = (jd - 2415020.5) / 29.5307 - 1;
	 
	flmoon(nlast,0,&lastnewjd);
	nlast++;
	flmoon(nlast,0,&newjd);
	while((newjd < jd) && (kount < 40)) {
		lastnewjd = newjd;
		nlast++;
		flmoon(nlast,0,&newjd);
	}
 	if(kount > 35) {
		oprntf("Didn't find phase in lun_age!\n");
		x = -10.;
                *nlun = 0;
	}
	else {
	  x = jd - lastnewjd;
	  *nlun = nlast - 1;
        }

	return(x);
}       

void print_phase(jd)
	double jd; 

{
	/* prints a verbal description of moon phase, given the
	   julian date.  */

	int n; /* appropriate lunation */
	int nlast, noctiles;
	double newjd, lastnewjd;
	double fqjd, fljd, lqjd;  /* jds of first, full, and last in this lun.*/
	short kount=0;
	float x;

	nlast = (jd - 2415020.5) / 29.5307 - 1;  /* find current lunation */
	 
	flmoon(nlast,0,&lastnewjd);
	nlast++;
	flmoon(nlast,0,&newjd);
	while((newjd < jd) && (kount < 40)) {
		lastnewjd = newjd;
		nlast++;
		flmoon(nlast,0,&newjd);
	}
	if(kount > 35) {  /* oops ... didn't find it ... */
		oprntf("Didn't find phase in print_phase!\n");
		x = -10.;
	}
	else {     /* found lunation ok */
		x = jd - lastnewjd;
		nlast--;
		noctiles = x / 3.69134;  /* 3.69134 = 1/8 month; truncate. */
		if(noctiles == 0) oprntf("%3.1f days since new moon",x);
		else if (noctiles <= 2) {  /* nearest first quarter */
			flmoon(nlast,1,&fqjd);
			x = jd - fqjd;
			if(x < 0.)
			  oprntf("%3.1f days before first quarter",(-1.*x));
			else
			  oprntf("%3.1f days since first quarter",x);
		}
		else if (noctiles <= 4) {  /* nearest full */
			flmoon(nlast,2,&fljd);
			x = jd - fljd;
			if(x < 0.) 
			  oprntf("%3.1f days until full moon",(-1.*x));
			else
			  oprntf("%3.1f days after full moon",x);
		}
		else if (noctiles <= 6) {  /* nearest last quarter */
			flmoon(nlast,3,&lqjd);
			x = jd - lqjd;
			if(x < 0.)
			  oprntf("%3.1f days before last quarter",(-1.*x));
			else
			  oprntf("%3.1f days after last quarter",x);
		}
		else oprntf("%3.1f days before new moon",(newjd - jd));
	}
}       

double lunskybright(alpha,rho,kzen,altmoon,alt, moondist) 

	double alpha,rho,kzen,altmoon,alt,moondist; 
		   
/* Evaluates predicted LUNAR part of sky brightness, in 
   V magnitudes per square arcsecond, following K. Krisciunas
   and B. E. Schaeffer (1991) PASP 103, 1033.

   alpha = separation of sun and moon as seen from earth,
   converted internally to its supplement,
   rho = separation of moon and object,
   kzen = zenith extinction coefficient, 
   altmoon = altitude of moon above horizon,
   alt = altitude of object above horizon 
   moondist = distance to moon, in earth radii

   all are in decimal degrees. */

{

    double istar,Xzm,Xo,Z,Zmoon,Bmoon,fofrho,rho_rad,test;

    rho_rad = rho/DEG_IN_RADIAN;
    alpha = (180. - alpha); 
    Zmoon = (90. - altmoon)/DEG_IN_RADIAN;
    Z = (90. - alt)/DEG_IN_RADIAN;
    moondist = moondist/(60.27);  /* divide by mean distance */

    istar = -0.4*(3.84 + 0.026*fabs(alpha) + 4.0e-9*pow(alpha,4.)); /*eqn 20*/
    istar =  pow(10.,istar)/(moondist * moondist);
    if(fabs(alpha) < 7.)   /* crude accounting for opposition effect */
	istar = istar * (1.35 - 0.05 * fabs(istar));
	/* 35 per cent brighter at full, effect tapering linearly to 
	   zero at 7 degrees away from full. mentioned peripherally in 
	   Krisciunas and Scheafer, p. 1035. */
    fofrho = 229087. * (1.06 + cos(rho_rad)*cos(rho_rad));
    if(fabs(rho) > 10.)
       fofrho=fofrho+pow(10.,(6.15 - rho/40.));            /* eqn 21 */
    else if (fabs(rho) > 0.25)
       fofrho= fofrho+ 6.2e7 / (rho*rho);   /* eqn 19 */
    else fofrho = fofrho+9.9e8;  /*for 1/4 degree -- radius of moon! */
    Xzm = sqrt(1.0 - 0.96*sin(Zmoon)*sin(Zmoon));
    if(Xzm != 0.) Xzm = 1./Xzm;  
	  else Xzm = 10000.;     
    Xo = sqrt(1.0 - 0.96*sin(Z)*sin(Z));
    if(Xo != 0.) Xo = 1./Xo;
	  else Xo = 10000.; 
    Bmoon = fofrho * istar * pow(10.,(-0.4*kzen*Xzm)) 
	  * (1. - pow(10.,(-0.4*kzen*Xo)));   /* nanoLamberts */
    if(Bmoon > 0.001) 
      return(22.50 - 1.08574 * log(Bmoon/34.08)); /* V mag per sq arcs-eqn 1 */
    else return(99.);                                     
}

void accusun(jd,lst,geolat,ra,dec,dist,topora,topodec,x,y,z)  

	double jd,lst,geolat,*ra,*dec,*dist,*topora,*topodec;
 	double *x, *y, *z;  
{
      /*  implemenataion of Jean Meeus' more accurate solar
	  ephemeris.  For ultimate use in helio correction! From
	  Astronomical Formulae for Calculators, pp. 79 ff.  This
	  gives sun's position wrt *mean* equinox of date, not
	  *apparent*.  Accuracy is << 1 arcmin.  Positions given are
	  geocentric ... parallax due to observer's position on earth is 
	  ignored. This is up to 8 arcsec; routine is usually a little 
	  better than that. 
          // -- topocentric correction *is* included now. -- //
	  Light travel time is apparently taken into
	  account for the ra and dec, but I don't know if aberration is
	  and I don't know if distance is simlarly antedated. 

	  x, y, and z are heliocentric equatorial coordinates of the
	  EARTH, referred to mean equator and equinox of date. */

	double L, T, Tsq, Tcb;
	double M, e, Cent, nu, sunlong;
	double Lrad, Mrad, nurad, R;
	double A, B, C, D, E, H;
	double xtop, ytop, ztop, topodist, l, m, n, xgeo, ygeo, zgeo;

	jd = jd + etcorr(jd)/SEC_IN_DAY;  /* might as well do it right .... */
	T = (jd - 2415020.) / 36525.;  /* 1900 --- this is an oldish theory*/
	Tsq = T*T;
	Tcb = T*Tsq;
	L = 279.69668 + 36000.76892*T + 0.0003025*Tsq;
	M = 358.47583 + 35999.04975*T - 0.000150*Tsq - 0.0000033*Tcb;
	e = 0.01675104 - 0.0000418*T - 0.000000126*Tsq;

	L = circulo(L);
	M = circulo(M);
/*      printf("raw L, M: %15.8f, %15.8f\n",L,M); */

	A = 153.23 + 22518.7541 * T;  /* A, B due to Venus */
	B = 216.57 + 45037.5082 * T;
	C = 312.69 + 32964.3577 * T;  /* C due to Jupiter */
		/* D -- rough correction from earth-moon 
			barycenter to center of earth. */
	D = 350.74 + 445267.1142*T - 0.00144*Tsq;  
	E = 231.19 + 20.20*T;    /* "inequality of long period .. */
	H = 353.40 + 65928.7155*T;  /* Jupiter. */
	
	A = circulo(A) / DEG_IN_RADIAN;
	B = circulo(B) / DEG_IN_RADIAN;
	C = circulo(C) / DEG_IN_RADIAN;
	D = circulo(D) / DEG_IN_RADIAN;
	E = circulo(E) / DEG_IN_RADIAN;
	H = circulo(H) / DEG_IN_RADIAN;

	L = L + 0.00134 * cos(A) 
	      + 0.00154 * cos(B)
	      + 0.00200 * cos(C)
	      + 0.00179 * sin(D)
	      + 0.00178 * sin(E);

	Lrad = L/DEG_IN_RADIAN;
	Mrad = M/DEG_IN_RADIAN;
	
	Cent = (1.919460 - 0.004789*T -0.000014*Tsq)*sin(Mrad)
	     + (0.020094 - 0.000100*T) * sin(2.0*Mrad)
	     + 0.000293 * sin(3.0*Mrad);
	sunlong = L + Cent;


	nu = M + Cent;
	nurad = nu / DEG_IN_RADIAN;
	
	R = (1.0000002 * (1 - e*e)) / (1. + e * cos(nurad));
	R = R + 0.00000543 * sin(A)
	      + 0.00001575 * sin(B)
	      + 0.00001627 * sin(C)
	      + 0.00003076 * cos(D)
	      + 0.00000927 * sin(H);
/*      printf("solar longitude: %10.5f  Radius vector %10.7f\n",sunlong,R);
	printf("eccentricity %10.7f  eqn of center %10.5f\n",e,Cent);   */
	
	sunlong = sunlong/DEG_IN_RADIAN;

	*dist = R;
	*x = cos(sunlong);  /* geocentric */
	*y = sin(sunlong);
	*z = 0.;
	eclrot(jd, x, y, z);
	
/*      --- code to include topocentric correction for sun .... */

	geocent(lst,geolat,0.,&xgeo,&ygeo,&zgeo);
	
	xtop = *x - xgeo*EQUAT_RAD/ASTRO_UNIT;
	ytop = *y - ygeo*EQUAT_RAD/ASTRO_UNIT;
	ztop = *z - zgeo*EQUAT_RAD/ASTRO_UNIT;

	topodist = sqrt(xtop*xtop + ytop*ytop + ztop*ztop);
	
	l = xtop / (topodist);
	m = ytop / (topodist);
	n = ztop / (topodist);

	*topora = atan_circ(l,m) * HRS_IN_RADIAN;
	*topodec = asin(n) * DEG_IN_RADIAN; 

	*ra = atan_circ(*x,*y) * HRS_IN_RADIAN;
	*dec = asin(*z) * DEG_IN_RADIAN; 
	
	*x = *x * R * -1;  /* heliocentric */
	*y = *y * R * -1;
	*z = *z * R * -1;

}

double jd_moon_alt(alt,jdguess,lat,longit,elevsea)
                  
	double alt,jdguess,lat,longit,elevsea;
                  
{
	/* returns jd at which moon is at a given 
	altitude, given jdguess as a starting point. In current version
	uses high-precision moon -- execution time does not seem to be
	excessive on modern hardware.  If it's a problem on your machine,
	you can replace calls to 'accumoon' with 'lpmoon' and remove
	the 'elevsea' argument. */

	double jdout;
	double deriv, err, del = 0.002;
	double ra,dec,dist,geora,geodec,geodist,sid,ha,alt2,alt3,az,par;
	short i = 0;

	/* first guess */
	
	sid=lst(jdguess,longit);
	accumoon(jdguess,lat,sid,elevsea,&geora,&geodec,&geodist,
				&ra,&dec,&dist);
	ha = lst(jdguess,longit) - ra;
	alt2 = altit(dec,ha,lat,&az,&par);
	jdguess = jdguess + del;
	sid = lst(jdguess,longit);
	accumoon(jdguess,lat,sid,elevsea,&geora,&geodec,&geodist,
				&ra,&dec,&dist);
	alt3 = altit(dec,(sid - ra),lat,&az,&par);
	err = alt3 - alt;
	deriv = (alt3 - alt2) / del;
	while((fabs(err) > 0.1) && (i < 10)) {
		jdguess = jdguess - err/deriv;
		sid=lst(jdguess,longit);
		accumoon(jdguess,lat,sid,elevsea,&geora,&geodec,&geodist,
				&ra,&dec,&dist);
		alt3 = altit(dec,(sid - ra),lat,&az,&par);
		err = alt3 - alt;
		i++;
		if(i == 9) oprntf("Moonrise or -set calculation not converging!!...\n");
	}       
	if(i >= 9) jdguess = -1000.;
	jdout = jdguess;
	return(jdout);
}

double jd_sun_alt(alt,jdguess,lat,longit)
                 
	double alt,jdguess,lat,longit; 

{
	/* returns jd at which sun is at a given 
	altitude, given jdguess as a starting point. Uses
	low-precision sun, which is plenty good enough. */

	double jdout;
	double deriv, err, del = 0.002;
	double ra,dec,ha,alt2,alt3,az,par;
	short i = 0;

	/* first guess */
	
	lpsun(jdguess,&ra,&dec);
	ha = lst(jdguess,longit) - ra;
	alt2 = altit(dec,ha,lat,&az,&par);
	jdguess = jdguess + del;
	lpsun(jdguess,&ra,&dec);
	alt3 = altit(dec,(lst(jdguess,longit) - ra),lat,&az,&par);
	err = alt3 - alt;
	deriv = (alt3 - alt2) / del;
	while((fabs(err) > 0.1) && (i < 10)) {
		jdguess = jdguess - err/deriv;
		lpsun(jdguess,&ra,&dec);
		alt3 = altit(dec,(lst(jdguess,longit) - ra),lat,&az,&par);
		err = alt3 - alt;
		i++;
		if(i == 9) oprntf("Sunrise, set, or twilight calculation not converging!\n");
	}       
	if(i >= 9) jdguess = -1000.;
	jdout = jdguess;
	return(jdout);
}

float ztwilight(alt) 
	double alt;
{

/* evaluates a polynomial expansion for the approximate brightening
   in magnitudes of the zenith in twilight compared to its 
   value at full night, as function of altitude of the sun (in degrees).
   To get this expression I looked in Meinel, A.,
   & Meinel, M., "Sunsets, Twilight, & Evening Skies", Cambridge U.
   Press, 1983; there's a graph on p. 38 showing the decline of 
   zenith twilight.  I read points off this graph and fit them with a
   polynomial; I don't even know what band there data are for! */
/* Comparison with Ashburn, E. V. 1952, JGR, v.57, p.85 shows that this
   is a good fit to his B-band measurements.  */

	float y, val;
	
	y = (-1.* alt - 9.0) / 9.0;  /* my polynomial's argument...*/
	val = ((2.0635175 * y + 1.246602) * y - 9.4084495)*y + 6.132725;
	return(val);
}


void find_dst_bounds(int yr, double stdz, int use_dst,
			double *jdb, double *jde) 

{
	/* finds jd's at which daylight savings time begins 
	    and ends.  The parameter use_dst allows for a number
	    of conventions, namely:
		0 = don't use it at all (standard time all the time)
		1 = use USA convention (1st Sun in April to
		     last Sun in Oct 1986 - 2006; last Sun in April before;
                     2nd sunday in March to first in Nov from 2007 on.)
		2 = use Spanish convention (for Canary Islands)
		-1 = use Chilean convention (CTIO).
		-2 = Australian convention (for AAT).
	    Negative numbers denote sites in the southern hemisphere,
	    where jdb and jde are beginning and end of STANDARD time for
	    the year. 
	    It's assumed that the time changes at 2AM local time; so
	    when clock is set ahead, time jumps suddenly from 2 to 3,
	    and when time is set back, the hour from 1 to 2 AM local 
	    time is repeated.  This could be changed in code if need be. */

	struct date_time trial;
	int n_sundays = 0;

	if((use_dst == 1) || (use_dst == 0)) { 
	    /* USA Convention, and including no DST to be defensive */
	    /* Note that this ignores various wrinkles such as the
		brief Nixon administration flirtation with year-round DST,
		the extended DST of WW II, and so on. */
		trial.y = yr;
		trial.mo = 4;
		if(yr >= 1986) trial.d = 1;
		else trial.d = 30; 
		trial.h = 2;
		trial.mn = 0;
		trial.s = 0;

		/* Find first Sunday in April for 1986 to 2006 ... */
		if((yr >= 1986) && (yr <= 2006)) {
			while(day_of_week(date_to_jd(trial)) != 6) trial.d++;
		}

		/* Find second Sunday in MARCH for 2007+ (new energy bill) */
		else if(yr >= 2007) {
			trial.mo = 3;
			trial.d = 0;
			n_sundays = 0;
			while(n_sundays < 2) {
				trial.d++;
				if(day_of_week(date_to_jd(trial)) == 6) 
					n_sundays++;
			}
			
		}
			
		/* Find last Sunday in April for pre-1986 .... */
		else while(day_of_week(date_to_jd(trial)) != 6) 
				trial.d--;

		*jdb = date_to_jd(trial) + stdz/24.;    

		/* Find last Sunday in October ... for pre-2007*/
		if(yr <= 2006) { 
			trial.mo = 10;
			trial.d = 31;
			while(day_of_week(date_to_jd(trial)) != 6) trial.d--;
		}
		/* for post-2006, find first Sunday in November. */
		else {
			trial.mo = 11;
			trial.d = 1;
			while(day_of_week(date_to_jd(trial)) != 6) trial.d++;
		}

		*jde = date_to_jd(trial) + (stdz - 1.)/24.;             
	}
 
	else if (use_dst == 2) {  /* Spanish, for Canaries */
		trial.y = yr;
		trial.mo = 3;
		trial.d = 31; 
		trial.h = 2;
		trial.mn = 0;
		trial.s = 0;

		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d--;
		}
		*jdb = date_to_jd(trial) + stdz/24.;    
		trial.mo = 9;
		trial.d = 30;
		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d--;
		}
		*jde = date_to_jd(trial) + (stdz - 1.)/24.;             
	}              
	else if (use_dst == -1) {  /* Chilean, for CTIO, etc.  */
	   /* off daylight 2nd Sun in March, onto daylight 2nd Sun in October */
		trial.y = yr;
		trial.mo = 3;
		trial.d = 8;  /* earliest possible 2nd Sunday */
		trial.h = 2;
		trial.mn = 0;
		trial.s = 0;

		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d++;
		}
		*jdb = date_to_jd(trial) + (stdz - 1.)/24.;
			/* note jdb is beginning of STANDARD time in south,
				hence use stdz - 1. */  
		trial.mo = 10;
		trial.d = 8;
		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d++;
		}
		*jde = date_to_jd(trial) + stdz /24.;           
	}                       
	else if (use_dst == -2) {  /* For Anglo-Australian Telescope  */
	   /* off daylight 1st Sun in March, onto daylight last Sun in October */
		trial.y = yr;
		trial.mo = 3;
		trial.d = 1;  /* earliest possible 1st Sunday */
		trial.h = 2;
		trial.mn = 0;
		trial.s = 0;

		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d++;
		}
		*jdb = date_to_jd(trial) + (stdz - 1.)/24.;
			/* note jdb is beginning of STANDARD time in south,
				hence use stdz - 1. */  
		trial.mo = 10;
		trial.d = 31;
		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d--;
		}
		*jde = date_to_jd(trial) + stdz /24.;           
	}               
}

double zone(use_dst,stdz,jd,jdb,jde) 
           
	short use_dst;
	double stdz,jd,jdb,jde; 
           
{
	/* Returns zone time offset when standard time zone is stdz,
	   when daylight time begins (for the year) on jdb, and ends
	   (for the year) on jde.  This is parochial to the northern
	   hemisphere.  */
	/* Extension -- specifying a negative value of use_dst reverses
	   the logic for the Southern hemisphere; then DST is assumed for
	   the Southern hemisphere summer (which is the end and beginning
	   of the year. */

	if(use_dst == 0) return(stdz);
	else if((jd > jdb) && (jd < jde) && (use_dst > 0)) return(stdz-1.);
	   /* next line .. use_dst < 0 .. for Southern Hemisphere sites. */
	else if(((jd < jdb) || (jd > jde)) && (use_dst < 0)) return(stdz-1.);
	else return(stdz);
}

double true_jd(date, use_dst, enter_ut, night_date, stdz) 
       
/* takes the values in the date-time structure, the standard time
   zone (in hours west), the prevailing conventions for date and
   time entry, and returns the value of the true julian date. */

	struct date_time date;
	short use_dst, enter_ut, night_date;
	double stdz;
{
	double jd, jdb, jde, test;

	if(enter_ut == 0) {
           find_dst_bounds(date.y,stdz,use_dst,&jdb,&jde);
	   jd = date_to_jd(date);
	   if((night_date == 1)  && (date.h < 12)) jd = jd + 1.;
	   if(use_dst != 0)  {  /* check at time changes */
		test = jd + stdz/24. - jdb;
		if((test > 0.) && (test < 0.041666666))   {
			/* 0.0416 = 1 hr; nonexistent time */
			oprntf("Error in true_jd -- nonexistent input time during std->dst change.\n");
			oprntf("Specify as 1 hour later!\n");
			return(-1.); /* signal of nonexistent time */
		}
		test = jd + stdz/24. - jde;
		if((test > 0.) && (test < 0.041666666))   {
			oprntf("WARNING ... ambiguous input time during dst->std change!\n");
		}
	   }
	   jd = jd + zone(use_dst,stdz,(jd+stdz/24.),jdb,jde)/24.;
			/* effect should be to default to standard time. */
        }
	else jd = date_to_jd(date);
	
	return(jd);
}


void print_tz(jd,use,jdb,jde,zabr) 
             
	double jd;
	short use;
	double jdb,jde;
	char zabr; 
             
{
	/* prints correct time abbreviation, given zabr as the
	   single character abbreviation for the time zone,
	   "D" or "S" depending on daylight or standard (dst 
	    begins at jdb, ends at jde) and current jd. */
    
	oprntf(" %c",zabr);
	if((jd > jdb) && (jd < jde) && (use > 0)) oprntf("D");
	  else if(((jd < jdb) || (jd > jde)) && (use < 0)) oprntf("D");  
	  else oprntf("S");
	oprntf("T");
}

void xyz_cel(x,y,z,ra,dec) 
	double x,y,z;  /* cartesian coordinate triplet */
	double *ra, *dec;  /* corresponding right ascension and declination,
                returned in decimal hours and decimal degrees. */

/* A much cleaner rewrite of the original skycalc code for this,
   which was transcribed from a PL/I routine .... */

{

   double mod;    /* modulus */
   double xy;     /* component in xy plane */

   /* normalize explicitly and check for bad input */

   mod = sqrt(x*x + y*y + z*z);
   if(mod > 0.) {
	   x = x / mod; y = y / mod; z = z / mod;  
   }
   else {   /* this has never happened ... */
	printf("Bad data in xyz_cel .... zero modulus position vector.\n");
	*ra = 0.; *dec = 0.;
	return;
   }
   
   xy = sqrt(x*x + y*y);

   if(xy < 1.0e-11) {   /* practically on a pole -- limit is arbitrary ...  */
      *ra = 0.;   /* degenerate anyway */
      *dec = PI_OVER_2; 
      if(z < 0.) *dec *= -1.;
   }
   else { /* in a normal part of the sky ... */
      *dec = asin(z);
      *ra = atan_circ(x,y);
   }

   *ra *= HRS_IN_RADIAN;
   *dec *= DEG_IN_RADIAN;
}

aberrate(epoch, vec, from_std) 
	double epoch,  /* decimal year ...  */
	vec[];  /* celestial unit vector ...  */
        int from_std;  /* 1 = apply aberration, -1 = take aberration out. */

/* corrects celestial unit vector for aberration due to earth's motion. 
   Uses accurate sun position ... replace with crude one for more speed if
   needed. */

{
	double jd, jd1, jd2, g, lambda, L, Xdot, Ydot, Zdot;   /* page C24 */
	double ras, decs, dists, topora, topodec; /* throwaways */
	double x, y, z, x1, y1, z1, x2, y2, z2;
	double norm;
	double rarad, decrad;
	
	/* find heliocentric velocity of earth as a fraction of the speed of light ... */

	jd = J2000 + (epoch - 2000.) * 365.25;  
        jd1 = jd - EARTH_DIFF;
	jd2 = jd + EARTH_DIFF;
	
	accusun(jd1,0.,0.,&ras,&decs,&dists,&topora,&topodec,&x1,&y1,&z1);
	accusun(jd2,0.,0.,&ras,&decs,&dists,&topora,&topodec,&x2,&y2,&z2);
	accusun(jd,0.,0.,&ras,&decs,&dists,&topora,&topodec,&x,&y,&z);

	Xdot = KMS_AUDAY*(x2 - x1)/(2.*EARTH_DIFF * SPEED_OF_LIGHT);  /* numerical differentiation */
	Ydot = KMS_AUDAY*(y2 - y1)/(2.*EARTH_DIFF * SPEED_OF_LIGHT);  /* crude but accurate */
	Zdot = KMS_AUDAY*(z2 - z1)/(2.*EARTH_DIFF * SPEED_OF_LIGHT);

	/* approximate correction ... non-relativistic but very close.  */

	vec[1] += from_std * Xdot;
        vec[2] += from_std * Ydot;
	vec[3] += from_std * Zdot;

	norm = pow((vec[1] * vec[1] + vec[2] * vec[2] + vec[3] * vec[3]), 0.5);

	vec[1] = vec[1] / norm; vec[2] = vec[2] / norm; vec[3] = vec[3] / norm;

}

void nutation_params(date_epoch, del_psi, del_ep) 
	double date_epoch, *del_psi, *del_ep;

/* computes the nutation parameters delta psi and
   delta epsilon at julian epoch (in years) using approximate
   formulae given by Jean Meeus, Astronomical Formulae for
   Calculators, Willman-Bell, 1985, pp. 69-70. Accuracy
   appears to be a few hundredths of an arcsec or better
   and numerics have been checked against his example. 
   Nutation parameters are returned in radians. */

{

	double T, jd, L, Lprime, M, Mprime, Omega;

	jd = (date_epoch - 2000.) * 365.25 + J2000;
	T = (jd - 2415020.0) / 36525.;
	
	L = 279.6967 + (36000.7689  + 0.000303 * T) * T;
	Lprime = 270.4342 + (481267.8831 - 0.001133 * T ) * T;
	M = 358.4758 + (35999.0498 - 0.000150 * T) * T;
	Mprime = 296.1046 + (477198.8491 + 0.009192 * T ) * T;
	Omega = 259.1833 - (1934.1420 - 0.002078 * T) * T;

	L = L / DEG_IN_RADIAN;
	Lprime = Lprime / DEG_IN_RADIAN;
	M = M / DEG_IN_RADIAN;
	Mprime = Mprime / DEG_IN_RADIAN;
	Omega = Omega / DEG_IN_RADIAN;
        	
	
	*del_psi = -1. * (17.2327 + 0.01737 * T) * sin(Omega) 
	   - (1.2729 + 0.00013 * T) * sin(2. * L)
	   + 0.2088 * sin(2 * Omega)
	   - 0.2037 * sin(2 * Lprime)
	   + (0.1261 - 0.00031 * T) * sin(M)
	   + 0.0675 * sin(Mprime)
	   - (0.0497 - 0.00012 * T) * sin(2 * L + M)
	   - 0.0342 * sin(2 * Lprime - Omega)
	   - 0.0261 * sin(2 * Lprime + Mprime)
	   + 0.0214 * sin(2 * L - M)
	   - 0.0149 * sin(2 * L - 2 * Lprime + Mprime)
	   + 0.0124 * sin(2 * L - Omega) 
	   + 0.0114 * sin(2 * Lprime - Mprime);

	*del_ep = (9.2100 + 0.00091 * T) * cos(Omega)
	   + (0.5522 - 0.00029 * T) * cos(2 * L)
	   - 0.0904 * cos(2 * Omega)
	   + 0.0884 * cos(2. * Lprime)
	   + 0.0216 * cos(2 * L + M)
	   + 0.0183 * cos(2 * Lprime - Omega)
	   + 0.0113 * cos(2 * Lprime + Mprime)
	   - 0.0093 * cos(2 * L - M)
	   - 0.0066 * cos(2 * L - Omega);
	
/*	printf("del_psi = %f, del_ep = %f\n",*del_psi,*del_ep);
*/
	*del_psi = *del_psi / ARCSEC_IN_RADIAN;
	*del_ep  = *del_ep  / ARCSEC_IN_RADIAN;

/*	printf("del_psi = %f, del_ep = %f\n",*del_psi,*del_ep);
*/
}


void cooxform(rin, din, std_epoch, 
  date_epoch, rout, dout, just_precess, from_std)
	double rin, din;  /* input ra and dec */
	double std_epoch;
	double date_epoch;        
	double *rout, *dout;  /* output */
	int just_precess;  /* flag ... 1 does just precession, 0 
			includes aberration and nutation. */
	int from_std;    /* flag ... 1 --> from std to date,
				    -1 --> from date to std. */ 

   /* General routine for precession and apparent place. Either
      transforms from current epoch (given by jd) to a standard
      epoch or back again, depending on value of the switch 
      "from_std"; 1 transforms from standard to current, -1 goes
      the other way.  Optionally does apparent place including
      nutation and annual aberration
      (but neglecting diurnal aberration, parallax, proper motion,
      and GR deflection of light); switch for this is "just_precess",
      1 does only precession, 0 includes other aberration & nutation. */

   /* Precession uses a matrix procedures
      as outlined in Taff's Computational Spherical Astronomy book.
      This is the so-called 'rigorous' method which should give very
      accurate answers all over the sky over an interval of several
      centuries.  Naked eye accuracy holds to ancient times, too. 
      Precession constants used are the new IAU1976 -- the 'J2000'
      system. 
 
      Nutation is incorporated into matrix formalism by constructing an 
      approximate nutation matrix and taking a matrix product with 
      precession matrix.  

      Aberration is done by adding the vector velocity of the earth to 
      the velocity of the light ray .... not kosher relativistically,
      but empirically correct to a high order for the angle.  */

{

   /* all the 3-d stuff is declared as [4] 'cause I'm not using the
	 zeroth element. */
	
   double ti, tf, zeta, z, theta;  /* all as per  Taff */
   double cosz, coszeta, costheta, sinz, sinzeta, sintheta;  /* ftns */
   double p[4][4];
      /* elements of the rotation matrix */
   double n[4][4]; 
      /* elements of the nutation matrix */
   double r[4][4];
      /* their product */
   double t[4][4];  /* temporary matrix for inversion .... */
   double radian_ra, radian_dec;
   double del_psi, del_eps, eps;  /* nutation angles in radians */
   double orig[4];   /* original unit vector */
   double fin[4];   /* final unit vector */
   int i, j, k;
   double out;    /* for diagnostics */
 
 
   ti = (std_epoch - 2000.) / 100.;
   tf = (date_epoch  - 2000. - 100. * ti) / 100.;

   zeta = (2306.2181 + 1.39656 * ti + 0.000139 * ti * ti) * tf +
    (0.30188 - 0.000344 * ti) * tf * tf + 0.017998 * tf * tf * tf;
   z = zeta + (0.79280 + 0.000410 * ti) * tf * tf + 0.000205 * tf * tf * tf;
   theta = (2004.3109 - 0.8533 * ti - 0.000217 * ti * ti) * tf
     - (0.42665 + 0.000217 * ti) * tf * tf - 0.041833 * tf * tf * tf;
 
   /* convert to radians */

   zeta = zeta / ARCSEC_IN_RADIAN;
   z = z / ARCSEC_IN_RADIAN;
   theta = theta / ARCSEC_IN_RADIAN;
  
   /* compute the necessary trig functions for speed and simplicity */
 
   cosz = cos(z);
   coszeta = cos(zeta);
   costheta = cos(theta);
   sinz = sin(z);
   sinzeta = sin(zeta);
   sintheta = sin(theta);

   /* compute the elements of the precession matrix -- set up
      here as *from* standard epoch *to* input jd. */

   p[1][1] = coszeta * cosz * costheta - sinzeta * sinz;
   p[1][2] = -1. * sinzeta * cosz * costheta - coszeta * sinz;
   p[1][3] = -1. * cosz * sintheta;

   p[2][1] = coszeta * sinz * costheta + sinzeta * cosz;
   p[2][2] = -1. * sinzeta * sinz * costheta + coszeta * cosz;
   p[2][3] = -1. * sinz * sintheta;
 
   p[3][1] = coszeta * sintheta;
   p[3][2] = -1. * sinzeta * sintheta;
   p[3][3] = costheta;

   if(just_precess == XFORM_DOAPPAR) {  /* if apparent place called for */
	
   	/* do the same for the nutation matrix. */

   	nutation_params(date_epoch, &del_psi, &del_eps); 
   	eps = 0.409105;  /* rough obliquity of ecliptic in radians */	

   	n[1][1] = 1.; n[2][2] = 1.; n[3][3] = 1.;
  	n[1][2] = -1. * del_psi * cos(eps);
   	n[1][3] = -1. * del_psi * sin(eps);
   	n[2][1] = -1. * n[1][2];
   	n[2][3] = -1. * del_eps;
   	n[3][1] = -1. * n[1][3];
   	n[3][2] = -1. * n[2][3];
 
	/* form product of precession and nutation matrices ... */
   	for(i = 1; i <= 3; i++) {
		for(j = 1; j <= 3; j++) {
			r[i][j] = 0.;
			for(k = 1; k <= 3; k++) 
				r[i][j] += p[i][k] * n[k][j];
        	}
   	}
   }
   else {  /* if you're just precessing .... */
	for(i = 1; i <= 3; i++) {
		for(j = 1; j <=3; j++) 
			r[i][j] = p[i][j];  /* simply copy precession matrix */
	}
   }

   /*   Commented out code ...
	for checking against Almanac tables ... turns out we're
   	ok to a few parts in 10^7. 

   for(i = 1; i <= 3; i++) {
	for(j = 1; j<= 3; j++) {
		printf("%d %d ... ",i,j);
		if(i == j) out = (r[i][j] - 1.) * 1.0e8;
		else out = r[i][j] * 1.0e8;
		printf("%7.0f\n",out);
	}
   }
    ........ end of commented out code. */

   /* The inverse of a rotation matrix is its transpose ... */

   if(from_std == XFORM_TOSTDEP) {    /* if you're transforming back to std
					epoch, rather than forward from std */
	for(i = 1; i <= 3; i++) {
		for(j = 1; j <= 3; j++)
			t[i][j] = r[j][i];  /* store transpose ... */		                               
	}
	for(i = 1; i <= 3; i++) {
		for(j = 1; j <= 3; j++)
			r[i][j] = t[i][j];  /* replace original w/ transpose.*/
	}
   }

   /* finally, transform original coordinates */
 
   radian_ra = rin / HRS_IN_RADIAN;
   radian_dec = din / DEG_IN_RADIAN;

   orig[1] = cos(radian_dec) * cos(radian_ra);
   orig[2] = cos(radian_dec) * sin(radian_ra);
   orig[3] = sin(radian_dec);


   if(from_std == XFORM_TOSTDEP && just_precess == XFORM_DOAPPAR) 
     /* if you're transforming from jd to std epoch, and doing apparent place,
	first step is to de-aberrate while still in epoch of date ... */	
  	 aberrate(date_epoch, orig, from_std); 
 
	
   for(i = 1; i<=3; i++) {
	fin[i] = 0.;
	for(j = 1; j<=3; j++) {
		fin[i] += r[i][j] * orig[j];
	}
   }

   if(from_std == XFORM_FROMSTD && just_precess == XFORM_DOAPPAR) 
	/* if you're transforming from std epoch to jd,
         last step is to apply aberration correction once you're in 
         equinox of that jd. */
  	 aberrate(date_epoch, fin, from_std); 

   /* convert back to spherical polar coords */

   xyz_cel(fin[1], fin[2], fin[3], rout, dout);

   return;	
}

double near_hor_refr(app_alt, pressure) 
	double app_alt, pressure;
{

	/* Almanac 1992, p. B62 -- ignores temperature variation */
	/* formula for near horizon, function-ized for iteration ... */
	
	double r;
        r = pressure * 
		(0.1594 + 0.0196 * app_alt + 0.00002 * app_alt * app_alt) /
		(293. * (1. + 0.505 * app_alt + 0.0845 * app_alt * app_alt));
	return(r);
}

double refract_size(alt, elev) 
	double alt;   /* altitude in degrees */
	double elev;  /* meters */

{
 	/* Almanac for 1992, p. B 62.  Ignores variation in temperature
           and just assumes T = 20 celsius.  */
	
	double pressure, r, altrad, crit_alt, app_alt;
	int i;

	altrad = alt / DEG_IN_RADIAN;

	pressure = 1013. * exp(-1. * elev/8620.);
		/* exponential atmosphere
		at T = 20 C, g = 980 cm/s^2,  and mean molecular wgt 28.8 */
	
	if(alt > 89.9) return(0.); /* avoid blowing up */
	else if (alt >= 15.0) {
		r = 0.00452 * pressure / (293. * tan(altrad));
		return(r);
	}
	else {  /* here have to start worrying about distinction between
                apparent and true altitude ... a pity as true altitude is
		what is handed in ... */
		crit_alt = -1. * pressure * 0.1594 / 293.;
		/* that's the *true* altitude corresponding to an
			*apparent* altitude of zero ... forget it if
			    it's below this. */
		if (alt > crit_alt) { /* go ahead and get it ... */ 
			app_alt = alt;  /* initial ... */
			for(i = 1; i <= 3; i++) {
			/* tests show 3 iterations is good to < 0.5 arcmin
			   for objects below geom horizon just barely rising ..
			   further accuracy is spurious. */
		           r = near_hor_refr(app_alt, pressure); 
			   app_alt = alt + r;
			}
			r = near_hor_refr(app_alt, pressure);
			return(r);
		}
	
		else {
			return(-1.);  /* below horizon. */
		}
	}
}

refract_corr(ha , dec, lat, elev, size, sense)
	double *ha, *dec, *size, lat, elev;
	int sense;
	
/* if sense == 1 , applies refraction to a true ha and dec; if
   == -1, de-corrects already refracted coordinates. Uses elevation of
   observatory above sea level to estimate a mean atmospheric pressure. */

{
	double x,y,z, xpr, ypr, localdec, localha, alt, az, norm,
		sinlat, coslat;
	localdec = *dec / DEG_IN_RADIAN;
	localha = *ha / HRS_IN_RADIAN;
	lat = lat / DEG_IN_RADIAN;  /* thank heavens for pass-by-value */
	sinlat = sin(lat);  coslat = cos(lat);

	/* The calculation is done by computing xyz coordinates in the
           horizon system, adding to the vertical component, renormalizing
           back to a unit vector, rotating to the polar system, and
           transforming back to ha and dec .... a long way around the
           barn, but completely general. */

	x =  cos(localdec)*cos(localha)*coslat + sin(localdec)*sinlat;
		/* vertical component */
	y =  sin(localdec)*coslat - cos(localdec)*cos(localha)*sinlat; 
		/* due N comp. */
	z =  -1. * cos(localdec)*sin(localha); 
		/* due east comp. */

	*size = refract_size(DEG_IN_RADIAN * asin(x), 0.); 
		/* (this gives zero for zenith) */
		
	if(*size > 0.) {  /* guard against singular result at zenith */
     		norm = pow((y * y + z * z), 0.5);  /* in-ground component */
	 
		x = norm * tan(atan(x/norm) + sense * *size / DEG_IN_RADIAN);
		norm = pow((x*x + y*y + z*z),0.5);
		
		x = x / norm; y = y / norm; z = z / norm;		
	
		xpr = x * coslat - y * sinlat;
		ypr = x * sinlat + y * coslat;
	
		*dec = asin(ypr) * DEG_IN_RADIAN;
		*ha = -1. * atan2(z, xpr) * HRS_IN_RADIAN;
	}
}

void mass_precess() {
  
    double rorig = 1., dorig, orig_epoch, final_epoch, rf, df, 
            mura = 0., mudec = 0., dt;
    short do_proper = 0;

    printf("Mass precession.  The '=' command does precessions in a pinch, but\n");
    printf("the present command is convenient for doing several (all with same\n");
    printf("starting and ending epochs).  This routine does not affect parameters\n");
    printf("in the rest of the program.\n\n");
    printf("Type 1 if you need proper motions (They're a pain!), or 0:");
    scanf("%hd",&do_proper);
    if(do_proper == 1) {
	printf("\nA proper motion correction will be included -- it's a simple\n");
	printf("linear correction (adds mu * dt to coordinate).  Proper motion\n");
	printf("itself is NOT rigorously transformed (as in B1950->J2000, which\n");
	printf("involves a small change of inertial reference frame!)\n\n");
    }
    printf("Give epoch to precess from:");
    scanf("%lf",&orig_epoch);
    printf("Give epoch to precess to:");
    scanf("%lf",&final_epoch);
    dt = final_epoch - orig_epoch;
    while(rorig >= 0.) {
	printf("\nGive RA for %7.2f (h m s, -1 0 0 to exit):",orig_epoch);
        rorig = get_coord();  
        if(rorig < 0.) {
		printf("Exiting precession routine.  Type '?' if you want a menu.\n");
		return;
        }
        printf("Give dec (d m s):");
        dorig = get_coord();
	if(do_proper == 1) get_pm(dorig,&mura,&mudec);
        cooxform(rorig, dorig, orig_epoch, final_epoch, &rf, &df,
	     XFORM_JUSTPRE, XFORM_FROMSTD);	
	rf = rf + mura * dt / 3600.;
        df = df + mudec * dt / 3600.;
        oprntf("\n\n %7.2f : RA = ",orig_epoch);
        put_coords(rorig,4,0);
        oprntf(", dec = ");
        put_coords(dorig,3,1);
        oprntf("\n %7.2f : RA = ",final_epoch);
        put_coords(rf,4,0);
        oprntf(", dec = ");
        put_coords(df,3,1);
        if(do_proper == 1) oprntf("\n RA p.m. = %8.4f sec/yr, dec = %8.3f arcsec/yr",mura,mudec);
	oprntf("\n");
    }
}

void print_apparent(rain,decin,epochin,mura_sec,mudec,jd,lat,longit,elev) 
	double rain, decin, epochin, mura_sec, mudec, jd, lat, longit, elev;
{
	double ra, dec, ep, dt, sid, ha, refr, alt, az, parang;

	ep = 2000. + (jd - J2000) / 365.25;
	oprntf("\nStandard epoch: %6.1f : ",epochin);
	put_coords(rain,4,0);
	oprntf(",  ");
	put_coords(decin,3,1);
	oprntf("\nPrecessed Only: %6.1f : ",ep);
	cooxform(rain,decin,epochin,ep,&ra,&dec,XFORM_JUSTPRE,XFORM_FROMSTD);
	put_coords(ra,4,0);
	oprntf(",  ");
	put_coords(dec,3,1);
	if(mura_sec != 0. || mudec != 0.) {
		oprntf("\n\nAnnual proper motions: mu-alpha = %8.5f sec\n",mura_sec);
		oprntf("                       mu-delta = %8.4f arcsec",mudec);
		dt = ep - epochin;
		rain = rain + mura_sec * dt / 3600.;
		decin = decin + mudec * dt / 3600.;
		oprntf("\nW/ Proper Motn: %6.1f : ",ep);
		cooxform(rain,decin,epochin,ep,&ra,&dec,XFORM_JUSTPRE,XFORM_FROMSTD);
		put_coords(ra,4,0);
		oprntf(",  ");
		put_coords(dec,3,1);
	}
	else oprntf("\n(Proper motions are zero.)");
	oprntf("\n\nApparent place: %6.1f : ",ep);
	cooxform(rain,decin,epochin,ep,&ra,&dec,XFORM_DOAPPAR,XFORM_FROMSTD);
	put_coords(ra,4,0);
	oprntf(",  ");
	put_coords(dec,3,1);
	oprntf("\n(includes nutation and aberration)");

	/* refraction part ... */

	sid = lst(jd, longit);  
    	ha = adj_time(sid - ra);
	refract_corr(&ha , &dec, lat, elev, &refr, 1);
	if(refr > 0.) {
	  altit(dec,ha,lat,&az,&parang);
	  oprntf("\n\nApplying est. refraction %7.1f arcsec in P.A. %6.1f degr:\n",
		3600. * refr, parang);
	  ra = sid - ha;
	  while(ra < 0.) ra += 24.;
	  while(ra > 24.) ra -= 24.;
	  oprntf("Refraction cor: %6.1f : ",ep);
	  put_coords(ra,4,0);
          oprntf(",  ");
	  put_coords(dec,3,1);
	}
	else if (refr < 0.) 
		oprntf("\n\nObject below level horizon ... no refraction correction.");
	else if (refr == 0.) oprntf("\n\nObject near zenith ... no refraction correction.");
	oprntf("\n");
}

/* 
Constellation part, adapted from a code by 
Francois Ochsenbein [CDS], which in turn is an implementation of 
algorithms invented by Barry Rappaport and Nancy Roman, 
\bibitem[Roman(1987)]{1987PASP...99..695R} Roman, N.~G.\ 1987, 
\pasp, 99, 695.  The original code is available from CDS, at
ftp://cdsarc.u-strasbg.fr/cats/VI/42/ .  It is used here with
permission from F. Ochsenbein.
 
Most of this is one long lookup table designed to allow 
efficient finding of the correct constellation from the 
RA and dec.  Note that the constellation boundaries are
given for 1875 coordinates, so the coords need to be 
tranformed to 1875 internally before the search is done. */

#define ITEMS(a)  (sizeof(a)/sizeof(a[0])) 

typedef struct row { float ral, rau, del; char cst[4]; } ROW ;

/* The Constellation Boundaries as extracted from Cat.#6042,
   paper by N.G. Roman, 1987PASP...99..695R
*/
static ROW condata[] = {
   { 0, 360, 90, "" },			/* North Pole */
 {  0.0000,360.0000, 88.0000,"UMi"},
 {120.0000,217.5000, 86.5000,"UMi"},
 {315.0000,345.0000, 86.1667,"UMi"},
 {270.0000,315.0000, 86.0000,"UMi"},
 {  0.0000,120.0000, 85.0000,"Cep"},
 {137.5000,160.0000, 82.0000,"Cam"},
 {  0.0000, 75.0000, 80.0000,"Cep"},
 {160.0000,217.5000, 80.0000,"Cam"},
 {262.5000,270.0000, 80.0000,"UMi"},
 {302.5000,315.0000, 80.0000,"Dra"},
 {  0.0000, 52.6250, 77.0000,"Cep"},
 {172.5000,203.7500, 77.0000,"Cam"},
 {248.0000,262.5000, 75.0000,"UMi"},
 {302.5000,310.0000, 75.0000,"Cep"},
 {119.5000,137.5000, 73.5000,"Cam"},
 {137.5000,170.0000, 73.5000,"Dra"},
 {195.0000,248.0000, 70.0000,"UMi"},
 { 46.5000, 51.2500, 68.0000,"Cas"},
 {306.2500,310.0000, 67.0000,"Dra"},
 {170.0000,180.0000, 66.5000,"Dra"},
 {  0.0000,  5.0000, 66.0000,"Cep"},
 {210.0000,235.0000, 66.0000,"UMi"},
 {353.7500,360.0000, 66.0000,"Cep"},
 {180.0000,202.5000, 64.0000,"Dra"},
 {202.5000,216.2500, 63.0000,"Dra"},
 {347.5000,353.7500, 63.0000,"Cep"},
 { 91.5000,105.0000, 62.0000,"Cam"},
 {300.0000,306.2500, 61.5000,"Dra"},
 {308.0500,309.0000, 60.9167,"Cep"},
 {105.0000,119.5000, 60.0000,"Cam"},
 {119.5000,126.2500, 60.0000,"UMa"},
 {296.5000,300.0000, 59.5000,"Dra"},
 {300.0000,308.0500, 59.5000,"Cep"},
 {343.0000,347.5000, 59.0833,"Cep"},
 {  0.0000, 36.5000, 58.5000,"Cas"},
 {291.2500,296.5000, 58.0000,"Dra"},
 { 25.5000, 28.6250, 57.5000,"Cas"},
 { 36.5000, 46.5000, 57.0000,"Cas"},
 { 46.5000, 47.5000, 57.0000,"Cam"},
 {334.7500,343.0000, 56.2500,"Cep"},
 { 75.0000, 91.5000, 56.0000,"Cam"},
 {210.5000,216.2500, 55.5000,"UMa"},
 {216.2500,291.2500, 55.5000,"Dra"},
 { 47.5000, 50.0000, 55.0000,"Cam"},
 {332.0000,334.7500, 55.0000,"Cep"},
 {309.0000,329.5000, 54.8333,"Cep"},
 {  0.0000, 25.5000, 54.0000,"Cas"},
 { 91.5000, 97.5000, 54.0000,"Lyn"},
 {181.2500,202.5000, 53.0000,"UMa"},
 {228.7500,236.2500, 53.0000,"Dra"},
 {329.5000,332.0000, 52.7500,"Cep"},
 { 50.0000, 75.0000, 52.5000,"Cam"},
 {343.0000,350.0000, 52.5000,"Cas"},
 {236.2500,255.0000, 51.5000,"Dra"},
 { 30.6250, 37.7500, 50.5000,"Per"},
 {255.0000,273.5000, 50.5000,"Dra"},
 {  0.0000, 20.5000, 50.0000,"Cas"},
 { 20.5000, 25.0000, 50.0000,"Per"},
 { 97.5000,102.0000, 50.0000,"Lyn"},
 {350.0000,360.0000, 50.0000,"Cas"},
 {202.5000,210.5000, 48.5000,"UMa"},
 {  0.0000, 16.7500, 48.0000,"Cas"},
 {353.7500,360.0000, 48.0000,"Cas"},
 {272.6250,273.5000, 47.5000,"Her"},
 {273.5000,286.2500, 47.5000,"Dra"},
 {286.2500,287.5000, 47.5000,"Cyg"},
 { 25.0000, 30.6250, 47.0000,"Per"},
 {126.2500,137.5000, 47.0000,"UMa"},
 {  2.5000, 13.0000, 46.0000,"Cas"},
 {180.0000,181.2500, 45.0000,"UMa"},
 {102.0000,110.5000, 44.5000,"Lyn"},
 {328.6250,329.5000, 44.0000,"Cyg"},
 {328.1250,328.6250, 43.7500,"Cyg"},
 {287.5000,291.0000, 43.5000,"Cyg"},
 {137.5000,152.5000, 42.0000,"UMa"},
 {152.5000,161.7500, 40.0000,"UMa"},
 {231.5000,236.2500, 40.0000,"Boo"},
 {236.2500,245.0000, 40.0000,"Her"},
 {138.7500,143.7500, 39.7500,"Lyn"},
 {  0.0000, 37.7500, 36.7500,"And"},
 { 37.7500, 38.5000, 36.7500,"Per"},
 {290.3750,291.0000, 36.5000,"Lyr"},
 { 67.5000, 70.3750, 36.0000,"Per"},
 {326.0000,328.1250, 36.0000,"Cyg"},
 {328.1250,330.0000, 36.0000,"Lac"},
 { 98.0000,110.5000, 35.5000,"Aur"},
 {110.5000,116.2500, 35.5000,"Lyn"},
 {  0.0000, 30.0000, 35.0000,"And"},
 {330.0000,342.2500, 35.0000,"Lac"},
 {342.2500,343.0000, 34.5000,"Lac"},
 {343.0000,352.5000, 34.5000,"And"},
 { 38.5000, 40.7500, 34.0000,"Per"},
 {161.7500,165.0000, 34.0000,"UMa"},
 {180.0000,185.0000, 34.0000,"CVn"},
 {116.2500,138.7500, 33.5000,"Lyn"},
 {138.7500,148.2500, 33.5000,"LMi"},
 { 10.7500, 21.1250, 33.0000,"And"},
 {227.7500,231.5000, 33.0000,"Boo"},
 {352.5000,356.2500, 32.0833,"And"},
 {185.0000,198.7500, 32.0000,"CVn"},
 {356.2500,360.0000, 31.3333,"And"}, 
 {209.3750,210.5000, 30.7500,"CVn"},
 { 36.2500, 40.7500, 30.6667,"Tri"},
 { 40.7500, 67.5000, 30.6667,"Per"},
 { 67.5000, 71.2500, 30.0000,"Aur"},
 {272.6250,290.3750, 30.0000,"Lyr"},
 {165.0000,180.0000, 29.0000,"UMa"},
 {295.0000,313.7500, 29.0000,"Cyg"},
 { 71.2500, 88.2500, 28.5000,"Aur"},
 {148.2500,157.5000, 28.5000,"LMi"},
 {198.7500,209.3750, 28.5000,"CVn"},
 {  0.0000,  1.0000, 28.0000,"And"},
 { 21.1250, 25.0000, 28.0000,"Tri"},
 { 88.2500, 98.0000, 28.0000,"Aur"},
 {118.2500,120.0000, 28.0000,"Gem"},
 {313.7500,326.0000, 28.0000,"Cyg"},
 {288.8750,295.0000, 27.5000,"Cyg"},
 { 28.7500, 36.2500, 27.2500,"Tri"},
 {242.5000,245.0000, 27.0000,"CrB"},
 {226.2500,227.7500, 26.0000,"Boo"},
 {227.7500,242.5000, 26.0000,"CrB"},
 {275.5000,283.0000, 26.0000,"Lyr"},
 {161.2500,165.0000, 25.5000,"LMi"},
 {283.0000,288.8750, 25.5000,"Lyr"},
 { 25.0000, 28.7500, 25.0000,"Tri"},
 { 10.7500, 12.7500, 23.7500,"Psc"},
 {157.5000,161.2500, 23.5000,"LMi"},
 {318.7500,321.2500, 23.5000,"Vul"},
 { 85.5000, 88.2500, 22.8333,"Tau"},
 {  1.0000,  2.1250, 22.0000,"And"},
 {238.7500,240.5000, 22.0000,"Ser"},
 { 88.2500, 93.2500, 21.5000,"Gem"},
 {297.5000,303.7500, 21.2500,"Vul"},
 {283.0000,288.7500, 21.0833,"Vul"},
 {  2.1250, 12.7500, 21.0000,"And"},
 {303.7500,308.5000, 20.5000,"Vul"},
 {117.1250,118.2500, 20.0000,"Gem"},
 {308.5000,318.7500, 19.5000,"Vul"},
 {288.7500,297.5000, 19.1667,"Vul"},
 { 49.2500, 50.5000, 19.0000,"Ari"},
 {283.0000,285.0000, 18.5000,"Sge"},
 { 85.5000, 86.5000, 18.0000,"Ori"},
 { 93.2500, 94.6250, 17.5000,"Gem"},
 {285.0000,297.5000, 16.1667,"Sge"},
 { 74.5000, 80.0000, 16.0000,"Tau"},
 {238.7500,241.2500, 16.0000,"Her"},
 {297.5000,303.7500, 15.7500,"Sge"},
 { 69.2500, 74.5000, 15.5000,"Tau"},
 { 80.0000, 84.0000, 15.5000,"Tau"},
 {192.5000,202.5000, 15.0000,"Com"},
 {258.7500,273.7500, 14.3333,"Her"},
 {178.0000,192.5000, 14.0000,"Com"},
 {112.5000,117.1250, 13.5000,"Gem"},
 {251.2500,258.7500, 12.8333,"Her"},
 {  0.0000,  2.1250, 12.5000,"Peg"},
 { 84.0000, 86.5000, 12.5000,"Tau"},
 {105.0000,112.5000, 12.5000,"Gem"},
 {316.7500,320.0000, 12.5000,"Peg"},
 { 94.6250,104.0000, 12.0000,"Gem"},
 {273.7500,283.0000, 12.0000,"Her"},
 {313.1250,315.7500, 11.8333,"Del"},
 {315.7500,316.7500, 11.8333,"Peg"},
 {172.7500,178.0000, 11.0000,"Leo"},
 { 93.6250, 94.6250, 10.0000,"Ori"},
 {104.0000,105.0000, 10.0000,"Gem"},
 {117.1250,118.8750, 10.0000,"Cnc"},
 {357.5000,360.0000, 10.0000,"Peg"},
 { 25.0000, 49.2500,  9.9167,"Ari"},
 {302.1250,304.5000,  8.5000,"Del"},
 {202.5000,226.2500,  8.0000,"Boo"},
 {341.2500,357.5000,  7.5000,"Peg"},
 {118.8750,138.7500,  7.0000,"Cnc"},
 {138.7500,161.2500,  7.0000,"Leo"},
 {273.7500,279.9333,  6.2500,"Oph"},
 {279.9333,283.0000,  6.2500,"Aql"},
 {312.5000,313.1250,  6.0000,"Del"},
 {105.0000,105.2500,  5.5000,"CMi"},
 {273.7500,276.3750,  4.5000,"Ser"},
 {241.2500,251.2500,  4.0000,"Her"},
 {273.7500,276.3750,  3.0000,"Oph"},
 {322.0000,325.0000,  2.7500,"Peg"},
 {  0.0000, 30.0000,  2.0000,"Psc"},
 {278.7500,283.0000,  2.0000,"Ser"},
 {304.5000,312.5000,  2.0000,"Del"},
 {312.5000,320.0000,  2.0000,"Equ"},
 {320.0000,322.0000,  2.0000,"Peg"},
 {330.0000,341.2500,  2.0000,"Peg"},
 {325.0000,330.0000,  1.7500,"Peg"},
 {105.2500,108.0000,  1.5000,"CMi"},
 { 53.7500, 69.2500,  0.0000,"Tau"},
 { 69.2500, 70.0000,  0.0000,"Ori"},
 {108.0000,121.2500,  0.0000,"CMi"},
 {220.0000,226.2500,  0.0000,"Vir"},
 {267.5000,273.7500,  0.0000,"Oph"},
 { 39.7500, 49.2500, -1.7500,"Cet"},
 { 49.2500, 53.7500, -1.7500,"Tau"},
 {226.2500,244.0000, -3.2500,"Ser"},
 { 70.0000, 76.2500, -4.0000,"Ori"},
 { 87.5000, 93.6250, -4.0000,"Ori"},
 {267.5000,269.5000, -4.0000,"Ser"},
 {273.7500,278.7500, -4.0000,"Ser"},
 {278.7500,283.0000, -4.0000,"Aql"},
 {341.2500,357.5000, -4.0000,"Psc"},
 {161.2500,172.7500, -6.0000,"Leo"},
 {172.7500,177.5000, -6.0000,"Vir"},
 {  0.0000,  5.0000, -7.0000,"Psc"},
 {357.5000,360.0000, -7.0000,"Psc"},
 {213.7500,220.0000, -8.0000,"Vir"},
 {238.7500,244.0000, -8.0000,"Oph"},
 {300.0000,308.0000, -9.0000,"Aql"},
 {320.0000,328.0000, -9.0000,"Aqr"},
 {257.5000,269.5000,-10.0000,"Oph"},
 { 87.5000,121.2500,-11.0000,"Mon"},
 { 73.7500, 76.2500,-11.0000,"Eri"},
 { 76.2500, 87.5000,-11.0000,"Ori"},
 {121.2500,125.5000,-11.0000,"Hya"},
 {143.7500,161.2500,-11.0000,"Sex"},
 {177.5000,192.5000,-11.0000,"Vir"},
 {263.7500,265.0000,-11.6667,"Oph"},
 {283.0000,300.0000,-12.0333,"Aql"},
 { 72.5000, 73.7500,-14.5000,"Eri"},
 {308.0000,320.0000,-15.0000,"Aqr"},
 {257.5000,273.7500,-16.0000,"Ser"},
 {273.7500,283.0000,-16.0000,"Sct"},
 {125.5000,128.7500,-17.0000,"Hya"},
 {244.0000,245.6250,-18.2500,"Oph"},
 {128.7500,136.2500,-19.0000,"Hya"},
 {161.2500,162.5000,-19.0000,"Crt"},
 {244.0000,245.6250,-19.2500,"Sco"},
 {235.0000,238.7500,-20.0000,"Lib"},
 {188.7500,192.5000,-22.0000,"Crv"},
 {192.5000,213.7500,-22.0000,"Vir"},
 {136.2500,146.2500,-24.0000,"Hya"},
 { 25.0000, 39.7500,-24.3833,"Cet"},
 { 39.7500, 56.2500,-24.3833,"Eri"},
 {162.5000,177.5000,-24.5000,"Crt"},
 {177.5000,188.7500,-24.5000,"Crv"},
 {213.7500,223.7500,-24.5000,"Lib"},
 {244.0000,251.2500,-24.5833,"Oph"},
 {  0.0000, 25.0000,-25.5000,"Cet"},
 {320.0000,328.0000,-25.5000,"Cap"},
 {328.0000,357.5000,-25.5000,"Aqr"},
 {357.5000,360.0000,-25.5000,"Cet"},
 {146.2500,153.7500,-26.5000,"Hya"},
 { 70.5000, 72.5000,-27.2500,"Eri"},
 { 72.5000, 91.7500,-27.2500,"Lep"},
 {300.0000,320.0000,-28.0000,"Cap"},
 {153.7500,158.7500,-29.1667,"Hya"},
 {188.7500,223.7500,-29.5000,"Hya"},
 {223.7500,235.0000,-29.5000,"Lib"},
 {235.0000,240.0000,-29.5000,"Sco"},
 { 68.7500, 70.5000,-30.0000,"Eri"},
 {251.2500,264.0000,-30.0000,"Oph"},
 {264.0000,267.5000,-30.0000,"Sgr"},
 {158.7500,162.5000,-31.1667,"Hya"},
 { 91.7500,110.5000,-33.0000,"CMa"},
 {183.7500,188.7500,-33.0000,"Hya"},
 {162.5000,183.7500,-35.0000,"Hya"},
 { 52.5000, 56.2500,-36.0000,"For"},
 {125.5000,140.5000,-36.7500,"Pyx"},
 { 64.0000, 68.7500,-37.0000,"Eri"},
 {267.5000,287.5000,-37.0000,"Sgr"},
 {320.0000,345.0000,-37.0000,"PsA"},
 {345.0000,350.0000,-37.0000,"Scl"},
 { 45.0000, 52.5000,-39.5833,"For"},
 {140.5000,165.0000,-39.7500,"Ant"},
 {  0.0000, 25.0000,-40.0000,"Scl"},
 { 25.0000, 45.0000,-40.0000,"For"},
 { 58.0000, 64.0000,-40.0000,"Eri"},
 {350.0000,360.0000,-40.0000,"Scl"},
 {212.5000,223.7500,-42.0000,"Cen"},
 {235.0000,240.0000,-42.0000,"Lup"},
 {240.0000,246.3125,-42.0000,"Sco"},
 { 72.5000, 75.0000,-43.0000,"Cae"},
 { 75.0000, 98.7500,-43.0000,"Col"},
 {120.0000,125.5000,-43.0000,"Pup"},
 { 51.2500, 58.0000,-44.0000,"Eri"},
 {246.3125,267.5000,-45.5000,"Sco"},
 {267.5000,287.5000,-45.5000,"CrA"},
 {287.5000,305.0000,-45.5000,"Sgr"},
 {305.0000,320.0000,-45.5000,"Mic"},
 { 45.0000, 51.2500,-46.0000,"Eri"},
 { 67.5000, 72.5000,-46.5000,"Cae"},
 {230.0000,235.0000,-48.0000,"Lup"},
 {  0.0000, 35.0000,-48.1667,"Phe"},
 { 40.0000, 45.0000,-49.0000,"Eri"},
 { 61.2500, 64.0000,-49.0000,"Hor"},
 { 64.0000, 67.5000,-49.0000,"Cae"},
 {320.0000,330.0000,-50.0000,"Gru"},
 { 90.0000,120.0000,-50.7500,"Pup"},
 {120.0000,122.5000,-50.7500,"Vel"},
 { 36.2500, 40.0000,-51.0000,"Eri"},
 { 57.5000, 61.2500,-51.0000,"Hor"},
 {  0.0000, 27.5000,-51.5000,"Phe"},
 { 90.0000, 92.5000,-52.5000,"Car"},
 {122.5000,126.7500,-53.0000,"Vel"},
 { 52.5000, 57.5000,-53.1667,"Hor"},
 { 57.5000, 60.0000,-53.1667,"Dor"},
 {  0.0000, 23.7500,-53.5000,"Phe"},
 { 32.5000, 36.2500,-54.0000,"Eri"},
 { 67.5000, 75.0000,-54.0000,"Pic"},
 {225.7500,230.0000,-54.0000,"Lup"},
 {126.7500,132.5000,-54.5000,"Vel"},
 { 92.5000, 97.5000,-55.0000,"Car"},
 {177.5000,192.5000,-55.0000,"Cen"},
 {212.5000,225.7500,-55.0000,"Lup"},
 {225.7500,230.0000,-55.0000,"Nor"},
 { 60.0000, 65.0000,-56.5000,"Dor"},
 {132.5000,165.0000,-56.5000,"Vel"},
 {165.0000,168.7500,-56.5000,"Cen"},
 {262.5000,270.0000,-57.0000,"Ara"},
 {270.0000,305.0000,-57.0000,"Tel"},
 {330.0000,350.0000,-57.0000,"Gru"},
 { 48.0000, 52.5000,-57.5000,"Hor"},
 { 75.0000, 82.5000,-57.5000,"Pic"},
 { 97.5000,102.5000,-58.0000,"Car"},
 {  0.0000, 20.0000,-58.5000,"Phe"},
 { 20.0000, 32.5000,-58.5000,"Eri"},
 {350.0000,360.0000,-58.5000,"Phe"},
 { 65.0000, 68.7500,-59.0000,"Dor"},
 {230.0000,246.3125,-60.0000,"Nor"},
 {305.0000,320.0000,-60.0000,"Ind"},
 { 82.5000, 90.0000,-61.0000,"Pic"},
 {227.5000,230.0000,-61.0000,"Cir"},
 {246.3125,248.7500,-61.0000,"Ara"},
 {223.7500,227.5000,-63.5833,"Cir"},
 {248.7500,251.2500,-63.5833,"Ara"},
 { 90.0000,102.5000,-64.0000,"Pic"},
 {102.5000,135.5000,-64.0000,"Car"},
 {168.7500,177.5000,-64.0000,"Cen"},
 {177.5000,192.5000,-64.0000,"Cru"},
 {192.5000,218.0000,-64.0000,"Cen"},
 {202.5000,205.0000,-65.0000,"Cir"},
 {251.2500,252.5000,-65.0000,"Ara"},
 { 32.5000, 48.0000,-67.5000,"Hor"},
 { 48.0000, 68.7500,-67.5000,"Ret"},
 {221.2500,223.7500,-67.5000,"Cir"},
 {252.5000,262.5000,-67.5000,"Ara"},
 {262.5000,270.0000,-67.5000,"Pav"},
 {330.0000,350.0000,-67.5000,"Tuc"},
 { 68.7500, 98.7500,-70.0000,"Dor"},
 {205.0000,221.2500,-70.0000,"Cir"},
 {221.2500,255.0000,-70.0000,"TrA"},
 {  0.0000, 20.0000,-75.0000,"Tuc"},
 { 52.5000, 68.7500,-75.0000,"Hyi"},
 { 98.7500,135.5000,-75.0000,"Vol"},
 {135.5000,168.7500,-75.0000,"Car"},
 {168.7500,205.0000,-75.0000,"Mus"},
 {270.0000,320.0000,-75.0000,"Pav"},
 {320.0000,350.0000,-75.0000,"Ind"},
 {350.0000,360.0000,-75.0000,"Tuc"},
 { 11.2500, 20.0000,-76.0000,"Tuc"},
 {  0.0000, 52.5000,-82.5000,"Hyi"},
 {115.0000,205.0000,-82.5000,"Cha"},
 {205.0000,270.0000,-82.5000,"Aps"},
 { 52.5000,115.0000,-85.0000,"Men"},
 {  0.0000,360.0000,-90.0000,"Oct"},
} ;

void radec_to_constel(double ra, double dec, double epoch, char *constname)  {

/* takes a position in the sky (ra, dec, and epoch (sic))  and returns 
 * the name of the constellation in 'constname'. */
	
  ROW *pr, *pe, *br ;
  float ral, rau; 
  double ra1875, de1875 ;
  char *pgm, *p, *form, *z ;
  int n ;

  pe = condata + ITEMS(condata) ;   /* End-of-data sentinel */

  cooxform(ra,dec,1875.,epoch,&ra1875,&de1875,XFORM_JUSTPRE,XFORM_TOSTDEP);
  /* 
  printf("1875: ");
  put_coords(ra1875,3,0);
  printf(" ");
  put_coords(de1875,3,0);
  */
  ra1875 *= 15.;   /* convert to degrees */
  /*
  printf(" -> %f %f (deg)\n",ra1875,de1875);
  */
  for (pr=condata+1; pr<pe; pr++) {
     /* printf("%f %f %f %s\n",pr->ral,pr->rau,pr->del,pr->cst); */
     if ((ra1875 >= pr->ral) && (ra1875 < pr->rau) && 
			    (de1875 >= pr->del)) break ;
     }
     if (pr < pe) p = pr->cst;
     else p = "???" ;
     strcpy(constname,p);
}


void galact(ra,dec,epoch,glong,glat)
           
	double ra,dec,epoch,*glong,*glat;
           
{
	/* Homebrew algorithm for 3-d Euler rotation into galactic. 
	   Perfectly rigorous, and with reasonably accurate input 
	   numbers derived from original IAU definition of galactic
	   pole (12 49, +27.4, 1950) and zero of long (at PA 123 deg
	   from pole.) */

	double  p11= -0.066988739415,
		p12= -0.872755765853,
		p13= -0.483538914631,
		p21=  0.492728466047,
		p22= -0.450346958025,
		p23=  0.744584633299,
		p31= -0.867600811168,
		p32= -0.188374601707,
		p33=  0.460199784759;  /* derived from Euler angles of 
		theta   265.610844031 deg (rotates x axis to RA of galact center),
		phi     28.9167903483 deg (rotates x axis to point at galact cent),
		omega   58.2813466094 deg (rotates z axis to point at galact pole) */

	double r1950,d1950,
		x0,y0,z0,x1,y1,z1,
		check;

/*   EXCISED CODE .... creates matrix from Euler angles. Resurrect if
     necessary to create new Euler angles for better precision. 
     Program evolved by running and calculating angles,then initializing
     them to the values they will always have, thus saving space and time.

	cosphi = cos(phi); and so on
	p11 = cosphi * costhet;
	p12 = cosphi * sinthet;
	p13 = -1. * sinphi;
	p21 = sinom*sinphi*costhet - sinthet*cosom;
	p22 = cosom*costhet + sinthet*sinphi*sinom;
	p23 = sinom*cosphi;
	p31 = sinom*sinthet + cosom*sinphi*costhet;
	p32 = cosom*sinphi*sinthet - costhet*sinom;
	p33 = cosom * cosphi;

	printf("%15.10f %15.10f %15.10f\n",p11,p12,p13);
	printf("%15.10f %15.10f %15.10f\n",p21,p22,p23);
	printf("%15.10f %15.10f %15.10f\n",p31,p32,p33);

	check = p11*(p22*p33-p32*p23) - p12*(p21*p33-p31*p23) +
		p13*(p21*p32-p22*p31);   
	printf("Check: %lf\n",check);  check determinant .. ok 

    END OF EXCISED CODE..... */

	/* precess to 1950 */

	cooxform(ra,dec,1950.,epoch,&r1950,&d1950,
		XFORM_JUSTPRE,XFORM_TOSTDEP);   /* transform *to* std epoch ... */
	r1950 = r1950 / HRS_IN_RADIAN;
	d1950 = d1950 / DEG_IN_RADIAN;
	
	/* form direction cosines */
	x0 = cos(r1950) * cos(d1950);
	y0 = sin(r1950) * cos(d1950);
	z0 = sin(d1950);

	/* rotate 'em */
	x1 = p11*x0 + p12*y0 + p13*z0;
	y1 = p21*x0 + p22*y0 + p23*z0;
	z1 = p31*x0 + p32*y0 + p33*z0;
	
	/* translate to spherical polars for Galactic coords. */
	*glong = atan_circ(x1,y1)*DEG_IN_RADIAN;
	*glat = asin(z1)*DEG_IN_RADIAN;
}  

void eclipt(ra,dec,epoch,jd,curep,eclong,eclat)

           double ra,dec,epoch,jd,*curep,*eclong,*eclat;
           
 /* ra in decimal hrs, other coords in dec. deg. */

/* converts ra and dec to ecliptic coords -- precesses to current
   epoch first (and hands current epoch back for printing.)  */

{
	double incl, T;
	double racur, decur;
	double x0,y0,z0,x1,y1,z1;
	
	T = (jd - J2000)/36525.;  /* centuries since J2000 */
	*curep = 2000. + (jd - J2000) / 365.25; 

	incl = (23.439291 + T * (-0.0130042 - 0.00000016 * T))/DEG_IN_RADIAN; 
		/* 1992 Astron Almanac, p. B18, dropping the 
		   cubic term, which is 2 milli-arcsec! */

	cooxform(ra,dec,epoch,*curep,&racur,&decur,
		XFORM_JUSTPRE,XFORM_FROMSTD);
	racur = racur / HRS_IN_RADIAN;
	decur = decur / DEG_IN_RADIAN;

	x0=cos(decur)*cos(racur);
	y0=cos(decur)*sin(racur);
	z0=sin(decur);

	x1=x0;  /* doesn't change */
	y1 = cos(incl)*y0 + sin(incl)*z0;
	z1 = -1 * sin(incl)*y0 + cos(incl)*z0;
	*eclong = atan_circ(x1,y1) * DEG_IN_RADIAN;
	*eclat = asin(z1) * DEG_IN_RADIAN;
}               


/* Planetary part, added 1992 August.  The intention of this is
   to compute low-precision planetary positions for general info
   and to inform user if observation might be interfered with by
   a planet -- a rarity, but it happens.  Also designed to make
   handy low-precision planet positions available for casual
   planning purposes.  Do not try to point blindly right at the
   middle of a planetary disk with these routines!  */

double jd_el;     /* ************** */

/* elements of planetary orbits */
struct elements {
	char name[9];
	double incl;
	double Omega;
	double omega;
	double a;
	double daily;
	double ecc;
	double L_0;
	double mass;
};

struct elements el[10];

void comp_el(jd) 

	double jd;
{

   double T, Tsq, Tcb, d;
   double ups, P, Q, S, V, W, G, H, zeta, psi; /* Meeus p. 110 ff. */
   double sinQ,sinZeta,cosQ,cosZeta,sinV,cosV,
	sin2Zeta,cos2Zeta;

   jd_el = jd;   /* true, but not necessarily; set explicitly */
   d = jd - 2415020.;   
   T = d / 36525.;
   Tsq = T * T;
   Tcb = Tsq * T;

/* computes and loads mean elements for planets.  */

/* Mercury, Venus, and Mars from Explanatory Suppl., p. 113 */

   strcpy(el[1].name,"Mercury");
   el[1].incl = 7.002880 + 1.8608e-3 * T - 1.83e-5 * Tsq;
   el[1].Omega = 47.14594 + 1.185208 * T + 1.74e-4 * Tsq;
   el[1].omega = 75.899697 + 1.55549 * T + 2.95e-4 * Tsq;
   el[1].a = .3870986;
   el[1].daily = 4.0923388;
   el[1].ecc = 0.20561421 + 0.00002046 * T;
   el[1].L_0 = 178.179078 + 4.0923770233 * d  +
	 0.0000226 * pow((3.6525 * T),2.);

   strcpy(el[2].name,"Venus  ");
   el[2].incl = 3.39363 + 1.00583e-03 * T - 9.722e-7 * Tsq;
   el[2].Omega = 75.7796472 + 0.89985 * T + 4.1e-4 * Tsq;
   el[2].omega = 130.16383 + 1.4080 * T + 9.764e-4 * Tsq;
   el[2].a = .723325;
   el[2].daily = 1.60213049;
   el[2].ecc = 0.00682069 - 0.00004774 * T;
   el[2].L_0 = 342.767053 + 1.6021687039 * 36525 * T +
	 0.000023212 * pow((3.6525 * T),2.);

/* Earth from old Nautical Almanac .... */

   strcpy(el[3].name,"Earth  ");
   el[3].ecc = 0.01675104 - 0.00004180*T + 0.000000126*Tsq;
   el[3].incl = 0.0;
   el[3].Omega = 0.0;
   el[3].omega = 101.22083 + 0.0000470684*d + 0.000453*Tsq + 0.000003*Tcb;
   el[3].a = 1.0000007;;
   el[3].daily = 0.985599;
   el[3].L_0 = 358.47583 + 0.9856002670*d - 0.000150*Tsq - 0.000003*Tcb +
	    el[3].omega;

   strcpy(el[4].name,"Mars   ");
   el[4].incl = 1.85033 - 6.75e-04 * T - 1.833e-5 * Tsq;
   el[4].Omega = 48.786442 + .770992 * T + 1.39e-6 * Tsq;
   el[4].omega = 334.218203 + 1.840758 * T + 1.299e-4 * Tsq;
   el[4].a = 1.5236915;
   el[4].daily = 0.5240329502 + 1.285e-9 * T;
   el[4].ecc = 0.09331290 - 0.000092064 * T - 0.000000077 * Tsq;
   el[4].L_0 = 293.747628 + 0.5240711638 * d  +
	 0.000023287 * pow((3.6525 * T),2.);

/* Outer planets from Jean Meeus, Astronomical Formulae for
   Calculators, 3rd edition, Willman-Bell; p. 100. */
   
   strcpy(el[5].name,"Jupiter");
   el[5].incl = 1.308736 - 0.0056961 * T + 0.0000039 * Tsq;
   el[5].Omega = 99.443414 + 1.0105300 * T + 0.0003522 * Tsq 
		- 0.00000851 * Tcb;
   el[5].omega = 12.720972 + 1.6099617 * T + 1.05627e-3 * Tsq
	- 3.43e-6 * Tcb;
   el[5].a = 5.202561;
   el[5].daily = 0.08312941782;
   el[5].ecc = .04833475  + 1.64180e-4 * T - 4.676e-7*Tsq -
	1.7e-9 * Tcb;
   el[5].L_0 = 238.049257 + 3036.301986 * T + 0.0003347 * Tsq -
	1.65e-6 * Tcb;

   /* The outer planets have such large mutual interactions that
      even fair accuracy requires lots of perturbations --- here
      are some of the larger ones, from Meeus' book. */

   ups = 0.2*T + 0.1;
   P = (237.47555 + 3034.9061 * T) / DEG_IN_RADIAN;
   Q = (265.91650 + 1222.1139 * T) / DEG_IN_RADIAN;
   S = (243.51721 + 428.4677 * T) / DEG_IN_RADIAN;
   V = 5*Q - 2*P;
   W = 2*P - 6*Q + 3*S;
   zeta = Q - P;
   psi = S - Q;
   sinQ = sin(Q);
   cosQ = cos(Q);
   sinV = sin(V);
   cosV = cos(V);
   sinZeta = sin(zeta);
   cosZeta = cos(zeta);
   sin2Zeta = sin(2*zeta);
   cos2Zeta = cos(2*zeta);

   el[5].L_0 = el[5].L_0 
	+ (0.331364 - 0.010281*ups - 0.004692*ups*ups)*sinV
	+ (0.003228 - 0.064436*ups + 0.002075*ups*ups)*cosV
	- (0.003083 + 0.000275*ups - 0.000489*ups*ups)*sin(2*V)
	+ 0.002472 * sin(W) + 0.013619 * sinZeta + 0.018472 * sin2Zeta 
	+ 0.006717 * sin(3*zeta) 
	+ (0.007275  - 0.001253*ups) * sinZeta * sinQ
	+ 0.006417 * sin2Zeta * sinQ  
	- (0.033839 + 0.001253 * ups) * cosZeta * sinQ 
	- (0.035681 + 0.001208 * ups) * sinZeta * sinQ;
	/* only part of the terms, the ones first on the list and
	   selected larger-amplitude terms from farther down. */

   el[5].ecc = el[5].ecc + 1e-7 * (
	  (3606 + 130 * ups - 43 * ups*ups) * sinV 
	+ (1289 - 580 * ups) * cosV - 6764 * sinZeta * sinQ 
	- 1110 * sin2Zeta * sin(Q) 
	+ (1284 + 116 * ups) * cosZeta * sinQ 
	+ (1460 + 130 * ups) * sinZeta * cosQ 
	+ 6074 * cosZeta * cosQ);

   el[5].omega = el[5].omega 
	+ (0.007192 - 0.003147 * ups) * sinV
	+ ( 0.000197*ups*ups - 0.00675*ups - 0.020428) * cosV
	+ 0.034036 * cosZeta * sinQ + 0.037761 * sinZeta * cosQ;

   el[5].a = el[5].a + 1.0e-6 * (
	205 * cosZeta - 263 * cosV + 693 * cos2Zeta + 312 * sin(3*zeta)
	+ 147 * cos(4*zeta) + 299 * sinZeta * sinQ 
	+ 181 * cos2Zeta * sinQ + 181 * cos2Zeta * sinQ
	+ 204 * sin2Zeta * cosQ + 111 * sin(3*zeta) * cosQ 
	- 337 * cosZeta * cosQ - 111 * cos2Zeta * cosQ
	);
   
   strcpy(el[6].name,"Saturn ");
   el[6].incl = 2.492519 - 0.00034550*T - 7.28e-7*Tsq;
   el[6].Omega = 112.790414 + 0.8731951*T - 0.00015218*Tsq - 5.31e-6*Tcb ;
   el[6].omega = 91.098214 + 1.9584158*T + 8.2636e-4*Tsq;
   el[6].a = 9.554747;
   el[6].daily = 0.0334978749897;
   el[6].ecc = 0.05589232 - 3.4550e-4 * T - 7.28e-7*Tsq;
   el[6].L_0 = 266.564377 + 1223.509884*T + 0.0003245*Tsq - 5.8e-6*Tcb 
	+ (0.018150*ups - 0.814181 + 0.016714 * ups*ups) * sinV 
	+ (0.160906*ups - 0.010497 - 0.004100 * ups*ups) * cosV
	+ 0.007581 * sin(2*V) - 0.007986 * sin(W) 
	- 0.148811 * sinZeta - 0.040786*sin2Zeta 
	- 0.015208 * sin(3*zeta) - 0.006339 * sin(4*zeta) 
	- 0.006244 * sinQ
	+ (0.008931 + 0.002728 * ups) * sinZeta * sinQ 
	- 0.016500 * sin2Zeta * sinQ
	- 0.005775 * sin(3*zeta) * sinQ 
	+ (0.081344 + 0.003206 * ups) * cosZeta * sinQ 
	+ 0.015019 * cos2Zeta * sinQ 
	+ (0.085581 + 0.002494 * ups) * sinZeta * cosQ
	+ (0.025328 - 0.003117 * ups) * cosZeta * cosQ
	+ 0.014394 * cos2Zeta * cosQ;   /* truncated here -- no
		      terms larger than 0.01 degrees, but errors may
		      accumulate beyond this.... */
   el[6].ecc = el[6].ecc + 1.0e-7 * (
	  (2458 * ups - 7927.) * sinV + (13381. + 1226. * ups) * cosV
	+ 12415. * sinQ + 26599. * cosZeta * sinQ 
	- 4687. * cos2Zeta * sinQ - 12696. * sinZeta * cosQ 
	- 4200. * sin2Zeta * cosQ +(2211. - 286*ups) * sinZeta*sin(2*Q)
	- 2208. * sin2Zeta * sin(2*Q) 
	- 2780. * cosZeta * sin(2*Q) + 2022. * cos2Zeta*sin(2*Q) 
	- 2842. * sinZeta * cos(2*Q) - 1594. * cosZeta * cos(2*Q)
	+ 2162. * cos2Zeta*cos(2*Q) );  /* terms with amplitudes
	    > 2000e-7;  some secular variation ignored. */
   el[6].omega = el[6].omega 
	+ (0.077108 + 0.007186 * ups - 0.001533 * ups*ups) * sinV
	+ (0.045803 - 0.014766 * ups - 0.000536 * ups*ups) * cosV
	- 0.075825 * sinZeta * sinQ - 0.024839 * sin2Zeta*sinQ
	- 0.072582 * cosQ - 0.150383 * cosZeta * cosQ +
	0.026897 * cos2Zeta * cosQ;  /* all terms with amplitudes 
	    greater than 0.02 degrees -- lots of others! */
   el[6].a = el[6].a + 1.0e-6 * (
	2933. * cosV + 33629. * cosZeta - 3081. * cos2Zeta 
	- 1423. * cos(3*zeta) + 1098. * sinQ - 2812. * sinZeta * sinQ 
	+ 2138. * cosZeta * sinQ  + 2206. * sinZeta * cosQ 
	- 1590. * sin2Zeta*cosQ + 2885. * cosZeta * cosQ 
	+ 2172. * cos2Zeta * cosQ);  /* terms with amplitudes greater
	   than 1000 x 1e-6 */

   strcpy(el[7].name,"Uranus ");
   el[7].incl = 0.772464 + 0.0006253*T + 0.0000395*Tsq;
   el[7].Omega = 73.477111 + 0.4986678*T + 0.0013117*Tsq;
   el[7].omega = 171.548692 + 1.4844328*T + 2.37e-4*Tsq - 6.1e-7*Tcb;
   el[7].a = 19.21814;
   el[7].daily = 1.1769022484e-2;
   el[7].ecc = 0.0463444 - 2.658e-5 * T;
   el[7].L_0 = 244.197470 + 429.863546*T + 0.000316*Tsq - 6e-7*Tcb;
   /* stick in a little bit of perturbation -- this one really gets
      yanked around.... after Meeus p. 116*/
   G = (83.76922 + 218.4901 * T)/DEG_IN_RADIAN;
   H = 2*G - S;
   el[7].L_0 = el[7].L_0 + (0.864319 - 0.001583 * ups) * sin(H)
	+ (0.082222 - 0.006833 * ups) * cos(H)
	+ 0.036017 * sin(2*H);
   el[7].omega = el[7].omega + 0.120303 * sin(H) 
	+ (0.019472 - 0.000947 * ups) * cos(H)
	+ 0.006197 * sin(2*H);
   el[7].ecc = el[7].ecc + 1.0e-7 * (
	20981. * cos(H) - 3349. * sin(H) + 1311. * cos(2*H));
   el[7].a = el[7].a - 0.003825 * cos(H);

   /* other corrections to "true longitude" are ignored. */

   strcpy(el[8].name,"Neptune");
   el[8].incl = 1.779242 - 9.5436e-3 * T - 9.1e-6*Tsq;
   el[8].Omega = 130.681389 + 1.0989350 * T + 2.4987e-4*Tsq - 4.718e-6*Tcb;
   el[8].omega = 46.727364 + 1.4245744*T + 3.9082e-3*Tsq - 6.05e-7*Tcb;
   el[8].a = 30.10957;
   el[8].daily = 6.020148227e-3;
   el[8].ecc = 0.00899704 + 6.33e-6 * T;
   el[8].L_0 = 84.457994 + 219.885914*T + 0.0003205*Tsq - 6e-7*Tcb;
   el[8].L_0 = el[8].L_0 
	- (0.589833 - 0.001089 * ups) * sin(H) 
	- (0.056094 - 0.004658 * ups) * cos(H)
	- 0.024286 * sin(2*H);
   el[8].omega = el[8].omega + 0.024039 * sin(H) 
	- 0.025303 * cos(H);
   el[8].ecc = el[8].ecc + 1.0e-7 * (
	4389. * sin(H) + 1129. * sin(2.*H)
	+ 4262. * cos(H) + 1089. * cos(2.*H));
   el[8].a = el[8].a + 8.189e-3 * cos(H); 

/* crummy -- osculating elements a la Sept 15 1992 */

   d = jd - 2448880.5;  /* 1992 Sep 15 */       
   T = d / 36525.;
   strcpy(el[9].name,"Pluto  ");
   el[9].incl = 17.1426;
   el[9].Omega = 110.180;
   el[9].omega = 223.782;
   el[9].a = 39.7465;
   el[9].daily = 0.00393329;
   el[9].ecc = 0.253834;
   el[9].L_0 = 228.1027 + 0.00393329 * d;
/*   printf("inc Om om : %f %f %f\n",el[9].incl,el[9].Omega,el[9].omega);
   printf("a  dail ecc: %f %f %f\n",el[9].a,el[9].daily,el[9].ecc);
   printf("L_0 %f\n",el[9].L_0);
*/
   el[1].mass = 1.660137e-7;  /* in units of sun's mass, IAU 1976 */
   el[2].mass = 2.447840e-6;  /* from 1992 *Astron Almanac*, p. K7 */
   el[3].mass = 3.040433e-6;  /* earth + moon */
   el[4].mass = 3.227149e-7;
   el[5].mass = 9.547907e-4;
   el[6].mass = 2.858776e-4; 
   el[7].mass = 4.355401e-5;  
   el[8].mass = 5.177591e-5;
   el[9].mass = 7.69e-9;  /* Pluto+Charon -- ? */

}
	
void planetxyz(p, jd, x, y, z)
             
	int p; 
	double jd, *x, *y, *z;
              
/* produces ecliptic x,y,z coordinates for planet number 'p'
   at date jd. */

{
	double M, omnotil, nu, r;
	double e, LL, Om, om, nuu, ii;

/* see 1992 Astronomical Almanac, p. E 4 for these formulae. */

	ii = el[p].incl/DEG_IN_RADIAN;
	e = el[p].ecc;

	LL = (el[p].daily * (jd - jd_el) + el[p].L_0) / DEG_IN_RADIAN;
	Om = el[p].Omega / DEG_IN_RADIAN;
	om = el[p].omega / DEG_IN_RADIAN;
	
	M = LL - om;
	omnotil = om - Om;
	nu = M + (2.*e - 0.25 * pow(e,3.)) * sin(M) +
	     1.25 * e * e * sin(2 * M) +
	     1.08333333 * pow(e,3.) * sin(3 * M);
	r = el[p].a * (1. - e*e) / (1 + e * cos(nu));

	*x = r * 
	     (cos(nu + omnotil) * cos(Om) - sin(nu +  omnotil) * 
		cos(ii) * sin(Om));
	*y = r * 
	     (cos(nu +  omnotil) * sin(Om) + sin(nu +  omnotil) * 
		cos(ii) * cos(Om));
	*z = r * sin(nu +  omnotil) * sin(ii);
}


void planetvel(p, jd, vx, vy, vz)
	int p; 
	double jd, *vx, *vy, *vz;
{
/* numerically evaluates planet velocity by brute-force
numerical differentiation. Very unsophisticated algorithm. */

	double dt; /* timestep */
	double x1,y1,z1,x2,y2,z2,r1,d1,r2,d2,ep1;

	dt = 0.1 / el[p].daily; /* time for mean motion of 0.1 degree */
	planetxyz(p, (jd - dt), &x1, &y1, &z1);
	planetxyz(p, (jd + dt), &x2, &y2, &z2);
	*vx = 0.5 * (x2 - x1) / dt;    
	*vy = 0.5 * (y2 - y1) / dt;
	*vz = 0.5 * (z2 - z1) / dt;   
	/* answer should be in ecliptic coordinates, in AU per day.*/
}

void xyz2000(jd,x,y,z) 
	double jd, x, y, z;

/* simply transforms a vector x, y, and z to 2000 coordinates
   and prints -- for use in diagnostics. */
{ 
	double r1, d1, ep1, r2, d2, mod; 

	mod = sqrt(x*x + y*y + z*z);
	xyz_cel(x,y,z,&r1,&d1);
	ep1 = 2000. + (jd - J2000)/365.25;
	cooxform(r1,d1,2000.,ep1,&r2,&d2,XFORM_JUSTPRE,XFORM_TOSTDEP);
	x = mod * cos(r2/HRS_IN_RADIAN) * cos(d2/DEG_IN_RADIAN);
	y = mod * sin(r2/HRS_IN_RADIAN) * cos(d2/DEG_IN_RADIAN);
	z = mod * sin(d2/DEG_IN_RADIAN);
	printf("%f to 2000 -->  %f %f %f \n",ep1,x,y,z);

}

void xyz2000xf(jd,x,y,z) 
	double jd, *x, *y, *z;

/* simply transforms a vector x, y, and z to 2000 coordinates
   and hands it back through pointers.  For really transforming it! */
{ 
	double r1, d1, ep1, r2, d2, mod; 

	mod = sqrt(*x * *x + *y * *y + *z * *z);
	xyz_cel(*x,*y,*z,&r1,&d1);
	ep1 = 2000. + (jd - J2000)/365.25;
	cooxform(r1,d1,2000.,ep1,&r2,&d2,XFORM_JUSTPRE,XFORM_TOSTDEP);
	*x = mod * cos(r2/HRS_IN_RADIAN) * cos(d2/DEG_IN_RADIAN);
	*y = mod * sin(r2/HRS_IN_RADIAN) * cos(d2/DEG_IN_RADIAN);
	*z = mod * sin(d2/DEG_IN_RADIAN);
}

void earthview(x, y, z, i, ra, dec) 
 	
	double *x, *y, *z;
	int i; 
	double *ra, *dec;

/* given computed planet positions for planets 1-10, computes
   ra and dec of i-th planet as viewed from earth (3rd) */

{
	double dx, dy, dz;

	dx = x[i] - x[3];
	dy = y[i] - y[3];
	dz = z[i] - z[3];


	xyz_cel(dx,dy,dz,ra,dec);

}

void pposns(jd,lat,sid,print_option,planra,plandec) 

	double jd,lat,sid;
	short print_option;
	double *planra, *plandec; 
                      
/* computes and optionally prints positions for all the planets. */
/*  print_option 1 = print positions, 0 = silent */

{
	int i;
	double x[10], y[10], z[10], ha, alt, az, secz, par;
	double rasun, decsun, distsun, topora,topodec;
	double georamoon,geodecmoon,geodistmoon,toporamoon,topodecmoon,
	      topodistmoon; 

	accusun(jd,0.,0.,&rasun,&decsun,&distsun,&topora,&topodec,x+3,y+3,z+3);
/*      planetxyz(3,jd,x+3,y+3,z+3);   get the earth first (EarthFirst!?) 
	eclrot(jd,x+3,y+3,z+3);  */

	accumoon(jd,lat,sid,0.,&georamoon,&geodecmoon,&geodistmoon,
			 &toporamoon,&topodecmoon,&topodistmoon);         

	if(print_option == 1) {  /* set up table header */
		oprntf("\n\nPlanetary positions (epoch of date), accuracy about 0.1 deg:\n");
		oprntf("\n             RA       dec       HA");
		oprntf("       sec.z     alt   az\n\n");

	/* Throw in the sun and moon here ... */
		oprntf("Sun    : ");
		put_coords(topora,1,0);
		oprntf("   ");
		put_coords(topodec,0,1);
		ha = adj_time(sid - topora);
		oprntf("   ");
		put_coords(ha,0,1);
		alt = altit(topodec,ha,lat,&az,&par);
		secz = secant_z(alt);
 		if(fabs(secz) < 100.) oprntf("   %8.2f  ",secz);
		else oprntf("  (near horiz)");
		oprntf(" %5.1f  %5.1f\n",alt,az);
		oprntf("Moon   : ",el[i].name);
		put_coords(toporamoon,1,0);
		oprntf("   ");
		put_coords(topodecmoon,0,1);
		ha = adj_time(sid - toporamoon);
		oprntf("   ");
		put_coords(ha,0,1);
		alt=altit(topodecmoon,ha,lat,&az,&par);
		secz = secant_z(alt);
 		if(fabs(secz) < 100.) oprntf("   %8.2f  ",secz);
		else oprntf("  (near horiz)");
		oprntf(" %5.1f  %5.1f\n\n",alt,az);
        }
	for(i = 1; i <= 9; i++) {
		if(i == 3) goto SKIP;  /* skip the earth */
		planetxyz(i,jd,x+i,y+i,z+i);
		eclrot(jd,x+i,y+i,z+i);
		earthview(x,y,z,i,planra+i,plandec+i);
		if(print_option == 1) {
			oprntf("%s: ",el[i].name);
			put_coords(planra[i],1,0);
			oprntf("   ");
			put_coords(plandec[i],0,1);
			ha = adj_time(sid - planra[i]);
			oprntf("   ");
			put_coords(ha,0,1);
			alt=altit(plandec[i],ha,lat,&az,&par);
			secz = secant_z(alt);
 			if(fabs(secz) < 100.) oprntf("   %8.2f  ",secz);
			else oprntf("  (near horiz)");
			oprntf(" %5.1f  %5.1f",alt,az);
			if(i == 9) oprntf(" <-(least accurate)\n");
			else oprntf("\n");
		}
		SKIP: ;
	}
	if(print_option == 1) printf("Type command, or ? for menu:");
}

void barycor(jd,x,y,z,xdot,ydot,zdot)

	double jd,*x,*y,*z;
	double *xdot,*ydot,*zdot;
            
/* This routine takes the position
   x,y,z and velocity xdot,ydot,zdot, assumed heliocentric,
   and corrects them to the solar system barycenter taking into
   account the nine major planets.  Routine evolved by inserting
   planetary data (given above) into an earlier, very crude
   barycentric correction.  */

{

	int p;
	double xp, yp, zp, xvp, yvp, zvp;
	double xo, yo, zo;  /* for diagn */
	
	double xc=0.,yc=0.,zc=0.,xvc=0.,yvc=0.,zvc=0.;

	comp_el(jd);

	for(p=1;p<=9;p++) { /* sum contributions of the planets */
		planetxyz(p,jd,&xp,&yp,&zp);
		xc = xc + el[p].mass * xp;  /* mass is fraction of solar mass */
		yc = yc + el[p].mass * yp;
		zc = zc + el[p].mass * zc;
		planetvel(p,jd,&xvp,&yvp,&zvp);
		xvc = xvc + el[p].mass * xvp;
		yvc = yvc + el[p].mass * yvp;
		zvc = zvc + el[p].mass * zvc;
	/* diagnostic commented out ..... nice place to check planets if needed
		printf("%d :",p);
		xo = xp;
		yo = yp;
		zo = zp;
		eclrot(jd,&xo,&yo,&zo);
		xyz2000(jd,xo,yo,zo);
		printf("    ");
		xo = xvp;
		yo = yvp;
		zo = zvp;
		eclrot(jd,&xo,&yo,&zo);
		xyz2000(jd,xo,yo,zo);    */
	}
	/* normalize properly and rotate corrections to equatorial coords */
	xc = xc / SS_MASS;
	yc = yc / SS_MASS;
	zc = zc / SS_MASS;     /* might as well do it right ... */
	
	eclrot(jd, &xc, &yc, &zc);

/*      printf("posn corrn:");   diagnostic commented out
	xyz2000(jd,xc,yc,zc);
	printf(" vel corrn:");
	xyz2000(jd,(1.0e9 * xvc),(1.0e9 * yvc),(1.0e9 * zvc));  */

	xvc = xvc * KMS_AUDAY / SS_MASS;  
	yvc = yvc * KMS_AUDAY / SS_MASS;
	zvc = zvc * KMS_AUDAY / SS_MASS;
	eclrot(jd, &xvc, &yvc, &zvc);

	/* add them in */
	*x = *x - xc;  /* these are in AU -- */
	*y = *y - yc;
	*z = *z - zc;
/*      xyz2000(jd,*x,*y,*z);   */
	*xdot = *xdot - xvc;
	*ydot = *ydot - yvc;
	*zdot = *zdot - zvc;
	
/* diagnostic -- trash variables 
	-- scale for direct comparison with almanac
	xp = 1.0e9 * *xdot / KMS_AUDAY;    
	yp = 1.0e9 * *ydot / KMS_AUDAY;
	zp = 1.0e9 * *zdot / KMS_AUDAY;
	xyz2000(jd,xp,yp,zp);
					*/      
}

void helcor(jd,ra,dec,ha,lat,elevsea,tcor,vcor)

	double jd,ra,dec,ha;
  	double lat,elevsea,*tcor,*vcor;
           
/* finds heliocentric correction for given jd, ra, dec, ha, and lat.
   tcor is time correction in seconds, vcor velocity in km/s, to 
   be added to the observed values.  
   Input ra and dec assumed to be at current epoch */

{
	double M, ec, E, obliq, omeg, x, y, z, xdot, ydot, zdot;
	double xobj,yobj,zobj;
	double ras, decs, dists, jd1, jd2, x1, x2, y1, y2, z1, z2;
	double topora, topodec;
	double x_geo, y_geo, z_geo;  /* geocentric coords of observatory */
	double diff;
	double a=499.0047837;  /* light travel time for 1 AU, sec  */

	dec=dec/DEG_IN_RADIAN; /* pass by value! */
	ra = ra/HRS_IN_RADIAN;
	ha = ha/HRS_IN_RADIAN;
	
	xobj = cos(ra) * cos(dec);
	yobj = sin(ra) * cos(dec);
	zobj = sin(dec);
	
/* diagnostic -- temporarily trash jd1 
	jd1 = jd1 + etcorr(jd1) / SEC_IN_DAY;
	printf("TDT: %20f\n",jd1);  */

	jd1 = jd - EARTH_DIFF;
	jd2 = jd + EARTH_DIFF;
	
	accusun(jd1,0.,0.,&ras,&decs,&dists,&topora,&topodec,&x1,&y1,&z1);
	accusun(jd2,0.,0.,&ras,&decs,&dists,&topora,&topodec,&x2,&y2,&z2);
	accusun(jd,0.,0.,&ras,&decs,&dists,&topora,&topodec,&x,&y,&z);

/*      printf("ra dec distance:");  diagnostic -- commented out
	put_coords(ras,3,0);
	printf(" ");
	put_coords(decs,2,1);
	printf(" %f\n",dists);   */

	xdot = KMS_AUDAY*(x2 - x1)/(2.*EARTH_DIFF);  /* numerical differentiation */
	ydot = KMS_AUDAY*(y2 - y1)/(2.*EARTH_DIFF);  /* crude but accurate */
	zdot = KMS_AUDAY*(z2 - z1)/(2.*EARTH_DIFF);

/*      printf("Helio earth:");  diagnostic -- commmented out
	xyz2000(jd,x,y,z);
	xyz2000(jd,xdot,ydot,zdot);   */

	barycor(jd,&x,&y,&z,&xdot,&ydot,&zdot);
	*tcor = a * (x*xobj + y*yobj + z*zobj);
	*vcor = xdot * xobj + ydot * yobj + zdot * zobj;
/* correct diurnal rotation for elliptical earth including obs. elevation */
	geocent(0., lat, elevsea, &x_geo, &y_geo, &z_geo);
/* longitude set to zero arbitrarily so that x_geo = perp. distance to axis */
	*vcor = *vcor - 0.4651011 * x_geo * sin(ha) * cos(dec);
/* 0.4651011 = 6378.137 km radius * 2 pi / (86164.1 sec per sidereal day) */
/* could add time-of-flight across earth's radius here -- but rest of
   theory is not good to 0.02 seconds anyway. */
}

void lsrcor(ra,dec,epoch,vcor) 
	double ra,dec,epoch;
  	double *vcor;
{

/* Computes the correction from bary/helio centric to local standard
 of rest, i.e. (v wrt lsr) = (v helio) + vcor.  */
/* velocity of the sun is taken to be 13 km/s toward 1900 coords
   18 0 0 +30 0 0 */

	double ra1900, dec1900;
	double x, y, z;
	double xdotsun = 0., ydotsun = -11.258, zdotsun = 6.5;

	cooxform(ra,dec,1900.,epoch,&ra1900,&dec1900,
		XFORM_JUSTPRE,XFORM_TOSTDEP);   /* transform *to* 1900 */
	/*
	x = cos(dec1900/DEG_IN_RADIAN) * cos(ra1900/HRS_IN_RADIAN);
	 ... no x-component in std solar motion */
	y = cos(dec1900/DEG_IN_RADIAN) * sin(ra1900/HRS_IN_RADIAN);
	z = sin(dec1900/DEG_IN_RADIAN);

	*vcor = y * ydotsun + z * zdotsun;
}


void parellipse(double jd, double ra, double dec, double epoch, double lat, 
                double longit, double *dra, 
		double *ddec, double *aberra, double *aberdec)

/* Computes the normalized annual parallax shift at a given jd for 
   a star at ra (decimal hours) and dec (decimal degrees) (coords are for
   epoch).  Shifts are expressed as arcseconds in the ra direction and
   arcseconds in the dec direction; note that the dra is the actual 
   change of ra times cosine delta.  Uses skycalc routines to do a 
   vastly better job than needed.  Approximations break down within a 
   very small arc of the pole. 

   Observatory latitude and WEST longitude are lat (deg) and longit
    (hours).

   Also returned is the aberration in the ra and dec directions.

   Answers returned are on the equinox J2000.

*/
 
{
	double ra2000, dec2000;
	double epnow;
	double mod;
	double rarad, decrad;
	double rapar, decpar;
	double raaber, decaber;
	double xstar, ystar, zstar;
	double xpar, ypar, zpar;
	double xearth, yearth, zearth;
	double sid, ha; 
	double rasun, decsun, distsun;  /* not used */
	double toporasun, topodecsun;   /* also not used */
	double xdot, ydot, zdot;        /* also not used ... */
	double vec[4];  /* don't ask */
	
	epnow = 2000. + (jd - J2000) / 365.25;   /* julian epoch */
	
	cooxform(ra, dec, 2000., epnow, &ra2000, &dec2000, XFORM_JUSTPRE, 
			XFORM_TOSTDEP);
	rarad = ra2000 / HRS_IN_RADIAN;
	decrad = dec2000 / DEG_IN_RADIAN;
	
	xstar = cos(rarad) * cos(decrad);
	ystar = sin(rarad) * cos(decrad);
	zstar = sin(decrad);
	
	sid = lst(jd, longit);
	ha = sid - ra2000;   /* close enough for our purposes here */
	accusun(jd,sid,lat,&rasun,&decsun,&distsun,&toporasun,
			&topodecsun,&xearth,&yearth,&zearth);
	xdot = 1.; ydot = 1.; zdot = 1.;  /* to avoid fp error */
	barycor(jd,&xearth,&yearth,&zearth,&xdot,&ydot,&zdot);
	xyz2000xf(jd,&xearth,&yearth,&zearth);
	
	printf("Earth XYZ (AU, J2000 frame): %7.4f  %7.4f  %7.4f\n",xearth,yearth,zearth);
	
	xpar = xstar - xearth / ARCSEC_IN_RADIAN;
	ypar = ystar - yearth / ARCSEC_IN_RADIAN;
	zpar = zstar - zearth / ARCSEC_IN_RADIAN;
	
	xyz_cel(xpar,ypar,zpar,&rapar,&decpar);
	
	*dra = (rapar - ra2000) * cos(decrad) * 54000.;
	*ddec = (decpar - dec2000) * 3600.;
	
	printf("Parallax factors (arcsec, J2000): dRA = %7.4f, ddec = %7.4f\n",*dra,*ddec);
	
	vec[1] = xstar;  /* duuuh */
	vec[2] = ystar; 
	vec[3] = zstar;
	
	aberrate(epnow, vec, 1);
	
	xyz_cel(vec[1], vec[2], vec[3], &raaber, &decaber);
	*aberra  = (raaber - ra2000) * cos(decrad) * 54000.;
	*aberdec = (decaber - dec2000) * 3600.;
	
	printf("Annual aberration (arcsec): dra = %8.3f, ddec = %8.3f\n",*aberra,*aberdec);
	
}

/* A couple of eclipse predictors.... */

float overlap(r1, r2, sepn)

	double r1, r2, sepn; 
{
             
             
/* for two circles of radii r1 and r2,
   computes the overlap area in terms of the area of r1
   if their centers are separated by
   sepn. */

	float result, a, rlarge, rsmall;

	if(sepn > r1+r2) return(0.); /* They don't overlap */
	if((r2 < 0.) || (r1 < 0.)) return(-1.);  /* oops -- r's > 0.*/
	if(r2 > r1) {
		rlarge = r2;
		rsmall = r1;
	}
	else {
		rlarge = r1;
		rsmall = r2;
	}
	if(sepn < (rlarge - rsmall)) { /* small circle completely contained */
		if(r2 > r1) return(1.);  /* r1 is completely covered */
		else return((r2*r2) / (r1*r1)); /* ratio of areas */
	}
	a = (rlarge*rlarge + sepn*sepn - rsmall*rsmall) / (2.*sepn);
	   /* a is the distance along axis of centers to where circles cross */
	result = PI * rlarge * rlarge / 4. 
	     - a * sqrt(rlarge * rlarge - a*a) / 2.
	     - rlarge * rlarge * asin(a/rlarge) / 2;
	a = sepn - a;
	result = result + PI * rsmall * rsmall / 4.
	     - a * sqrt(rsmall * rsmall - a*a) / 2.
	     - rsmall * rsmall * asin(a/rsmall) / 2.;
	result = result * 2.;
	return(result / (PI*r1*r1)); /* normalize to circle 1's area. */  
}

void solecl(sun_moon,distmoon,distsun) 

	double sun_moon,distmoon,distsun; 
           
           
{
	double ang_sun, ang_moon; /* angular sizes */
	float magnitude; /* fraction of sun covered */

	ang_moon = DEG_IN_RADIAN * asin(RMOON / (distmoon * EQUAT_RAD));
	ang_sun = DEG_IN_RADIAN * asin(RSUN / (distsun * ASTRO_UNIT));

	if(sun_moon >= (ang_sun + ang_moon)) return ;

	else if(ang_sun >= ang_moon) {                   /* annular */
		magnitude = overlap(ang_sun,ang_moon,sun_moon);
		if(sun_moon > (ang_sun - ang_moon))
			oprntf("PARTIAL ECLIPSE OF THE SUN, %4.2f covered.\n",
				magnitude);
		else oprntf("ANNULAR ECLIPSE OF THE SUN, %4.2f covered.\n",
				magnitude);
	}
	else {
		magnitude = overlap(ang_sun,ang_moon,sun_moon);
		if(sun_moon > (ang_moon - ang_sun))
			oprntf("PARTIAL ECLIPSE OF THE SUN, %4.2f covered.\n",
				magnitude);
		else oprntf("TOTAL ECLIPSE OF THE SUN!\n");
	}
}

short lunecl(georamoon,geodecmoon,geodistmoon,rasun,decsun,distsun)

	double georamoon,geodecmoon,geodistmoon,rasun,decsun,distsun;
            
            
{
	/* quickie lunar eclipse predictor -- makes a number of 
	   minor assumptions, e. g. small angle approximations, plus
	   projects phenomena onto a plane at distance = geocentric 
	   distance of moon . */

	double ang_sun, ang_moon; /* angular sizes */
	double radius_um, radius_penum;
	double ra_shadow, dec_shadow;
	double lun_to_shadow;  /* angular separation of centerline of shadow
		  from center of moon ... */
	float magnitude;  /* portion covered */

	ang_sun = asin(RSUN / (distsun * ASTRO_UNIT));  /* radians */
	
	ra_shadow = rasun + 12.;
	if(ra_shadow > 24.) ra_shadow = ra_shadow - 24.; /* not strictly 
	   necessary, perhaps, but wise */
	dec_shadow = -1. * decsun;
	radius_um = (1./geodistmoon - ang_sun) * DEG_IN_RADIAN;  
	radius_penum = (1./geodistmoon + ang_sun) * DEG_IN_RADIAN;
	
		
	ang_moon = DEG_IN_RADIAN * asin(RMOON / (geodistmoon * EQUAT_RAD));
	lun_to_shadow = DEG_IN_RADIAN * 
		subtend(georamoon,geodecmoon,ra_shadow,dec_shadow);
	if(lun_to_shadow > (radius_penum + ang_moon)) return (0);
	else if(lun_to_shadow >= (radius_um + ang_moon))  {
		if(lun_to_shadow >= (radius_penum - ang_moon)) {
			oprntf("PARTIAL PENUMBRAL (BRIGHT) ECLIPSE OF THE MOON.\n");
			return (1);
		}
		else  {
			oprntf("PENUMBRAL (BRIGHT) ECLIPSE OF THE MOON.\n");
			return (2);
		}
	}
	else if (lun_to_shadow >= (radius_um - ang_moon)) {
		magnitude = overlap(ang_moon,radius_um,lun_to_shadow);
		oprntf("PARTIAL UMBRAL (DARK) ECLIPSE OF THE MOON, %4.2f covered.\n",
		   magnitude);
		return(3);
	}
	else {
		oprntf("TOTAL ECLIPSE OF THE MOON!\n");
		return(4);
	}
}
 
void planet_alert(jd,ra,dec,tolerance)
                 
	double jd,ra,dec,tolerance;
                 
/* given a jd, ra, and dec, this computes rough positions
   for all the planets, and alerts the user if any of them
   are within less than a settable tolerance of the ra and dec. */

{
	double pra[10],pdec[10], angle;
	int i;

	comp_el(jd);
	pposns(jd,0.,0.,0,pra,pdec);
	for(i = 1; i<=9 ; i++) {
		if(i == 3) goto SKIP;
		angle = subtend(pra[i],pdec[i],ra,dec) * DEG_IN_RADIAN;
		if(angle < tolerance) {
			oprntf("-- CAUTION -- proximity to %s -- low-precision calculation shows\n ",
				el[i].name);
			oprntf("this direction as %5.2f deg away from %s ---\n",
				angle,el[i].name);
		}       
		SKIP: ;
	}
}        

short setup_time_place(date,longit,lat,stdz,use_dst,zone_name,
        zabr,site_name,enter_ut,night_date,jdut,jdlocal,jdb,jde,sid,
	curepoch)
  
struct date_time date;
double lat, longit, stdz, *jdut, *jdlocal, *jdb, *jde, *sid, *curepoch;
short use_dst, enter_ut, night_date;
char zabr;
char *site_name;
char *zone_name;
       
/* This takes the date (which contains the time), and the site parameters,
   and prints out a banner giving the various dates and times; also
   computes and returns various jd's, the sidereal time, and the epoch. 
   Returns negative number to signal error if date is out of range of
   validity of algorithms, or if you specify a bad time during daylight-time
   change; returns zero if successful.  */

{
	
	double jd, jdloc, jdtest, yearday;

	if((date.y <= 1900) | (date.y >= 2100)) {
		oprntf("Date out of range - 1901 -> 2099\n");
		return(-1);
	}
	find_dst_bounds(date.y,stdz,use_dst,jdb,jde);
	oprntf("\nW Long (hms):");
	put_coords(longit,3,1);   
	oprntf(", lat (dms): ");
	put_coords(lat,2,1);
	oprntf(", std time zone %3.0f hr W\n",stdz);
	oprntf("\n");

        /* Establish the correct time from the input values ... 
             This is done more compactly in the true_jd function, but 
             the present code prints more warnings, etc, so leave it.*/

	if(enter_ut==0) {
		/* "night date" convention -- adjust date up by one if local
		      time is in morning. */
		if((night_date == 1) && (date.h < 12))
		      date.d = date.d + 1;  /* pass-by-value, I hope */
		jdloc = date_to_jd(date); /* local */
		jdtest = jdloc + stdz/24.; /* almost */
		if(use_dst != 0) {
		   if((fabs(jdtest-*jdb)<0.5) | (fabs(jdtest-*jde)<0.5)) {
		   oprntf("** be careful -- near standard and daylight change. ** \n");
		   }
		}
		if(use_dst > 0) {  /* handle exceptions during dst -- std change */
			if((jdtest - *jdb > 0.) 
				&& (jdtest - *jdb < 0.041666667)) {
				/* this time doesn't really exist */
			  oprntf("Nonexistent input time during change to DST!\n");
			  oprntf("Please respecify time as 1 hour later.\n");
			  return(-1);  /* return error flag. */
			}
			else if((jdtest - *jde >= 0.0) &&
				(jdtest - *jde < 0.04166667)) {
			  oprntf("\nCAUTION - DST --> STD -- ambiguous input time defaults to STD!\n"); 
			  oprntf("Use 'g' option to specify UT unambiguously.\n");
			  jd = jdloc + stdz / 24.;
			}
			else jd = jdloc + 
				zone(use_dst,stdz,jdtest,*jdb,*jde) / 24.;
		}
		else if(use_dst < 0) {  /* south -- reverse jdb & jde */
			if((jdtest - *jde > 0.) 
				&& (jdtest - *jde < 0.041666667)) {
				/* this time doesn't really exist */
			  oprntf("Nonexistent input time during change to DST!\n");
			  oprntf("Please respecify time as 1 hour later.\n");
			  return(-1);  /* return error flag. */
			}
			else if((jdtest - *jdb >= 0.0) &&
				(jdtest - *jdb < 0.04166667)) {
			  oprntf("\nCAUTION - DST --> STD -- ambiguous input time defaults to STD!\n"); 
			  oprntf("Use 'g' option to specify UT unambiguously.\n");
			  jd = jdloc + stdz / 24.;
			}
			else jd = jdloc + 
				zone(use_dst,stdz,jdtest,*jdb,*jde) / 24.;
		}
		else jd = jdtest;  /* not using dst ... no problem. */
		oprntf("Local Date and time: ");
		print_all(jdloc);
		oprntf(" ");
		print_tz(jd,use_dst,*jdb,*jde,zabr); 
/*              if(fabs(frac_part(jdloc) - 0.5) < MIDN_TOL)
		   oprntf("\n(Watch for possible day error near midnight.)");
 .... this really shouldn't happen now. */
		if(use_dst != 0) {
		}

		oprntf("\n   UT Date and time: ");
		print_all(jd);
		if(fabs(frac_part(jd) - 0.5) < 0.00001)
		   oprntf("\n(Watch for possible day error near midnight.)");
	}
	else {  /* times entered are ut already */

		jd = date_to_jd(date);
		oprntf("   UT Date and time: ");
		print_all(jd);
		if(fabs(frac_part(jd) - 0.5) < 0.00001)
		   oprntf("\n(Watch for possible day error near midnight.)");
		jdloc = jd - zone(use_dst,stdz,jd,*jdb,*jde)/24.; 
		if(use_dst != 0) {
		   if((fabs(jd-*jdb)<0.05) | (fabs(jd-*jde)<0.05))
			oprntf("\nCAUTION..dst/std time changing!");
		}
		oprntf("\nlocal Date and time: ");
		print_all(jdloc);
		oprntf(" ");
		print_tz(jd,use_dst,*jdb,*jde,zabr); 
		if(fabs(frac_part(jdloc) - 0.5) < 0.00001)
		   oprntf("\n(Watch for possible day error near midnight.)");
				
	}
	if(use_dst != 0) {
		if(stdz==zone(use_dst,stdz,jd,*jdb,*jde))
			oprntf("\nDST assumed not in effect.\n");
		else oprntf("\nDST assumed in effect.\n");
	}
	else printf("\n");
	oprntf("\nJulian date: %lf",jd);
	*sid=lst(jd,longit);    
	oprntf("   LMST: ");
	put_coords(*sid,3,0);
	*jdut = jd;
	*jdlocal = jdloc;
	*curepoch = 2000.+(jd-J2000)/365.25;
	return(0);
}

void print_tonight(date,lat,longit,elevsea,elev,horiz,site_name,stdz,
	zone_name,zabr,use_dst,jdb,jde,short_long)

struct date_time date;
double lat, longit, elevsea, elev, horiz, stdz, *jdb, *jde;
char *site_name, *zone_name, zabr;
short use_dst, short_long;     /* short_long is a fossil argument which 
                    allows a slightly shorter version to be printed. */

/* Given site and time information, prints a summary of
   the important phenomena for a single night. 

   The coding in this routine is extremely tortuous ... I even use
   the dreaded goto statement!  It's inelegant, but (a) the logic is
   indeed somewhat complicated and (b) it works.  */

{

	double jd, jdmid, jdcent, stmid, ramoon, decmoon, distmoon;
	double geora, geodec, geodist;  /* geocent for moon, not used here.*/
	double locjdb, locjde;
	double rasun, decsun, min_alt, max_alt;
	double hasunset, jdsunset, jdsunrise, sid;
	double hamoonrise, hamoonset, tmoonrise, tmoonset,
		jdmoonrise, jdmoonset;
	double hatwilight, jdtwilight, jdetw, jdmtw;
	double ill_frac;
	short dow; /* day of week */
	float set_to_rise, twi_to_twi, hrs,
		   moonless_hrs, moony_hrs, moon_print;
	
	find_dst_bounds(date.y,stdz,use_dst,jdb,jde);
	locjdb = *jdb-stdz/24.;
	locjde = *jde-(stdz-1)/24.;
	date.h = 18;  /* local afternoon */
	date.mn = 0;
	date.s = 0;  /* afternoon */
	
	oprntf("\nAlmanac for %s:\nlong. ",site_name);
	put_coords(longit,2,1);
	oprntf(" (h.m.s) W, lat. ");
	put_coords(lat,1,1);
	oprntf(" (d.m), elev. %5.0f m\n",elevsea);
	jd = date_to_jd(date); /* not really jd; local equivalent */
	if(use_dst > 0) {
		oprntf("%s Daylight Savings Time assumed from 2 AM on\n",zone_name);
		print_calendar(locjdb,&dow);
		oprntf(" to 2 AM on ");
		print_calendar(locjde,&dow);
		oprntf("; standard zone = %3.0f hrs W",stdz);
		if((fabs(jd - locjdb) < 0.45) || (fabs(jd - locjde) < 0.45)) 
			oprntf("\n   ** TIME CHANGE IS TONIGHT! **");
	}
	else if (use_dst < 0) {
		oprntf("%s Daylight Savings Time used before 2 AM \n",zone_name);
		print_calendar(*jdb-stdz/24.,&dow);
		oprntf(" and after 2 AM ");
		print_calendar(*jde-(stdz-1)/24.,&dow);
		oprntf("; standard zone = %3.0f hrs W",stdz);
		if((fabs(jd - locjdb) < 0.45) || (fabs(jd - locjde) < 0.45)) 
			oprntf("\n   ** TIME CHANGE IS TONIGHT! **");
	}
	else oprntf("%s Standard Time (%3.0f hrs W) in use all year.",zone_name,stdz);

	oprntf("\n\n");

	oprntf("For the night of: ");
	print_day(day_of_week(jd));
	oprntf(", ");
	print_calendar(jd,&dow); /* translate back e.g. 11/31-12/1 */
	jd = jd + .5;  /* local morning */
	oprntf(" ---> ");
	print_day(day_of_week(jd));
	oprntf(", ");
	print_calendar(jd,&dow);
	oprntf("\n");
	jd = jd - 0.25;  /* local midnight */
	jdmid = jd + zone(use_dst,stdz,jd,*jdb,*jde) / 24.;  
					/* corresponding ut */
	oprntf("Local midnight = ");
	print_calendar(jdmid,&dow);
	oprntf(", ");
	print_time(jdmid,-1); /* just the hours! */ 
	oprntf(" UT, or JD %11.3f\n",jdmid);    
	oprntf("Local Mean Sidereal Time at midnight = ");
	stmid = lst(jdmid,longit);
	put_coords(stmid,3,0);            
	oprntf("\n\n");
	
	accumoon(jdmid,lat,stmid,elevsea,
	   &geora,&geodec,&geodist,&ramoon,&decmoon,&distmoon);    
	lpsun(jdmid,&rasun,&decsun);
	hasunset = ha_alt(decsun,lat,-(0.83+horiz));
	if(hasunset > 900.) {  /* flag for never sets */
		oprntf("Sun up all night!\n");
		set_to_rise = 0.;
		twi_to_twi = 0.;
		jdcent = -1.;
		goto DO_MOON_ETC;  /* aargh! Ugly flow of control!! But
                   twilight certainly irrelevant if sun up all night, so
		   skip it. */
	}
	if(hasunset < -900.) {
		oprntf("Sun down all day!\n");
		set_to_rise = 24.;
		jdcent = -1.;
		goto CHECK_TWI18;  /* More ugly flow  - 
			checks for twilight even if sun down. */
	}               
	jdsunset = jdmid + adj_time(rasun+hasunset-stmid)/24.;  
		/* initial guess */
	jdsunset = jd_sun_alt(-(0.83+horiz),jdsunset,lat,longit);
	if(jdsunset > 0.) {
		oprntf("Sunset (%5.0f m horizon):  ",elev);
		print_time((jdsunset-zone(use_dst,stdz,jdsunset,*jdb,*jde)/24.),0);
		print_tz(jdsunset,use_dst,*jdb,*jde,zabr); 
	}
	else oprntf("Sunset not correctly computed; ");
	jdsunrise = jdmid + adj_time(rasun-hasunset-stmid)/24.;
	jdsunrise = jd_sun_alt(-(0.83+horiz),jdsunrise,lat,longit);
	if(jdsunrise > 0.) {
		oprntf("; Sunrise: ");
		print_time((jdsunrise-zone(use_dst,stdz,jdsunrise,*jdb,*jde)/24.),0);
		print_tz(jdsunrise,use_dst,*jdb,*jde,zabr);
	}
	if((jdsunrise > 0.) && (jdsunset > 0.)) {
		set_to_rise = (jdsunrise - jdsunset) * 24.;
		jdcent = (jdsunrise + jdsunset) / 2.;
	}
	else {
		oprntf(" Sunrise not correctly computed.");
		jdcent = -1.;
	}

	/* This block checks for and prints out 18-degree twilight
           as appropriate. */

	CHECK_TWI18: hatwilight = ha_alt(decsun,lat,-18.);
	if(hatwilight < -900.) {
		oprntf("\nFull darkness all day (sun below -18 deg).\n");
		twi_to_twi = 24.;
		goto DO_MOON_ETC;  /* certainly no 12-degree twilight */
	}
	if(hatwilight > 900.) {
		oprntf("\nSun higher than 18-degree twilight all night.\n");
		twi_to_twi = 0.;
		goto CHECK_TWI12;  /* but maybe 12-degree twilight occurs ...*/
	}               
	
	/* compute & print evening twilight and LST at eve. twilight */
	jdtwilight = jdmid + adj_time(rasun+hatwilight-stmid)/24.;  /* rough */
	jdtwilight = jd_sun_alt(-18.,jdtwilight,lat,longit);  /* accurate */
	jdetw = jdtwilight;
	if(jdtwilight > 0.) {
		oprntf("\nEvening twilight: ");
		print_time((jdtwilight-zone(use_dst,stdz,jdtwilight,*jdb,*jde)/24.),0);
		sid = lst(jdtwilight,longit);
		print_tz(jdtwilight,use_dst,*jdb,*jde,zabr);
		oprntf(";  LMST at evening twilight: ");
		put_coords(sid,0,0);
		oprntf("\n");
	}
	else oprntf("Evening twilight incorrectly computed.\n");                
	
        /* Now do morning twilight */
	jdtwilight = jdmid + adj_time(rasun-hatwilight-stmid)/24.;
	jdtwilight = jd_sun_alt(-18.,jdtwilight,lat,longit);
	jdmtw = jdtwilight;
	if(jdtwilight > 0.) {
		oprntf("Morning twilight: ");
		print_time((jdtwilight-zone(use_dst,stdz,jdtwilight,*jdb,*jde)/24.),0);         
		sid = lst(jdtwilight,longit);
		print_tz(jdtwilight,use_dst,*jdb,*jde,zabr);
		oprntf(";  LMST at morning twilight: ");
		put_coords(sid,0,0);
	}
	else oprntf("Morning twilight incorrectly computed.");
	if((jdetw > 0.) && (jdmtw > 0.)) twi_to_twi = 24. * (jdmtw - jdetw); 
	
        /* Now do same for 12-degree twilight */
	CHECK_TWI12: hatwilight = ha_alt(decsun,lat,-12.);
	if(hatwilight < -900.) {
		oprntf("\nSun always below 12-degree twilight...\n");
		goto DO_MOON_ETC;
	}
	if(hatwilight > 900.) {
		oprntf("\nSun always above 12-degree twilight...\n");
		goto DO_MOON_ETC;
	}               
	jdtwilight = jdmid + adj_time(rasun+hatwilight-stmid)/24.;
	jdtwilight = jd_sun_alt(-12.,jdtwilight,lat,longit);
	if(jdtwilight > 0.) {
		oprntf("\n12-degr twilight:");
		print_time((jdtwilight-zone(use_dst,stdz,jdtwilight,*jdb,*jde)/24.),0);
		print_tz(jdtwilight,use_dst,*jdb,*jde,zabr);
	}
	else oprntf("Evening 12-degree twilight incorrectly computed.\n");  

	jdtwilight = jdmid + adj_time(rasun-hatwilight-stmid)/24.;
	jdtwilight = jd_sun_alt(-12.,jdtwilight,lat,longit);
	if(jdtwilight > 0.) {
		oprntf(" -->");
		print_time((jdtwilight-zone(use_dst,stdz,jdtwilight,*jdb,*jde)/24.),0);         
		print_tz(jdtwilight,use_dst,*jdb,*jde,zabr);
		oprntf("; ");
	}
	else oprntf("Morning 12-degree twilight incorrectly computed.");

        DO_MOON_ETC: 
        if(jdcent > 0.) {
		oprntf("night center: ");
		print_time((jdcent-zone(use_dst,stdz,jdcent,*jdb,*jde)/24.),0);         
		print_tz(jdcent,use_dst,*jdb,*jde,zabr);
	}
	oprntf("\n\n");
	min_max_alt(lat,decmoon,&min_alt,&max_alt);  /* rough check -- occurs? */       
	if(max_alt < -(0.83+horiz)) {
		oprntf("Moon's midnight position does not rise.\n");
		jdmoonrise = -1.;
		goto MORE_MOON;
	}
	if(min_alt > -(0.83+horiz)) {
		oprntf("Moon's midnight position does not set.\n");
		jdmoonrise = 1.;
		goto MORE_MOON;
	}                       

	/* compute moonrise and set if they're likely to occur */

	hamoonset = ha_alt(decmoon,lat,-(0.83+horiz)); /* rough approx. */
	tmoonrise = adj_time(ramoon-hamoonset-stmid);
	tmoonset = adj_time(ramoon+hamoonset-stmid);
	jdmoonrise = jdmid + tmoonrise / 24.;
	jdmoonrise = jd_moon_alt(-(0.83+horiz),jdmoonrise,lat,longit,elevsea);
	jdmoonset = jdmid + tmoonset / 24.;
	jdmoonset = jd_moon_alt(-(0.83+horiz),jdmoonset,lat,longit,elevsea);
	if(fabs(set_to_rise) > 10.) moon_print = 0.65*set_to_rise;
		 else moon_print = 6.5;  

	/* more lousy flow of control, but it's nice to see the event which happens
	   first printed first .... */
	
	if(jdmoonset < jdmoonrise) goto PRINT_MOONSET;  /* OUCH!! */
	PRINT_MOONRISE: if((jdmoonrise > 0.) && (fabs(tmoonrise) < moon_print)) { 
	  /* print it if computed correctly and more-or-less at night */
		oprntf("Moonrise: ");
		print_time((jdmoonrise-zone(use_dst,stdz,jdmoonrise,*jdb,*jde)/24.),0);
		print_tz(jdmoonrise,use_dst,*jdb,*jde,zabr);
		oprntf("   ");
	}
	else if (jdmoonrise < 0.) oprntf("Moonrise incorrectly computed. ");
	if(jdmoonset < jdmoonrise) goto MORE_MOON;
	PRINT_MOONSET: if((jdmoonset > 0.) && (fabs(tmoonset) < moon_print)) {
		oprntf("Moonset : ");
		print_time((jdmoonset-zone(use_dst,stdz,jdmoonset,*jdb,*jde)/24.),0);
		print_tz(jdmoonset,use_dst,*jdb,*jde,zabr);
		oprntf("   ");
	}
	else if (jdmoonset < 0.) oprntf("Moonset incorrectly computed. ");
	if(jdmoonset < jdmoonrise) goto PRINT_MOONRISE; 

     MORE_MOON: ill_frac=0.5*(1.-cos(subtend(ramoon,decmoon,rasun,decsun)));
	oprntf("\nMoon at civil midnight: ");
	oprntf("illuminated fraction %5.3f\n",ill_frac);
	print_phase(jdmid);
	oprntf(", RA and dec: ");       
	put_coords(ramoon,2,0);
	oprntf(", ");
	put_coords(decmoon,1,1);
	oprntf("\n\n");

     /* print more information if desired */
     if(short_long == 2) {  /* wacky indenting here ... */
	oprntf("The sun is down for %4.1f hr; %4.1f hr from eve->morn 18 deg twilight.\n",
		set_to_rise,twi_to_twi);
	if((jdmoonrise > 100.) && (jdmoonset > 100.) &&
	   (twi_to_twi > 0.) && (twi_to_twi < 24.)) { 
	  /* that is, non-pathological */
		if((jdmoonrise > jdetw) && (jdmoonrise < jdmtw)) /* rises at night */
			oprntf("%4.1f dark hours after end of twilight and before moonrise.\n",
			  (24.*(jdmoonrise - jdetw)));
		if((jdmoonset > jdetw) && (jdmoonset < jdmtw)) /* sets at night */
			oprntf("%4.1f dark hours after moonset and before beginning of twilight.\n",
			  (24.*(jdmtw - jdmoonset)));
		if((jdmoonrise < jdetw) && (jdmoonset > jdmtw)) 
			oprntf("Bright all night (moon up from evening to morning twilight).\n");       
		if((jdmoonrise > jdmtw) && (jdmoonset < jdetw))
			oprntf("Dark all night (moon down from evening to morning twilight).\n");
	}
     }  /* closes the wacky indent. */
}

void print_circumstances(objra,objdec,objepoch,jd,curep,
	mura_arcs,mura_sec,mudec,sid,lat,elevsea,horiz)

double objra,objdec,objepoch,curep,mura_arcs,mura_sec,mudec,lat,horiz;
double jd,sid,elevsea;

/* Given object, site, and time information, prints the circumstances
   of an observation.  The heart of the "calculator" mode. */

{
	double objra_adj,objdec_adj,curra,curdec,ha,alt,az,secz,airm,par;
	char constelname[5];
	float ill_frac,tanz,howbad;
	double ramoon,decmoon,distmoon,ang_moon,
		georamoon, geodecmoon, geodistmoon,
		rasun,decsun,distsun,x,y,z,
		toporasun,topodecsun, 
		moon_alt,sun_alt,obj_moon,obj_lunlimb,sun_moon,Vmoon,
		obj_alt,eclong,eclat,tcor,vcor;
	short luneclcode = 0;

	oprntf("\n\nStd epoch--> RA:");
	put_coords(objra,3,0);
	oprntf(", dec: ");
	put_coords(objdec,2,1);
	oprntf(", ep %7.2f  ",objepoch);
        radec_to_constel(objra, objdec, objepoch, constelname);
	oprntf(" (in %s)\n",constelname);
	if((mura_sec != 0.) | (mura_arcs != 0.) |(mudec != 0.)) {
		objra_adj = objra + (curep-objepoch)* mura_sec/3600.;
		objdec_adj = objdec + (curep-objepoch)*mudec/3600.;
		oprntf("Adj for p.m: RA:");
		put_coords(objra_adj,3,0);
		oprntf(", dec: ");
		put_coords(objdec_adj,2,1);
		oprntf(", epoch %7.2f, equinox %7.2f\n", curep,objepoch);
		oprntf("(Annual proper motions:");
		mura_arcs = mura_sec * 15. *
		     cos(objdec / DEG_IN_RADIAN);
      		oprntf(" RA: %8.4f sec //%7.3f arcsec, ",mura_sec,mura_arcs);
		oprntf("dec: %7.3f)\n",mudec); 
		cooxform(objra_adj,objdec_adj,objepoch,curep,&curra,&curdec,
			XFORM_JUSTPRE,XFORM_FROMSTD);
	}
	else cooxform(objra,objdec,objepoch,curep,&curra,&curdec,
		XFORM_JUSTPRE, XFORM_FROMSTD);
	oprntf("Current  --> RA:");
	put_coords(curra,3,0);
	oprntf(", dec: ");
	put_coords(curdec,2,1);
	oprntf(", ep %7.2f\n",curep);
	ha = adj_time(sid - curra);
	oprntf("HA: ");
	put_coords(ha,2,1);
	alt=altit(curdec,ha,lat,&az,&par);
	obj_alt = alt;

	/* test size of sec z to avoid overflowing space provided */
	secz = secant_z(alt);
	if(secz >= 0.) {
          if(secz < 12.) 
 	      oprntf("; airmass = %8.3f",true_airmass(secz));
	   else if(secz <= 99.) 
	      oprntf("; sec.z = %8.3f",secz);
	   else oprntf(" Obj very near horizon.");
           if(secz > 3.) oprntf(" -- Large airmass!");
	}
	else oprntf(" -- BELOW HORIZON.");
	oprntf("\n");

	oprntf("altitude %6.2f, azimuth %6.2f, ",alt,az);
	oprntf("parallactic angle %4.1f",par);
	/* also give +- 180 ..... */
	if((par <= 180.) && (par > 0.)) oprntf("  [%4.1f]\n",par - 180.);
	 else oprntf("  [%4.1f]\n",par + 180.);
	/* And, also give a measure of how bad the refraction effects will be
	 * if you *don't* set to the parallactic angle ... */
	tanz = tan((90. - alt)/DEG_IN_RADIAN);
	howbad = fabs(tanz * sin(par/DEG_IN_RADIAN));
	if(alt > 3.) {
	  if(howbad < 100.) {
	     oprntf("Parallactic non-correction penalty: %5.2f\n\n",howbad);
	  }
	  else oprntf("Parallactic angle is very critical.\n");
	}

	accumoon(jd,lat,sid,elevsea,&georamoon,&geodecmoon,&geodistmoon,
			&ramoon,&decmoon,&distmoon);
	accusun(jd,sid,lat,&rasun,&decsun,&distsun,
		&toporasun,&topodecsun,&x,&y,&z);
	alt=altit(topodecsun,(sid-toporasun),lat,&az,&par);
	sun_moon = subtend(ramoon,decmoon,toporasun,topodecsun);
	ill_frac= 0.5*(1.-cos(sun_moon)); /* ever so slightly inaccurate ...
	   basis of ancient Greek limit on moon/sun distance ratio! */
	sun_moon = sun_moon * DEG_IN_RADIAN;
/*      printf("geocentric moon: ");
	put_coords(georamoon,4,0);
	printf("  ");
	put_coords(geodecmoon,3,1);
	printf(" %f \n",geodistmoon);  */
	sun_alt = alt;
	if(alt < -18.) oprntf("The sun is down; there is no twilight.\n");
	else {
		if (alt < -(0.83+horiz)) 
			oprntf("In twilight, sun alt %4.1f, az %5.1f ",alt,az);
		else oprntf("The sun is up, alt %4.1f, az %4.1f",alt,az);
		oprntf("; Sun at");
		put_coords(toporasun,3,0);
		oprntf(", ");
		put_coords(topodecsun,2,1);
		oprntf("\n");
		if(sun_moon < 1.5) solecl(sun_moon,distmoon,distsun);
			/* check for solar eclipse if it's close */
		if (alt < -(0.83+horiz))
		  oprntf("Clear zenith twilight (blue) approx %4.1f  mag over dark night sky.\n",
				ztwilight(alt));
	}       
	moon_alt=altit(decmoon,(sid-ramoon),lat,&az,&par);
	if(moon_alt > -2.) {
		oprntf("Moon:");
		put_coords(ramoon,2,0);
		oprntf(", ");
		put_coords(decmoon,1,1);
		oprntf(", alt %5.1f, az %5.1f;",moon_alt,az);
		oprntf("%6.3f illum.\n",ill_frac);
		print_phase(jd);
		obj_moon = DEG_IN_RADIAN * subtend(ramoon,decmoon,curra,curdec);
		if(fabs(obj_moon) > 10.) {
		  oprntf(".  Object is %5.1f degr. from moon.\n",obj_moon);
		}       
		else  {
			ang_moon = DEG_IN_RADIAN * asin(RMOON / (distmoon * EQUAT_RAD));
			if((obj_lunlimb = obj_moon - ang_moon) > 0.) 
			  oprntf(" ** NOTE ** Object %4.1f degr. from lunar limb!\n",obj_lunlimb);
			else oprntf(" ** NOTE ** You're looking AT the moon!\n");
		}
		if(sun_moon > 176.) 
			luneclcode = lunecl(georamoon,geodecmoon,geodistmoon,
			rasun,decsun,distsun);
		if((moon_alt > 0.) && (obj_alt > 0.5) && (sun_alt < -9.)) {
		  /*if it makes sense to estimate a lunar sky brightness */
		  Vmoon =                        
		     lunskybright(sun_moon,obj_moon,KZEN,moon_alt,
				obj_alt,distmoon); 
		  
		     oprntf("Lunar part of sky bright. = %5.1f V mag/sq.arcsec (estimated).\n",Vmoon);
		     if(luneclcode != 0) 
			oprntf(" NOT including effect of LUNAR ECLIPSE ...!\n"); 
		}
	}
	else {
		print_phase(jd);
		oprntf(".  The moon is down.\n");
	}
	eclipt(objra,objdec,objepoch,jd,&curep,&eclong,&eclat);
	if(fabs(eclat)<10.) {
		oprntf("Ecliptic latitude %4.1f; ",eclat);
		oprntf("watch for low-flying minor planets.\n");
	} 
	planet_alert(jd,curra,curdec,PLANET_TOL); 
	/* tolerance set to 3 degrees earlier. */
	helcor(jd,curra,curdec,ha,lat,elevsea,&tcor,&vcor);
	oprntf("Barycentric corrections: add %6.1f sec, %5.2f",tcor,vcor);
	oprntf(" km/sec to observed values.\n");
	oprntf("Barycentric Julian date = %14.6f\n",(jd+tcor/SEC_IN_DAY)); 
	printf("\nType command, 'f' for fast tour, '?' for a menu:");
}

void hourly_airmass(date,stdz,lat,longit,horiz,use_dst,objra,objdec,
  objepoch, mura_sec,mura_arcs,mudec)

/* Given a slew of information, prints a table of hourly airmass, etc.
   for use in scheduling observations.  Also prints sun and moon
   altitude when these are relevant.  Precesses coordinates as well. */

struct date_time date;
double stdz,lat,longit,horiz,objra,objdec,objepoch,mura_sec,mura_arcs,mudec;
short use_dst;

{
	
	double jd, jdb, jde, jdmid, curra, curdec, curep;
	double sid, ha, alt, az, secz, par, rasun, decsun, ramoon,
		decmoon,distmoon;       
 	double ill_frac,sepn;
        long int jdclong;
	int nch;
	short i, hr_span, dow;
	char obj_name[40], errprompt[80];
	char constelname[5];
        double hasset, jdsset, jdsrise, jdcent;
			
	if((date.y <= 1900) | (date.y >= 2100)) {
		printf("Date out of range - 1901 -> 2099\n");
		return;
	}
	printf("Name of object:");
	nch = get_line(obj_name);

	find_dst_bounds(date.y,stdz,use_dst,&jdb,&jde);
	date.h = 24; /* local midn */
	date.mn = 0;
	date.s = 0;
	jdmid = date_to_jd(date) + stdz/24.;
	/* first approx.-imperfect near time change */
	jdmid = date_to_jd(date) +
		zone(use_dst,stdz,jdmid,jdb,jde)/24.;
	curep = 2000.+(jdmid - J2000)/365.25;
	oprntf("\n\n*** Hourly airmass for %s ***\n\n",obj_name);
	if((use_dst != 0) && 
		((fabs(jdmid - jdb) < 0.5) || (fabs(jdmid - jde) < 0.5))) 
		oprntf("*** NOTE STD/DAYLIGHT TIME CHANGE TONIGHT ***\n");
	oprntf("Epoch %7.2f: RA ",objepoch);
	put_coords(objra,3,0);
	oprntf(", dec ");
	put_coords(objdec,2,1);
	radec_to_constel(objra, objdec, objepoch, constelname);
	oprntf("  (in %s)\n",constelname);
	oprntf("Epoch %7.2f: RA ",curep);
	cooxform(objra,objdec,objepoch,curep,&curra,&curdec,
			XFORM_JUSTPRE,XFORM_FROMSTD);
	put_coords(curra,3,0);
	oprntf(", dec ");
	put_coords(curdec,2,1);
	if((mura_sec != 0.) | (mura_arcs != 0.) | (mudec != 0)) 
	   oprntf("\n Caution .. proper motion ignored\n\n");
	else oprntf("\n\n");
	oprntf("At midnight: UT date ");
	print_calendar(jdmid,&dow);
	lpsun(jdmid,&rasun,&decsun);
	sid=lst(jdmid,longit);
	lpmoon(jdmid,lat,sid,&ramoon,&decmoon,&distmoon); /* close enuf */      
	ill_frac=0.5*(1.-cos(subtend(ramoon,decmoon,rasun,decsun)));
	sepn = DEG_IN_RADIAN * subtend(ramoon,decmoon,curra,curdec);
	oprntf(", Moon %4.2f illum, %3.0f degr from obj\n",ill_frac,sepn);
	planet_alert(jdmid,curra,curdec,PLANET_TOL);   /* better know about it .... */
	oprntf("\n  Local      UT      LMST");
	oprntf("      HA     secz   par.angl. SunAlt MoonAlt\n\n");
       
        /* figure out how much to print ... */

        hasset = ha_alt(decsun,lat,-0.83);
        jdsset = jd_sun_alt(-0.83,(jdmid - (12. - hasset) / 24.),lat,longit);
        jdsrise = jd_sun_alt(-0.83,(jdmid + (12. - hasset) / 24.),lat,longit);
        jdcent = (jdsset + jdsrise) / 2.;   /* center of night .. not local mid */
        hr_span = (short) (12. * (jdsrise - jdsset) + 0.5);
 /*       printf("jdsset %f jdsrise %f jdcent %f hr_span %d\n",
                jdsset,jdsrise,jdcent,hr_span);   ... diagnostic */
        jdclong = (long) ((24. * jdcent) + 0.5);  /* round to nearest hour */
        jdcent = jdclong / 24. + 0.00001;  /* add a hair to prevent "24 00"
              rounding ugliness in time table. */
 /*       printf("jdcent .... %f\n",jdcent);  ... diagnostic */
	for(i=(-1 * hr_span);i<=hr_span;i++) {
		jd = jdcent + i/24.;
		sid=lst(jd,longit);
		lpsun(jd,&rasun,&decsun);
		if(altit(decsun,(sid-rasun),lat,&az,&par) > 0.) goto SKIP;
		oprntf(" ");
		print_time((jd-zone(use_dst,stdz,jd,jdb,jde)/24.),0);
		oprntf("   ");
		print_time(jd,0);
		oprntf("   ");
	     /*	sid=lst(jd,longit); */
		put_coords(sid,0,0);
		oprntf("   ");
		ha = adj_time(sid - curra);
		put_coords(ha,0,0);
		oprntf("  ");
		alt=altit(curdec,ha,lat,&az,&par);
	        if(alt < -(horiz)) oprntf(" (down)");   
		else if(alt < 1.0) oprntf("(v.low)"); 
		else {
		     secz=secant_z(alt);
      		     oprntf(" %6.3f",secz);
		}
		oprntf("  %6.1f ",par);
		/* lpsun(jd,&rasun,&decsun); */
		alt=altit(decsun,(sid-rasun),lat,&az,&par);
		if(alt < -18.) oprntf("    ... ");
		   else oprntf("   %5.1f",alt);
		lpmoon(jd,lat,sid,&ramoon,&decmoon,&distmoon); /* close enuf */
		alt=altit(decmoon,(sid-ramoon),lat,&az,&par);
		if(alt < -2.)  oprntf("    ... \n");
		   else oprntf("   %5.1f\n",alt);
                SKIP: ;
	}
}

void print_params(date,enter_ut,night_date,stdz,lat,longit,site_name,
    elevsea,elev,use_dst,objra,objdec,objepoch,mura_sec,mura_arcs,
    mudec)

    struct date_time date;
    short enter_ut;
    short night_date;
    double stdz;
    double lat;
    double longit;
    char *site_name;
    double elevsea;
    double elev;
    short use_dst;
    double objra;
    double objdec;
    double objepoch;
    double mura_sec;
    double mura_arcs;
    double mudec;
		 
/* This simply prints a nicely formatted list of the *input* parameters
   without doing any computations.  Helpful for the perplexed user, and
   for checking things. */
{
	short dow;
	double jd;
	
/* better reflect the updating staus in 'l' ... */

#if SYS_CLOCK_OK == 1
	if(update_on == 1) {   
       		if(get_sys_date(&date,use_dst,enter_ut,
			night_date,stdz,update_delta) != 0) 
		       printf("Can't get system date! \n");
	}
#endif

	jd = date_to_jd(date);
	oprntf("\nCurrent INPUT parameter values:\n\n");
	oprntf("      DATE: ");
	print_calendar(jd,&dow); 
	oprntf("\n      TIME:");
	print_time(jd,3);
	if(night_date == 1) 
		oprntf("\n\nNIGHT_DATE:  ON    -- date applies all evening & next morning.\n");
	else oprntf("\n\nNIGHT_DATE:  OFF   -- date is taken literally.\n");
	if(enter_ut == 1) 
		oprntf("  UT_INPUT:  ON    -- input times taken to be UT.\n");
	else oprntf("  UT_INPUT:  OFF   -- input times taken to be local.\n");
	if(use_dst == 0) 
		oprntf("   USE_DST:   0    -- Standard time in use all year.\n");
	else oprntf("   USE_DST: %3d    -- Daylight savings used part of year.\n",      
			use_dst);
#if SYS_CLOCK_OK == 1
	if(update_on == 1) {
		oprntf("AUTO UPDATE: ON   -- output synched to system clock, offset %6.1f min.\n",update_delta);
	}
	else oprntf("AUTO UPDATE: OFF  -- system clock not automatically read on output.\n");
#endif

	if(night_date == 1) 
	oprntf("\n");
	oprntf("            RA: ");
	put_coords(objra,4,0);
	oprntf("\n           DEC: ");
	put_coords(objdec,3,1);
	oprntf("\n   INPUT EPOCH:   %8.2f\n",objepoch);
	if((mura_sec == 0.) && (mura_arcs == 0.) && (mudec == 0.))
		oprntf("PROPER MOTIONS:  OFF\n");
	else    {
          mura_arcs = mura_sec * 15. * cos(objdec / DEG_IN_RADIAN);
	  oprntf("PROPER MOTIONS: RA(sec)=%6.3f, RA(arc)=%6.2f, DEC=%6.2f\n",
		mura_sec,mura_arcs,mudec);
        }
	oprntf("\nSITE: %s\n",site_name);
	oprntf("      E.longit. = "); 
	put_coords(-15. * longit,1,1);
	oprntf(", latit. = ");
	put_coords(lat,1,1); 
	oprntf(" (degrees)\n");
	oprntf("      Standard zone = %3.0f hrs ",fabs(stdz));
	if(stdz >= 0.) oprntf(" West\n");
	  else oprntf(" East\n");
	oprntf("      Elevation above horizon = %4.0f m, True elevation = %4.0f m\n",
		elev,elevsea);
}

void print_menu() 
	{
	if(sclogfl != NULL) fprintf(sclogfl,"\n\n *** Menu Choices *** \n");
	oprntf("Circumstance calculator, type '=' for output.\n");
	oprntf("Commands are SINGLE (lower-case!) CHARACTERS as follows:\n");
	oprntf(" ? .. prints this menu; other information options are:\n");
	oprntf("i,f . 'i' prints brief Instructions and examples, 'f' fast tour\n");
	oprntf(" w .. prints info on internal Workings, accuracy & LEGALITIES\n");
	oprntf("TO SET PARAMETERS & OPTIONS, use these (note formats!):\n");
	oprntf(" r .. enter object RA, in hr min sec,  example: r 3 12 12.43 or 3:12:12.43\n");
	oprntf(" d .. enter object Dec in deg min sec, example: d -0 18 0 or -0:18:0\n");
	oprntf(" y .. enter date, starting with Year   example: y 1994 10 12\n");
#if SYS_CLOCK_OK == 1
     	oprntf("t,T: t = enter time, e.g.: t 22 18 02 [see 'g' and 'n']; T = right now+\n");
#else
     	oprntf(" t .. enter time, example: t 22 18 02 ; see 'g' and 'n' below\n");
#endif
     	oprntf(" n .. *toggles* whether date is used as 'evening' (default) or literal\n");
	oprntf(" g .. *toggles* whether time is used as Greenwich or local\n");         
	oprntf(" e .. enter Epoch used to interpret input coords (default = 2000)\n");
	oprntf(" p .. enter object Proper motions (complicated, follow prompts).\n");
	oprntf(" s .. change Site (again, follow prompts).\n");
	oprntf(" l .. Look at current parameter values (no computation).\n");
	oprntf("TO CALCULATE AND SEE RESULTS, use these commands: \n");
	oprntf(" = .. type out circumstances for specified instant of time, ra, dec, etc.\n");
	oprntf(" a .. type out night's Almanac for specified (evening) date\n");
	oprntf(" h .. type out Hourly airmass table for specified date, ra, dec\n");
        oprntf(" o .. tabulate Observability at 2-week intervals (at full&new moon)\n");
	oprntf(" m .. Major planets -- print 0.1 deg positions for specified instant\n");
        oprntf(" x .. eXtra goodies: log file, object files, other calculators; try x? \n");
	oprntf(" Q .. QUIT .. STOPS PROGRAM. ---> ");
}

void print_tutorial()
{
	if(sclogfl != NULL) fprintf(sclogfl,"\n  *** Fast Guided Tour listing ***\n\n");
	oprntf("FAST GUIDED TOUR: (type 'f' to see this again).\n\n");
	oprntf("To explore this program quickly, try the following (in order):\n\n");
	oprntf("--> Specify an evening date (e.g., 2002 April 4) and then display\n"); 
	oprntf("    the almanac for that night by typing \n");
	oprntf("        y 2002 4 4  a \n");
	oprntf("--> Specify an ra and dec and then make an hourly airmass table\n");
	oprntf("    by typing, for instance \n");
	oprntf("        r 15 38 29.2  d -0 01 02  h \n");
	oprntf("    (Notice how the date you specified earlier is still in effect.)\n");
	oprntf("--> Specify a time (e.g., 5 10 AM local time) and display the\n");
	oprntf("    circumstances at that instant by typing\n");
	oprntf("        t 5 10 00  =\n");
	oprntf("--> Type 'o' and follow the prompts to show observability\n");
	oprntf("    through a season (tabulated at new & full moon).\n");
	oprntf("--> Type 'l' to list the current input parameters and options.\n");
	oprntf("--> Type 'm' to list planetary positions.\n");
#if SYS_CLOCK_OK == 1
 	oprntf("--> Type 'T' to set time & date to right now (+ settable offset).\n");
        oprntf("--> Type 'x?' for eXtra goodies menu; log file, object lists etc..\n");
	oprntf("--> Read the help texts; '?' prints a short menu, 'i' gives\n");
#else
	oprntf("--> Read the help texts; '?' prints a short menu, 'i' gives\n");
#endif
	oprntf("    more complete information, 'w' tells you about inner workings.\n");
	oprntf("    An upper-case 'Q' exits the program.\n");
	oprntf(" Give a command: ");
}

void print_examples() {

	char cdum;

	if(sclogfl != NULL) fprintf(sclogfl,"\n\n *** On-Line Documentation listing ***");
	oprntf("\n\nMost parameters are entered as a single character followed\n");
	oprntf("by a value; you then type an EQUALS SIGN to compute circumstances\n");
	oprntf("for current site, time, and celestial position. FOR EXAMPLE:\n");
	oprntf("You observe an object at RA 19 02 23.3, dec -5 18 17\n");
	oprntf("epoch 1950, on 1993 June 18, at 21 32 18 local time, from\n");
	oprntf("the site you've selected.  To calculate circumstances, type\n\n");
	oprntf("r 19 02 23.3 d -5 18 17 e 1950 y 1993 6 18 t 21 32 18 =\n\n");
	oprntf("Note: coords, dates and times are entered as TRIPLETS OF NUMBERS\n");
	oprntf("separated by blanks (or colons), and date format is y m d.\n");
	oprntf("Commands are CASE SENSITIVE (R is not r)!  Carriage returns\n");
	oprntf("can go almost anywhere. Parameters remain in effect until you\n");
	oprntf("change them, so in this example typing\n");
	oprntf("t 23 32 18 =  \n");
	oprntf("would give circumstances for the same object two hours later.\n\n"); 

	oprntf("Typing 'a' prints an almanac for the currently specified date -- be\n");
	oprntf("careful you have the correct night if you've used 'n' or 'g'. (below).\n\n");

	oprntf("********** (Hit return to continue listing...) **********");
	scanf("%c",&cdum);
	scanf("%c",&cdum);  /* first one eats fossil input.. */

	oprntf("\n\nTwo commands invoke options for the interpretation of times and dates.\n");
        oprntf("[N. B. -- This issue is a little confusing, but important .... ]\n");
	oprntf("By default, the program interprets dates as 'evening', so late PM and\n");
	oprntf("early AM times refer to the same night.  Typing 'n' (for 'night date')\n");
	oprntf("TOGGLES bewteen this and literal dates.  Also by default, input\n");
	oprntf("times are in local (zone) time, but input times can be interpreted as UT;\n");
	oprntf("typing 'g' TOGGLES between these.  Invoking UT mode disables night dates.\n");
	oprntf("Note that 'g', 'n' (and the epoch option 'e') do not cause calculation, but\n");
	oprntf("only affect interpretation of the relevant parameters at compute time.\n\n");

	oprntf("Typing 'l' (look) causes all the input parameters to be typed out without\n");
        oprntf("any computation, which is especially helpful if you're lost by 'g' or 'n'.\n\n");

	oprntf("Typing 'h' creates a table of airmass, etc. at hourly intervals through\n");
	oprntf("the night.  You may wish to redirect output (using system) to print.\n\n");

	oprntf("Typing 'o' computes observability of an object for each dark\n");
	oprntf("and bright run through a range of dates.\n\n");

	oprntf("********** (Hit return to continue listing...) **********");
	scanf("%c",&cdum);  

	oprntf("\nTyping 'm' displays rough (see 'w' option) positions of major planets.\n\n");

	oprntf("To change observing sites, type 's' and answer the prompts.\n\n");

#if SYS_CLOCK_OK == 1
	oprntf("Typing 'T' sets BOTH the TIME and DATE to 'now' using the system clock.\n");
	oprntf("You can also specify a number of minutes into the future to set.  The status\n");
	oprntf("of the 'n' and 'g' options are properly taken into account.\n\n");
	oprntf("Typing 'xU' (one of many 'x' subcommands) causes system clock to be\n");
	oprntf("read every time output is computed, auto-updating.\n\n");
#endif
	oprntf("To enable an approximate proper motion correction, type 'p' and\n");
	oprntf("answer the prompts *carefully*.  Re-enter as zero to turn off.\n\n");
	oprntf("Typing 'x' invokes many useful commands -- log file, object files, etc!\n\n");
	oprntf("Type command, Q (upper case) to stop, or ? for a menu:");
}

void print_accuracy() 

{

	char cdum;
	
        if(sclogfl != NULL) fprintf(sclogfl,"\n\n");
	oprntf("ACCURACY INFO:\n\n");
	oprntf("The distinctions between UTC, UT1, TDT (etc.) are ignored\n");
	oprntf("except that a rough correction to TDT is used for the moon.\n");
	oprntf("The solar ephemeris used is good to a few arcsec.  Moon positions\n");
	oprntf("positions are topocentric and +- about 30 arcsec, hence solar\n"); 
	oprntf("eclipse paths are +- 50 km and +- 1 min. All rise/set times are\n");
	oprntf("computed to about +-1 min; non-level horizon, site elevation,\n");
	oprntf("and refraction uncertainties are often larger than this.\n\n");
	oprntf("The lunar sky brightness model assumes ideal atmospheric conditions;\n");
	oprntf("true lunar contributions to sky will vary widely.  To compare\n");
	oprntf("a dark site has about V=21.5 mag/sq.arcsec (variable)! Twilight\n");
	oprntf("brightness prediction is for blue, and only very approximate.\n\n");
	oprntf("The planetary calculations are truncated, but the error should\n");
	oprntf("seldom exceed 0.1 degree; MV are best(1'), MJSU ok, Pluto worst.\n");
	oprntf("\nNote that the local sidereal time given is Mean, not true,\n");
	oprntf("and that it assumes the input is true UT, not UTC (< 1 second)\n");
        oprntf("********** (Hit return to continue listing...) **********");
	scanf("%c",&cdum);
	scanf("%c",&cdum);  /* first one usually eats fossil input.. */
	oprntf("\nDaylight savings time, if selected, is established using a\n");
	oprntf("site-specific convention (e.g., USA).  Beware of ambiguities\n");
	oprntf("and nonexistent times when the clock is reset.  If necessary,\n");
	oprntf("use the 'g' option and enter times and dates as Greenwich (UT),\n");
	oprntf("or disable DST in site params (see discussion under 'i').\n\n");
	oprntf("The *precession* routine used is a 'rigorous'");
	oprntf(" rotation matrix routine,\ntaken from L. Taff's Computational");
	oprntf(" Spherical Astronomy.\nIt uses IAU1976 constants, is good to < 1 arcsec in 50 years,\n");
	oprntf("and has no troubles near the pole.  Proper motion corrections\n");
	oprntf("are done crudely as x = x0 + mu * dt; this is inaccurate near\n");
	oprntf("the poles.  Use another routine if sub-arcsec accuracy is critical.\n");
	oprntf("'xa' --> apparent place(nutation, aberration, no parallax)+-0.1 arcsc\n");
        oprntf("********** (Hit return to continue listing...) **********");
	scanf("%c",&cdum); 
	oprntf("\n\nThe parallactic angle is described by Filippenko (1982, PASP 94, 715),\n");
	oprntf("This implementation avoids a root-choice problem in his formulae.\n");

	oprntf("\nThe barycentric ('heliocentric') corrections are computed using\n");
	oprntf("an elliptical earth orbit with a few periodic perturbations\n");
	oprntf("including lunar recoil.  The helio-to-barycentric transformation\n");
	oprntf("uses the same algorithms as the planetary postions.  Overall max error:\n");
	oprntf("< 0.2 sec and < 0.005 km/s.  Velocity corrn. includes earth rotation.\n");
	oprntf("\nThe galactic coordinate routine is rigorously accurate, and\n");
	oprntf("precesses to 1950 before transforming.  The ecliptic coord.\n");
	oprntf("routine is for coordinates of date and is good to < 1 arcsec.\n\n");
	oprntf("These routines are **not necessarily correct** at times very far from\n");
	oprntf("the present (2000).  The program rejects input outside 1900-2100.\n");
	oprntf("\nWhen porting to a new machine, run tests to ensure\n");
	oprntf("correctness and accuracy.  Experience shows that compiler\npeculiarities arise ");
	oprntf("surprisingly often.\n\n");
        oprntf("********** (Hit return to continue listing...) **********");
	scanf("%c",&cdum); 
	oprntf("CAUTIONS, LEGALITIES:\n");
	oprntf("Many routines take a time argument which is a double-precision\n");
	oprntf("floating-point julian date; on most workstations, this gives\n");
	oprntf("time resolution of < 0.1 second.  When porting\n");
	oprntf("to another machine or compiler, test that the accuracy is\n");
	oprntf("sufficient.  \n\n");
	oprntf("I cannot guarantee that this program is bug-free, and caution\n");
   /* although at this time I really don't know of any bugs ... ! */
	oprntf("that not all routines are thoroughly precise and rigorous.\n");
	oprntf("The user of this program is responsible for interpreting \n");
	oprntf("results correctly.  I disavow any legal liability for damages\n");
	oprntf("caused by use of this program.\n\n");
	oprntf("Program copyright 2000, John Thorstensen, Dartmouth College\n");
	oprntf("Permission hereby granted for scientific/educational use.\n");
	oprntf("For-profit users must license.  Please communicate problems or\n");
	oprntf("suggestions to the author, John Thorstensen, Dept. Physics\n");
	oprntf("and Astronomy, Dartmouth College, Hanover NH 03755\n");
	oprntf("John.Thorstensen@Dartmouth.edu\n");
	oprntf("\nType command, or ? to see menu again.:");
}

void print_legalities() 

{	
	printf("\nThis has been superseded -- this shouldn't print out.\n");
}


void ephemgen(ra, dec, ep, lat, longit)
	double ra, dec, ep, lat, longit;

{

/* Prompts for elements of an ephemeris, then 
   generates a printout of the *geocentric*
   ephemeris of an object.  Uses current values of ra and dec; 
   ignores observatory parameters, which are much more important
   for velocity than for time.  Not *strictly* correct in that
   earth's position in helio-to-geo conversion is computed (strictly
   speaking) for heliocentric jd, while ephemeris is for (presumably)
   geocentric jd.  This makes a difference of at most 30 km/s  over c,
   or about 1 in 10**4, for the heliocentric correction.  Pretty trivial.

   Observatory parameters are used in determining observability of 
   phenomenon.
*/

  double per; /* ephemeris -- t0 + per * cycle count */
  int ibeg, iend, cycles;
  double curra, curdec, curep;
  double max_secz, min_alt, max_alt, min_ok_alt, sun_alt, 
		max_ok_sun, secz, par;
  double sid, ha, altitude, az, rasun, decsun, hasun;
  double jd, jdstart, jdend, tcor, vcor, errbar;
  struct date_time tempdate;
  int valid_date = 0;
  short dow;

  printf("Prints geocentric times of repeating phenom, e.g. stellar eclipses.\n");
  printf("Not suitable for very short periods (seconds) -- too many events.\n");
  printf("You may optionally specify that only events observable from your site\n");
  printf("be printed.\n\n");  
  printf("Current period %f +- %f, epoch %f +- %f\n",
		star_period, star_perr, star_tzero, star_terr);
  printf("Give heliocentric period (days); negative keeps all the old values:");
  scanf("%lf",&per);
  if(per > 0.) { 
	star_period = per;
	printf("Give error bar for period:");
	scanf("%lf",&star_perr);
  	printf("Give heliocentric jd of time zero, *all* the digits!:");
  	scanf("%lf",&star_tzero);
  	if(star_tzero < 2000000.) {
		printf("Exiting -- JD has to be greater than 2 million!\n");
		return;
	}
	printf("Give error bar for time zero:");
	scanf("%lf",&star_terr);
  }
  printf("Listing will start and end at 0h UT ... \n");
  printf("Give starting UT date for listing, yyyy mm dd, neg year exits:");
  if((valid_date = get_date(&tempdate)) < 0) return;
/*  printf("Give starting ut, h m s:");
  get_time(&tempdate);   Too much of pain to enter .. force to ut = 0 */
  tempdate.h = 0.;
  tempdate.mn = 0.;
  tempdate.s = 0.;
  jdstart = date_to_jd(tempdate);
  printf("Give ending date for listing:");
  get_date(&tempdate);
/*  printf("Give ending time:");
  get_time(&tempdate);  */
  jdend = date_to_jd(tempdate);
  printf("Give maximum airmass for listing, negative to print all:");
  scanf("%lf",&max_secz);
  if(max_secz > 0.) min_ok_alt = 90. - DEG_IN_RADIAN * acos(1/max_secz);
  else min_ok_alt = -100.;
  curep = 2000. + (jdstart - J2000)/365.25;    /* precess to check observ. */
  cooxform(ra,dec,ep,curep,&curra,&curdec,
		XFORM_JUSTPRE,XFORM_FROMSTD);
  min_max_alt(lat,curdec,&min_alt,&max_alt);
  if(max_alt < min_ok_alt) {
  oprntf("Object never rises this high -- max alt %5.1f deg -- returning!\n",
           max_alt);
	return;
  }

  printf("Give maximum allowable sun altitude to print, 90 to print all:");
  scanf("%lf",&max_ok_sun);

  oprntf("\nObject RA = ");
  put_coords(ra,3,0);
  oprntf(", dec = ");
  put_coords(dec,2,1);
  oprntf(" epoch %6.1f\n",ep);
  oprntf("Site long = ");
  put_coords(longit,3,1);
  oprntf("W (hr.mn.sec), lat = ");
  put_coords(lat,2,1);
  oprntf("N \n");
  oprntf("\nEphemeris: %15.6f + %14.8f E (Heliocentric)\n\n",
		star_tzero, star_period);
  oprntf("Listing of events in interval from ");
  print_calendar(jdstart,&dow);
  oprntf(" UT to ");
  print_calendar(jdend,&dow);
  oprntf(" UT.\n");
  if(min_ok_alt > 0.) oprntf("Only events at sec(z) < %5.2f are shown.\n",
        max_secz);
  if(max_ok_sun < 89.) 
	oprntf("Only events with sun altit. < %6.1f degrees are shown.\n\n",
	max_ok_sun);
  oprntf("First columns show cycle number, GEOCENTRIC jd, date and time (UT)\n");
  oprntf("(including accumulated error in minutes)\n");
  oprntf("Then: HA, sec(z), 'ngt' if night, sun alt if twilit, 'day' if day.\n");
  ibeg = (jdstart - star_tzero) / star_period;  /* truncate */
  iend = (jdend - star_tzero) / star_period;    /* truncate again */

  for(cycles = ibeg; cycles <= iend; cycles++) {
 	jd = star_tzero + (double) cycles * star_period;
	errbar = (star_terr + (double) cycles * star_perr) * 1440.;
	curep = 2000. + (jd - J2000)/365.25;    /* precess to current */
        cooxform(ra,dec,ep,curep,&curra,&curdec,XFORM_JUSTPRE,XFORM_FROMSTD);
	helcor(jd, curra, curdec, 0., 0., 0., &tcor, &vcor);
/*      oprntf("Raw = %12.4f, Correction: %f\n",jd,tcor);   diagnostic */
        jd = jd - tcor / SEC_IN_DAY;  /* subtract to convert to geo .*/
        sid = lst(jd,longit);
        ha = adj_time(sid - curra);
        altitude = altit(curdec,ha,lat,&az,&par);
	secz = secant_z(altitude);
        if(altitude >= min_ok_alt) {
           lpsun(jd,&rasun,&decsun);  /* lpsun plenty good enough */
	   hasun = lst(jd,longit) - rasun;
	   sun_alt = altit(decsun,hasun,lat,&az,&par);
           if(sun_alt < max_ok_sun) {
              oprntf(" %d  %12.4f = ",cycles,jd);
	      print_calendar(jd,&dow);
	      oprntf(",");
              print_time(jd,3);
	      oprntf(" (+- %5.2f m) UT",errbar);
	      put_coords(ha,0,0);
	      if((secz < 10.) && (secz > 0.)) oprntf(" %5.2f",secz);
              else if (secz > 0.) oprntf(" v.low");
	      else oprntf(" down!");
	      if(sun_alt <= -18.) oprntf(" ngt");
	      else if(sun_alt <= -0.53) oprntf(" %3.0f",sun_alt);
	      else oprntf(" day ");
	      oprntf("\n");
	   }
       }
  }
  oprntf("Listing done.  Type ? if you need a menu.\n");
}

#define ALT_3  19.47  /* 19.47 degrees altitude => sec z = 3 */
#define ALT_2  30.
#define ALT_15 41.81
#define SID_RATE 1.0027379093  /* sidereal / solar rate */

double hrs_up(jdup, jddown, jdeve, jdmorn) 

	double jdup, jddown, jdeve, jdmorn;
{

   /* If an object comes up past a given point at jdup,
      and goes down at jddown, and evening and morning
      are at jdeve and jdmorn, computes how long
      object is up *and* it's dark.  ... Written as
      function 'cause it's done three times below. */

   double jdup2, jddown0;  /* for the next ... or previous ...
                            time around */

   if(jdup < jdeve) {
     if(jddown >= jdmorn)   /* up all night */
	return ((jdmorn - jdeve) * 24.);
     else if(jddown >= jdeve)  {
	/* careful here ... circumpolar objects can come back *up*
	   a second time before morning.  jdup and jddown are
	   the ones immediately preceding and following the upper
	   culmination nearest the center of the night, so "jdup"
	   can be on the previous night rather than the one we want. */
       jdup2 = jdup + 1.0/SID_RATE;
       if(jdup2 > jdmorn)  /* the usual case ... doesn't rise again */
	   return ((jddown - jdeve) * 24.);
       else return(((jddown - jdeve) + (jdmorn - jdup2)) * 24.);
     }
     else return(0.);   
   }
   else if(jddown > jdmorn) {
     if(jdup >= jdmorn) return(0.);
     else {
	/* again, a circumpolar object can be up at evening twilight
	    and come 'round again in the morning ... */
        jddown0 = jddown - 1.0/SID_RATE;
	if(jddown0 < jdeve) return((jdmorn - jdup) * 24.);
        else return(((jddown0 - jdeve) + (jdmorn - jdup)) * 24.);
     }
   }
   else return((jddown - jdup)*24.);  /* up & down the same night ...
      might happen a second time in pathological cases, but this will
      be extremely rare except at very high latitudes.  */
}

void print_air(secz,prec) 
    	double secz;
	short prec;
{
   if((secz > 0.) && (secz < 100.)) {
	if(prec == 0) oprntf("%5.1f ",secz);
	         else oprntf("%5.2f ",secz);
   }
   else if(secz > 0.) oprntf("v.low ");
   else if(secz < 0.) oprntf(" down ");
}

void print_ha_air(ha,secz,prec1,prec2) 
	double ha, secz;
	short prec1, prec2;
{
	oprntf("  ");
 	put_coords(ha,prec1,1);
	oprntf(" ");
	print_air(secz,prec2);
}

void obs_season(ra, dec, epoch, lat, longit)

   double ra, dec, epoch, lat, longit;	

/* prints a table of observability through an observing
   season.  The idea is to help the observer come up
   with an accurately computed "range of acceptable
   dates", to quote NOAO proposal forms ... */

{
   int valid_date, nlun, nph, nch;
   char obj_name[40];
   short dow;
   double sun_twi;
   double jdstart, jdend, jd, jdtrunc, jdevedate;
   double jdmid, jdcent, jdeve, jdmorn, midnfrac;
   double steve,stmid,stcent,stmorn;
   double hatwi, hasun, hacent, haeve, hamorn;
   double seczcent, seczeve, seczmorn, par;
   double dt15,dt2,dt3;
   double jdtrans, jd3_1, jd3_2, jd2_1, jd2_2,
	  jd15_1, jd15_2, rasun, decsun;
   double hrs_3, hrs_2, hrs_15;
   double min_alt, max_alt, altitude, az;
   double curep, curra, curdec;
   struct date_time tempdate;

   printf("\nThis types out a summary of the observability of your object\n");
   printf("through the observing season.  Observability is summarized\n");
   printf("at new & full moon through a range of dates you specify.\n");
   printf("(a 24-line display holds about a 6-month range.)\n\n");
   if(fabs(lat) > 70.) {

	oprntf("NOTE: The site's geographical latitude is so high that there's\n");
	oprntf("a possibility that 'darkness' will last all day (depending\n");
        oprntf("on how you define the end of twilight below.)  In this case\n");
	oprntf("the HA and sec z at 'twilight' will be for +-12 hours from\n");
	oprntf("sun's lower culmination.\n\n");
   }

   printf("Give approx. starting date, yyyy mm dd, neg year exits:");

  if((valid_date = get_date(&tempdate)) < 0) return;
  tempdate.h =  0.;
  tempdate.mn = 0.;
  tempdate.s = 0.;
  jdstart = date_to_jd(tempdate);
  printf("Give ending date for listing:");
  if((valid_date = get_date(&tempdate)) < 0) return;
  jdend = date_to_jd(tempdate);
  printf("You must now specify the altitude of the sun which defines twilight.\n");
  printf("Type -18. , or give a non-standard choice:");
  scanf("%lf",&sun_twi);

  printf("Name of object:");     /* for labeling redirected output */
  nch = get_line(obj_name);

  jd = jdstart - lun_age(jdstart,&nlun);
  nph = 0;  /* jd is adjusted to last previous new moon */

  jdcent = (jd + jdend) / 2.;  /* temporary use of var */
  curep = 2000. + (jdcent - J2000)/365.25;    /* precess to check observ. */
  cooxform(ra,dec,epoch,curep,&curra,&curdec,XFORM_JUSTPRE,XFORM_FROMSTD);
  min_max_alt(lat,curdec,&min_alt,&max_alt);

  oprntf("\n          *** Seasonal Observability of %s ***\n",obj_name);
  oprntf("\n     RA & dec: ");
  put_coords(ra,3,0);
  oprntf(", ");
  put_coords(dec,2,1);
  oprntf(", epoch %6.1f\n", epoch);
  /*  printf("\nObject RA&dec: ");
  put_coords(curra,3,0);
  oprntf(", ");
  put_coords(curdec,2,0);
  oprntf(", epoch %6.1f\n", curep);  */

  oprntf("Site long&lat: ");
  put_coords(longit,3,1);
  oprntf(" (h.m.s) West, ");
  put_coords(lat,2,1);
  oprntf(" North.\n\n");


  if(max_alt < ALT_15)
    oprntf("NOTE -- This object is always at sec.z > %5.2f\n",
	1./cos((90. - max_alt) / DEG_IN_RADIAN));
  if(min_alt > ALT_3)
    oprntf("NOTE -- This object is always at sec.z < %5.2f\n",
	1./cos((90. - min_alt) / DEG_IN_RADIAN));

  oprntf("Shown: local eve. date, moon phase, hr ang and sec.z at (1) eve. twilight,\n");
  oprntf("(2) natural center of night, and (3) morning twilight; then comes number of\n");
  oprntf("nighttime hours during which object is at sec.z less than 3, 2, and 1.5.\n");
  oprntf("Night (and twilight) is defined by sun altitude < %4.1f degrees.\n\n",sun_twi);

  oprntf(" Date (eve) moon      eve            cent           morn");
  oprntf("     night hrs@sec.z:\n");
  oprntf("                    HA  sec.z      HA  sec.z      HA  sec.z");
  oprntf("    <3   <2   <1.5\n");


  while(jd <= jdend) {
       if(nph == 0) nph = 2;
       else if(nph == 2) {
	  nlun++;
	  nph = 0;
       }
       flmoon(nlun,nph,&jd);

       /* Take care to compute for the date nearest full or new ...
	  people may use this as a lunar calendar, which it sort of
	  isn't. */

       midnfrac = 0.5 + longit/24.;  /* rough fractional part of jd
		       for local midnight ... */
       jdtrunc = (double)((long) jd);
       jdmid = jdtrunc + midnfrac;
       if((jd - jdmid) > 0.5) jdmid = jdmid + 1.;
       else if((jd - jdmid) < -0.5) jdmid = jdmid - 1.; 

       lpsun(jdmid,&rasun,&decsun);
       hasun = adj_time(lst(jdmid,longit) - rasun); /* at rough midn. */
       if(hasun > 0.) jdcent = jdmid + (12. - hasun) / 24.;
       else jdcent = jdmid - (hasun + 12.) / 24.;
	  /* jdcent is very close to sun's lower culmination
	      -- natural center of night */
   
       hatwi = ha_alt(decsun,lat,sun_twi);
       if(hatwi > 100.) {
	      jdeve = -1.;    /* signal -- no twilight */
              jdmorn = -1.;
       }

       else {
	  if(hatwi < -100.) {  /* always night */
	        jdeve = jdcent - 0.5; /* sensible, anyway  */
	        jdmorn = jdcent + 0.5;
          }
          else {
		jdmorn = jd_sun_alt(sun_twi,(jdcent +(12.-hatwi)/24.),
	           lat,longit);
          	jdeve = jd_sun_alt(sun_twi,(jdcent - (12.-hatwi)/24.),
	           lat,longit);      
          }
          haeve = adj_time(lst(jdeve,longit) - curra);
          altitude = altit(curdec,haeve,lat,&az,&par);
	  seczeve = secant_z(altitude);

          hamorn = adj_time(lst(jdmorn,longit) - curra);
          altitude = altit(curdec,hamorn,lat,&az,&par);
	  seczmorn = secant_z(altitude);
       }

       hacent = adj_time(lst(jdcent,longit) - curra);
       altitude = altit(curdec,hacent,lat,&az,&par);
       seczcent = secant_z(altitude);

       jdtrans = jdcent - hacent / (SID_RATE * 24.);
	  /* this will be the transit nearest midnight */
       if((min_alt < ALT_3) && (max_alt > ALT_3)) {
	 /* if it makes sense to compute the times when
	    this object passes three airmasses ... */
	  dt3 = ha_alt(curdec,lat,ALT_3) / (SID_RATE * 24.);
	  jd3_1 = jdtrans - dt3;	     
	  jd3_2 = jdtrans + dt3;
       }
       else {
	  jd3_1 = 0.;  jd3_2 = 0.;
       }
        if((min_alt < ALT_2) && (max_alt > ALT_2)) {
	 /* if it makes sense to compute the times when
	    this object passes two airmasses ... */
	  dt2 = ha_alt(curdec,lat,ALT_2) / (SID_RATE * 24.);
	  jd2_1 = jdtrans - dt2;	     
	  jd2_2 = jdtrans + dt2;
       }
       else {
	  jd2_1 = 0.;  jd2_2 = 0.;
       }
       if((min_alt < ALT_15) && (max_alt > ALT_15)) {
	 /* if it makes sense to compute the times when
	    this object passes 1.5 airmasses ... */
	  dt15 = ha_alt(curdec,lat,ALT_15) / (SID_RATE * 24.);
	  jd15_1 = jdtrans - dt15;	     
	  jd15_2 = jdtrans + dt15;
       }
       else {
	  jd15_1 = 0.;  jd15_2 = 0.;
       }

   /* Now based on times of twilight and times at which object
      passes the airmass points, figure out how long the object
      is up at night ... */

      if(jdeve > 0.) {  /* if twilight occurs ... */
         if(jd2_1 != 0.)
	   hrs_2 = hrs_up(jd2_1, jd2_2,jdeve,jdmorn);
         else if(min_alt > ALT_2)
	   hrs_2 = 24. * (jdmorn - jdeve);  /* always up ... */
         else hrs_2 = 0.;                   /* never up ... */

         if(jd3_1 != 0.)
	   hrs_3 = hrs_up(jd3_1, jd3_2,jdeve,jdmorn);
         else if(min_alt > ALT_3)
	   hrs_3 = 24. * (jdmorn - jdeve);
         else hrs_3 = 0.;

         if(jd15_1 != 0.)
           hrs_15 = hrs_up(jd15_1, jd15_2,jdeve,jdmorn);
         else if(min_alt > ALT_15)
	   hrs_15 = 24. * (jdmorn - jdeve);
         else hrs_15 = 0.;
      }
      else if(jdeve > -500.) {  /* twilight all night */
	  hrs_2 = 0.;
	  hrs_3 = 0.;
	  hrs_15 = 0.;
      }


       /* now print the table itself.  Take pains to
	  get the correct *evening* date .... */

       jdevedate = jdcent - longit/24. - 0.25;
       planet_alert(jdcent,curra,curdec,PLANET_TOL);
         /* better know about it ... */

       print_calendar(jdevedate,&dow);

       /* space table correctly -- what a pain ! */
       caldat(jdevedate,&tempdate,&dow);
       if(tempdate.d < 10) oprntf(" ");

       if(nph == 0) oprntf("   N");
       else oprntf("   F");
       if(jdeve > 0.) print_ha_air(haeve,seczeve,0,0);
       else oprntf(" twi.all.nght! ");
       print_ha_air(hacent,seczcent,0,0);       
       if(jdmorn > 0.) print_ha_air(hamorn,seczmorn,0,0);
       else oprntf(" twi.all.nght! ");

       oprntf(" %4.1f  %4.1f  %4.1f \n",hrs_3,hrs_2,hrs_15);

   }
   printf("Listing done.  'f' gives tutorial, '?' prints a menu.\n");
}

#if SYS_CLOCK_OK == 1
#include <time.h>
#endif

#if SYS_CLOCK_OK == 1

int get_sys_date(date, use_dst, enter_ut, night_date, stdz, toffset) 

	struct date_time *date;
	short use_dst, enter_ut, night_date;
        double stdz, toffset;
{
	/* Reads the system clock; loads up the date structure
           to conform to the prevailing conventions for the interpretation
           of times.  Optionally adds "toffset" minutes to the system
           clock, as in x minutes in the future. */

	time_t t, *tp;
	struct tm *stm;
        double jd, jdb, jde;
        short dow;

	tp = &t;  /* have to initialize pointer variable for it to
                     serve as an argument. */

	t = time(tp);
	if(t == -1) {
		oprntf("SYSTEM TIME UNAVAILABLE!  Date & time left alone.\n");
		return(-1);
        }
	stm = localtime(&t);
	date->y = (short) (stm->tm_year + 1900);
	date->mo = (short) (stm->tm_mon + 1);
	date->d = (short) (stm->tm_mday);
	date->h = (short) (stm->tm_hour);
	date->mn = (short) (stm->tm_min);
	date->s = (float) (stm->tm_sec);

	if(toffset != 0.) {
	   jd = date_to_jd(*date);
	   jd = jd + toffset / 1440.;
	   caldat(jd,date,&dow);
        }

        if(enter_ut == 1)  { /* adjust if needed */
           find_dst_bounds(date->y,stdz,use_dst,&jdb,&jde);
	   jd = date_to_jd(*date);
	   jd = jd + zone(use_dst,stdz,jd,jdb,jde)/24.;
           caldat(jd,date,&dow);
        }
	oprntf("\nDATE AND TIME SET AUTOMATICALLY USING SYSTEM CLOCK ... to\n");
	print_all(date_to_jd(*date));
	if(enter_ut == 0) oprntf(" local time\n");
        else oprntf(" Universal time\n"); 
      	if((night_date == 1) && (date->h < 12)) {
		date->d = date->d - 1;
		oprntf("(night-date is on, so internal value of date set to last evening)\n"); 		   
	}
#if SYS_CLOCK_OK == 1
	if(update_on == 1) 
	   oprntf("NOTE: time, date are automatically updated, offset %4.0f minutes.\n",update_delta);
	else oprntf("NOTE: time, date NOT auto-updated -- 'xU' toggles this, 'T' reads clock.\n");
#else
	oprntf("(NOTE: time does NOT automatically update; 'T' reads clock again.)\n");   
#endif
        return(0); /* success */
}

#endif     /* previous routine depends on sys clock */

/* Sorting package ...
  adapted from J. W. J. Williams, "Collected 
  Algorithms from Com. Assoc. Computing Machinery", 
  Algorithm # 232.  Heapsort. 

 Heaps are all inverted (using greater-than rather than
 less-than tests), since in last step the algorithm 
 strips off the first element of the array and copies it
 to the end ... this is then the greatest element, so sort
 comes out in correct order. 

 This algorithm coded from public-domain sources so as to
 avoid copyright restrictions from Numerical Recipes in 
 code originally distributed with skycalc.  It's a lot
 longer than then Numerical Recipes routine, and probably
 slower, but (a) it's way fast enough for our purposes,
 (b) it does work, and (c) it's completely legal.

*/

void inheap(float *a, int *ind, int *n) {

/* assumes that first n elements of a are already formed
   into a heap such that a[ind[1]] .... a[ind[n]] are a heap.
   Assumes that next element (which has ind[n+1] = n+1) is
   not yet in heap.  Rearranges ind array to place next 
   element on heap. */

	int i, j;
	
	*n = *n + 1;
	i = *n;

	SCAN: if(i > 1) {   
		j = i / 2;
		if(a[*n] > a[ind[j]]) {  
			ind[i] = ind[j];
			i = j;
			goto SCAN;
		}
	}
	ind[i] = *n;
}

void setheap(float *a, int *ind, int n) {
   
/* arranges a[1] through a[n] to form a heap. */

	int j, i;

	j = 1; 
	while (j < n) {
		inheap(a,ind,&j);
	}
}

void swapheap(float *a, int *ind, int n, int inputind, int *outind) {

/* for elements of array a, in which a[ind[1]] through a[ind[n]] already
  form a heap, returns the index of the largest (or smallest)
   element in *outind, and then adds a[inputind] to the heap */

	double temp, temp1;
	int i, j, intemp, intemp1;
	
	if(a[inputind] >= a[ind[1]]) *outind = inputind;
	
	else {
		i = 1;  
		ind[n+1] = inputind;
		*outind = ind[1];
		SCANSW: j = i + i; 
		if(j <= n) {
			temp = a[ind[j]];
			intemp = ind[j];
			temp1 = a[ind[j + 1]];
			intemp1 = ind[j+1];
			if(temp1 > temp) {
				temp = temp1;
				intemp = intemp1;
				j++;
			}
			if(temp > a[inputind]) {
				ind[i] = intemp;
				i = j;
				goto SCANSW;
			}
		}
		ind[i] = inputind;
	}		
}

void outheap(float *a, int *ind, int *n, int inputind, int *outputind)
{
	
	swapheap(a, ind, (*n - 1), ind[*n], outputind); 

	*n = *n - 1;
}


void indexx(int n, float *a, int *ind) {
	
/* for double-precision array a[1...n] and integer array ind[1...n],
   finds ind[1...n] so that a[ind[1]] ... a[ind[n]] are in
   ascending order.  Adapted from the Collected Algorithms of the
   CACM, #232, HEAPSORT, by J.W.J. Williams.  Adaptation and c code,
   J. R. Thorstensen, 1998 November 29. */

	int i, nremaining, outputind;
	double out;
	
	for(i = 1; i <= n ; i++) ind[i] = i; 

	setheap(a,ind,n);

	nremaining = n;
	while(nremaining > 1) {
		outheap(a, ind, &nremaining, ind[nremaining+1],&outputind);
			ind[nremaining+1] = outputind;  
	}
}

/* Object list section. */

struct objct {
	char name[20];
	double ra;
	double dec;
	double mura;
	double mudec;
	float ep;
	float xtra;  /* mag, whatever */
};

struct objct objs[MAX_OBJECTS];
int nobjects;

int read_obj_list() 

/* Reads a list of objects from a file.  Here's the rules:
     -- each line of file  must be formatted as follows:

name   rahr ramn rasec   decdeg decmin decsec   epoch   [optional number]

    -- the name cannot have any spaces in it.
    -- if the declination is negative, the first character of the
         declination degrees must be a minus sign, and there can
         be no blank space after the minus sign (e. g., -0 18 22 is ok;
	 - 12 13 14 is NOT ok).
    -- other than that, it doesn't matter how you space the lines. 

    I chose this format because it's the standard format for 
    pointing files at my home institution.  If you'd like to code
    another format, please be my guest!  */

{
        FILE *inf;
	char fname[60], resp[10];
	char buf[200];
	char decstr[10];
	double rah, ram, ras, ded, dem, des, ept;
        int i, nitems;
	float xtr1, xtr2, xtr3;

	printf("\nThis reads from a file of objects.  Basic Format:\n\n");
	printf("name_no_blanks<20char   hr mn sec  deg mn sec  epoch \n\n");
	printf(" -- one object per line, blanks between fields, otherwise free-form.\n"); 
	printf(" -- can optionally add [your own floating point number] or\n");
	printf(" --  [mu_ra] [mu_dec] or\n");
	printf(" --  [mu_ra] [mu_dec] [your own floating point number]\n");
	printf(" -- anything NON-NUMERIC after this information is ignored.\n");
        printf("Anything after the optional user-defined floating pt number is ignored.\n");
	printf("Error checking is unsophisticated; maximum of %d objects.\n\n",
		(MAX_OBJECTS-1));
	printf("Give name of file of objects (or QUIT):");
	scanf("%s",fname);
	
	if(strcmp(fname,"QUIT") == 0) {
		printf("Quitting ... did not attempt to open a file.\n");
		return(-1);
	}
	
	inf = fopen(fname,"r");

	if(inf == NULL) {
		printf("\nFILE DID NOT OPEN!  ... I quit.  Try again if you want.\n");
		return(-1);
	}
	else printf("\nopened ... \n");

	if(nobjects != 0) {
		printf("\nYou have %d objects already!\n",nobjects);
	        printf("Type a to append, or r to replace:");
		scanf("%s",resp);
		if(resp[0] == 'r') nobjects = 0;
	}
        /*  on first pass be sure xtra's have a value, just in case. */
        else for(i = 1; i < MAX_OBJECTS; i++) objs[i].xtra = 0.0;

	while((fgets(buf,200,inf) != NULL) && (nobjects < MAX_OBJECTS - 1)) {		
		nobjects++;   /* this will be 1-indexed */
		nitems = sscanf(buf,"%s %lf %lf %lf  %s %lf %lf  %lf %f",
			objs[nobjects].name,&rah,&ram,&ras,decstr,&dem,&des,
				&ept,&xtr1,&xtr2,&xtr3);			
		if(nitems >= 8) {  /* a little error checking here ... */
	             objs[nobjects].ra = rah + ram/60. + ras/3600.;
		     sscanf(decstr,"%lf",&ded);   /* careful with "-0" */
		     if(decstr[0] == '-') {
			if(ded <= 0.) 
				ded = ded * -1.;
			objs[nobjects].dec = 
					-1.* (ded + dem/60. + des/3600.);
		     }
	             else objs[nobjects].dec = ded + dem/60. + des/3600.;
                     objs[nobjects].ep = ept;
		     if(nitems == 9) {
			objs[nobjects].xtra = xtr1;
			objs[nobjects].mura = 0.; objs[nobjects].mudec = 0.;
		     } 	
		     else if(nitems > 9) {
			objs[nobjects].mura = xtr1;
			objs[nobjects].mudec = xtr2;
		        if(nitems == 11) objs[nobjects].xtra = xtr3;
		            else objs[nobjects].xtra = 99.9;
		     }
		     else if(nitems == 8) objs[nobjects].xtra = 99.9;
                }
	        else {
		     printf("Ignoring bad line: %s",buf);
		     nobjects--;
	        }
        }
	printf("\n .... %d objects read from file.\n",nobjects);
	if(nobjects == MAX_OBJECTS - 1)
	  printf("** WARNING ** AT MAX NUMBER OF OBJECTS. You may have missed some.\n");
	fclose(inf);
	return(0);  /* success */ 	
}

int find_by_name(ra, dec, epoch, date, use_dst, enter_ut, night_date, stdz,
		lat, longit)      	
	double *ra, *dec, epoch, stdz, lat, longit;
	struct date_time date;	
	short use_dst, enter_ut, night_date;
                                 
{

	/* finds object by name in list, and sets ra and dec to
           those coords if found.  Precesses to current value of
           epoch. */

	char objname[20];
	int i, found = 1;
	double jd, curep, curra, curdec, sid, ha, alt, az, 	
			par, secz, precra, precdec;

	if(nobjects == 0) {
		printf("No objects!\n");
		return;
	}

#if SYS_CLOCK_OK == 1
	if(update_on == 1) {   
       		if(get_sys_date(&date,use_dst,enter_ut,
			night_date,stdz,update_delta) != 0) 
		       printf("Can't get system date! \n");
	}
#endif

	jd = true_jd(date, use_dst, enter_ut, night_date, stdz);
	
	if(jd < 0.) {
		printf(" ... exiting because of time error!  Try again.\n");
		return(-1);
	}

	sid = lst(jd, longit);  
	curep = 2000. + (jd - J2000) / 365.25;

	printf("RA and dec will be set to list object you name.\n");
	printf("Give exact object name:");
	scanf("%s",objname);		

	i = 1;
	while((i <= nobjects) &&
		((found = strcmp(objs[i].name,objname)) != 0)) i++;
        if(found == 0) {
		if(objs[i].ep != epoch) {
		        cooxform(objs[i].ra,objs[i].dec,objs[i].ep,
                                epoch,&precra,&precdec,XFORM_JUSTPRE,XFORM_FROMSTD);
		}
		else {
			precra = objs[i].ra;
			precdec = objs[i].dec;
		}
		*ra = precra;
		*dec = precdec;
		printf("\nObject found -- name, coords, epoch, user#, HA, airmass --- \n\n");
		printf("%s  ",objs[i].name);
		put_coords(objs[i].ra,3,0);
		printf("  ");
		put_coords(objs[i].dec,2,1);
		printf("  %6.1f  %5.2f ",objs[i].ep,objs[i].xtra);
               	cooxform(objs[i].ra,objs[i].dec,
				objs[i].ep,curep,&curra,&curdec,
				XFORM_JUSTPRE, XFORM_FROMSTD);
    		ha = adj_time(sid - curra);
 		alt=altit(curdec,ha,lat,&az,&par);
		secz = secant_z(alt);
		print_ha_air(ha,secz,0,1); 
		printf("\n\n COORDINATES ARE NOW SET TO THIS OBJECT.\n"); 	
		if(objs[i].ep != epoch) 
		printf("(RA & dec have been precessed to %6.1f, your standard input epoch.)\n",
			epoch);
		return(0);
	}
	else {
		printf("Not found in %d entries.\n",i);
		return(-1);
	}
}

void type_list(date, use_dst, enter_ut, night_date, stdz,
		lat, longit)      	
     	double stdz, lat, longit;
	struct date_time date;	
	short use_dst, enter_ut, night_date;
                                  

{
	int i;
	double jd, curep, curra, curdec, sid, ha, alt, az, secz, par;
	short nstart = 1, nend, ok;
	char errprompt[40];
	
	if(nobjects == 0) {
		printf("No objects!\n");
		return;
	}

#if SYS_CLOCK_OK == 1
	if(update_on == 1) {   
       		if(get_sys_date(&date,use_dst,enter_ut,
			night_date,stdz,update_delta) != 0) 
		       printf("Can't get system date! \n");
	}
#endif
	jd = true_jd(date, use_dst, enter_ut, night_date, stdz);
	
	if(jd < 0.) {
		printf(" ... exiting because of time error!  Try again.\n");
		return;
	}

	sid = lst(jd, longit);  
	curep = 2000. + (jd - J2000) / 365.25;
	
	printf("(Listing will show name, coords, epoch, user#, HA and airmass.)\n");

	strcpy(errprompt,"ERROR IN INPUT ... ");
	oprntf("%d objects in list.\n",nobjects);
/*	while((nobjects > 0) && (nstart > 0)) { used to loop -- too tricky. */
	   /* get out the heavy input checking artillery to avoid
                running away here ... */
	  /*      printf("First and last (numbers) to list, -1 exits:"); */
	        printf("First and last (numbers) to list:");
		ok = getshort(&nstart,-1,(short)nobjects,errprompt);
		/* if(nstart < 0) break; */
		ok = getshort(&nend,nstart,(short)nobjects,errprompt);
		if(nend > nobjects) nend = nobjects;
		if(nstart > nend) nstart = nend;
		oprntf("\n\n");
		print_current(date,night_date,enter_ut);
		oprntf("\n");
                for(i = nstart; i <= nend; i++) {
 			oprntf("%20s ",objs[i].name);
			put_coords(objs[i].ra,3,0);
			oprntf("  ");
			put_coords(objs[i].dec,2,1);
			oprntf("   %6.1f  %5.3f ",objs[i].ep,objs[i].xtra);
               		cooxform(objs[i].ra,objs[i].dec,
				objs[i].ep,curep,&curra,&curdec,
				XFORM_JUSTPRE, XFORM_FROMSTD);
    			ha = adj_time(sid - curra);
			alt=altit(curdec,ha,lat,&az,&par);
			secz = secant_z(alt); 
			print_ha_air(ha,secz,0,1);
			oprntf("\n");	
		}
		oprntf("\n");
/*	}*/	
}

int find_nearest(ra, dec, epoch, date, use_dst, enter_ut, night_date, stdz,
		lat, longit)     
    
	double *ra, *dec, epoch, stdz, lat, longit;
	struct date_time date;	
	short use_dst, enter_ut, night_date;

/* given ra,dec, and epoch, sorts items in list with 
   respect to arc distance away, and queries user
   whether to accept.  */

{
 	int i, ind[MAX_OBJECTS];
	double precra, precdec, jd, sid, ha, alt, az, par,
		secz, seczob, curra, curdec, curep, hacrit, altcrit, aircrit;
	float arcs[MAX_OBJECTS];
	char resp[10];
	int found = 0;
        short sortopt,nprnt;

	if(nobjects == 0) {
		printf("No objects!\n");
		return(-1);
	}
#if SYS_CLOCK_OK == 1
	if(update_on == 1) {   
       		if(get_sys_date(&date,use_dst,enter_ut,
			night_date,stdz,update_delta) != 0) 
		       printf("Can't get system date! \n");
	}
#endif

	jd = true_jd(date, use_dst, enter_ut, night_date, stdz);
		
	if(jd < 0.) {
		printf(" ... exiting because of time error!  Try again.\n");
		return(-1);
	}
	
	printf("Select:\n");
	printf("  1 .. to rank objects by arc dist from current coords;\n");
	printf("  2 .. to rank objects by absolute value of hour angle:\n");
	printf("  3 .. sort by proximity in airmass to present coords:\n");
	printf("  4 .. sort by imminence of setting (reaching critical airmass):\n");
        printf("  5 .. sort by user-supplied optional number.\n");
	printf("  Type 1, 2, 3, 4, or 5 ---> ");
	scanf("%hd",&sortopt);

	if(sortopt == 4) {
		printf("Give critical airmass in west:");
	        scanf("%lf",&aircrit);
		if(aircrit < 1.) {
			oprntf("Airmass must be > 1. ... exiting!\n");
			return(-1);
		}
		altcrit = DEG_IN_RADIAN * asin(1.0 / aircrit);
	}
	sid = lst(jd, longit);  
	curep = 2000. + (jd - J2000) / 365.25;

	 /* compute present airmass for option 3 */
        if(sortopt == 3) {
	    cooxform(*ra,*dec,epoch,
                         curep,&curra,&curdec,XFORM_JUSTPRE,XFORM_FROMSTD); 
	    ha = adj_time(sid - curra);
	    alt=altit(curdec,ha,lat,&az,&par);
	    seczob = secant_z(alt);
        }

	for(i = 1; i <= nobjects; i++) {
		if(sortopt == 1) {   /* sort by arc distance */
		   if(objs[i].ep != epoch) 
 	  		cooxform(objs[i].ra,objs[i].dec,objs[i].ep,
                                epoch,&precra,&precdec,XFORM_JUSTPRE,XFORM_FROMSTD);
 		   else {
			precra = objs[i].ra;
			precdec = objs[i].dec;
		   }
	           arcs[i] = subtend(*ra,*dec,precra,precdec);
		}
		else if (sortopt == 2) {  /* sort by hour angle */
		   cooxform(objs[i].ra,objs[i].dec,objs[i].ep,
				curep,&curra,&curdec,XFORM_JUSTPRE,XFORM_FROMSTD);
	           arcs[i] = fabs(sid - curra);  
		}
		else if (sortopt == 3) {  /* sort by difference of airmass */
		   cooxform(objs[i].ra,objs[i].dec,objs[i].ep,
				curep,&curra,&curdec,XFORM_JUSTPRE,XFORM_FROMSTD);
    		   ha = adj_time(sid - curra);
		   alt=altit(curdec,ha,lat,&az,&par);
	    	   arcs[i] = fabs(secant_z(alt) - seczob);
		}
                else if (sortopt == 4) {  /* sort by proximity to critical airmass */
		   cooxform(objs[i].ra,objs[i].dec,objs[i].ep,
				curep,&curra,&curdec,XFORM_JUSTPRE,XFORM_FROMSTD);
		   hacrit = ha_alt(curdec,lat,altcrit);
                   if(fabs(hacrit) > 24.) arcs[i] = 100.;
    		   else {
	         	ha = adj_time(sid - curra);
                        arcs[i] = hacrit - ha;
		        if(arcs[i] < 0.) arcs[i] = 100.;
		   }
               }
               else   /* sort by user-supplied extra number */
		   arcs[i] = (double) objs[i].xtra;
	}
        
	indexx(nobjects,arcs,ind);

	printf("If you now select an object, RA & dec will be set to its coords.\n\n");
	if(sclogfl != NULL) fprintf(sclogfl,"\n\n *** Sorted object listing *** \n");
	oprntf("Listed for each: Name, ra, dec, epoch, user-defined #,\n");
	if(sortopt == 1) oprntf("arclength to coords (deg), ");
        if(sortopt == 4) oprntf("minutes til crit secz, ");
	oprntf("HA, secz, computed for the following date & time:\n\n");
	print_current(date,night_date,enter_ut);
	oprntf("\n\n");
	if((sortopt != 2) && (sortopt != 4)) { /* print relevant info */
	    oprntf("Current coords:  ");
	    put_coords(*ra,3,0);
	    oprntf("  ");
            put_coords(*dec,2,1);
	    oprntf("  %6.1f ",epoch);
	    cooxform(*ra,*dec,epoch,
                                curep,&curra,&curdec,XFORM_JUSTPRE,XFORM_FROMSTD); 
	    ha = adj_time(sid - curra);
	    alt=altit(curdec,ha,lat,&az,&par);
	    seczob = secant_z(alt);
	    print_ha_air(ha,seczob,0,1);
	    oprntf("\n\n");
	}

	i = 1;
	while(found == 0) {
	    for(nprnt=1;nprnt<=10;nprnt++) {		     
		cooxform(objs[ind[i]].ra,objs[ind[i]].dec,
				objs[ind[i]].ep,
                                   curep,&curra,&curdec,XFORM_JUSTPRE,XFORM_FROMSTD);
    		ha = adj_time(sid - curra);
		alt=altit(curdec,ha,lat,&az,&par);
		
 		oprntf("%2d %13s ",i,objs[ind[i]].name);
		put_coords(objs[ind[i]].ra,3,0);
		oprntf("  ");
		put_coords(objs[ind[i]].dec,2,1);
		oprntf("  %6.1f %6.2f ",objs[ind[i]].ep, objs[ind[i]].xtra);
		if(sortopt == 1) oprntf(" %6.3f",arcs[ind[i]] * DEG_IN_RADIAN);
                if(sortopt == 4) oprntf(" %5.0f",arcs[ind[i]] * 60.); 
		secz = secant_z(alt);
		print_ha_air(ha,secz,0,1);
                oprntf("\n");
                if(nprnt == 5) oprntf("\n");
		i++;
	        if(i > nobjects) break;
            }
	    printf("Type number to select an object, m to see more, q to quit:");

 	    scanf("%s",resp);
 	    if(resp[0] == 'q') {
			oprntf("Abandoning search.\n");
			return(found = -1);
	    }
	    else if((resp[0] == 'm') || (resp[0] == 'M')) {
		if(i > nobjects) {
			oprntf("Sorry -- that's all you have!\n");
			oprntf("Search abandoned.\n");
	        	return(found = -1);
                }
	    }
	    else if(isdigit(resp[0]) != 0) {
		sscanf(resp,"%d",&i);
		if((i < 0) || (i > nobjects)) {
		     	oprntf("BAD OBJECT INDEX -- %d -- start over!\n",i);
			return(-1);
		}
 		if(objs[ind[i]].ep != epoch) 
			     cooxform(objs[ind[i]].ra,objs[ind[i]].dec,
				objs[ind[i]].ep,
                                   epoch,&precra,&precdec,XFORM_JUSTPRE,XFORM_FROMSTD);
		else {
			precra = objs[ind[i]].ra;
			precdec = objs[ind[i]].dec;
		}
                *ra = precra;
		*dec = precdec;
		oprntf("\n%s  ",objs[ind[i]].name);
		put_coords(objs[ind[i]].ra,3,0);
		oprntf("  ");
		put_coords(objs[ind[i]].dec,2,1);
		oprntf("  %6.1f  %5.2f ",objs[ind[i]].ep,objs[ind[i]].xtra);
               	cooxform(objs[ind[i]].ra,objs[ind[i]].dec,
				objs[ind[i]].ep,curep,&curra,&curdec,XFORM_JUSTPRE,
				XFORM_FROMSTD);
    		ha = adj_time(sid - curra);
 		alt=altit(curdec,ha,lat,&az,&par);
		secz = secant_z(alt);
		print_ha_air(ha,secz,0,1); 
		oprntf("\n\n COORDINATES ARE NOW SET TO THIS OBJECT.\n"); 	
		if(objs[ind[i]].ep != epoch) 
   		oprntf("(RA & dec have been precessed to %6.1f, your current standard epoch.)\n",
				epoch);
		return(found = 1);
	    }
	    else {
		printf("Unrecognized response ... continuing ..\n");
		if(i > nobjects) {
			printf("That's all the objects .. abandoning search.\n");
			return(found = -1);
		}
            }
       }
}		

void set_zenith(date, use_dst, enter_ut, night_date, stdz, lat, 
	  longit, epoch, ra, dec) 
               
struct date_time date; 
short use_dst, enter_ut, night_date;
double stdz, lat, longit, epoch, *ra, *dec; 
               
/* sets RA and dec to zenith as defined by present time and date;
   coords are set to actual zenith but in currently used epoch.  */

{
	double jd, curep;

#if SYS_CLOCK_OK == 1
	if(update_on == 1) {   
       		if(get_sys_date(&date,use_dst,enter_ut,
			night_date,stdz,update_delta) != 0) 
		       printf("Can't get system date! \n");
	}
#endif
	
	jd = true_jd(date, use_dst, enter_ut, night_date, stdz); 

	if(jd < 0.) return;  /* nonexistent time. */

	*ra = lst(jd, longit); 

        *dec = lat;

	curep = 2000. + (jd - J2000) / 365.25;
	
	cooxform(*ra,*dec,epoch,curep,ra,dec,XFORM_JUSTPRE,XFORM_TOSTDEP);

	oprntf("\nCOORDINATES SET TO ZENITH: ");
	put_coords(*ra,3,0);
	oprntf(", ");
	put_coords(*dec,2,1);
	oprntf(", epoch %6.1f\n",epoch);
}



void print1phase(date, use_dst, enter_ut, night_date, stdz, lat, 
	  longit, epoch, ra, dec) 

/* prints phase of a repeating phenomenon at this instant. */

struct date_time date; 
short use_dst, enter_ut, night_date;
double stdz, lat, longit, epoch, ra, dec;

{
	double ans, phase, phaserr, jd, hjd, tcor, vcor, curep, curra, curdec;

	printf("Current period %f +- %f, epoch %f +- %f\n",
		star_period, star_perr, star_tzero, star_terr);
	printf("Give new period (d), or negative to keep old:");
	getdouble(&ans,-100.,1e6,
		"Answer out of range ... try again:");
	if(ans > 1.0e-5) {
		star_period = ans;
		printf("Give error bar for period:");	
		getdouble(&star_perr,0.,10000.);
	}
	else if(ans > 0. && ans < 1.0e-5) 
	 printf("That period is too short ... retained the old one!\n");
	printf("Give new epoch (full jd), or negative to keep old:");
	getdouble(&ans,-100.,LASTJD,
		"Answer out of range ... try again:");
	if(ans > 1.0e-5) {
		star_tzero = ans;
		printf("Give error bar for epoch:");
		getdouble(&star_terr,0.,10000.);
	}
	else if(ans > 0. && ans < FIRSTJD) 
          printf("That epoch is too early ... retained the old one!\n");
		
	jd = true_jd(date, use_dst, enter_ut, night_date, stdz); 

	if(jd < 0.) {
		printf("nonexistent time.\n");
	        return;  
	}

 	oprntf("\n\nRA & dec: ");
  	put_coords(ra,3,0);
  	oprntf(", ");
  	put_coords(dec,2,1);
  	oprntf(", epoch %6.1f\n", epoch);
	oprntf("Ephemeris: %lf (+- %lf) + \n     %14.9f (+- %12.9f) E (Heliocentric.\n",
		star_tzero, star_terr, star_period, star_perr);
	oprntf("Julian date %f --> UT ",jd);
	print_all(jd);
	oprntf("\n\n");
	
	curep = 2000. + (jd - J2000) / 365.25;
	
	cooxform(ra,dec,epoch,curep,&curra,&curdec,XFORM_JUSTPRE,XFORM_FROMSTD);

	helcor(jd,curra,curdec,0.,0.,0.,&tcor,&vcor);
		/* throwing away vcor .... tcor doesn't need ha, lat,
			or elevation. */
	hjd = jd + tcor/86400.;
	
	phase = (jd - star_tzero) / star_period; 

	/* timing error is error in period times number of cycles
		since epoch added in quadrature with error of epoch ... */
	phaserr = (phase * star_perr * phase * star_perr) 
		+ (star_terr * star_terr) ; 
	phaserr = sqrt(phaserr / (star_period * star_period));
	
	oprntf("If input time already helio, phase is %f\n",phase);
	if(phase < 0.) {
		phase = phase - (int) phase + 1;
		oprntf(" ... i.e., %f\n",phase);
	}
	phase = (hjd - star_tzero) / star_period;
	oprntf("Making helio correction for current coords, phase is %f\n",
		phase);
	if(phase < 0.) {
		phase = phase - (int) phase + 1;
		oprntf(" ... i.e., %f\n",phase);
	}
	printf("Phase uncertainty:  %10.3f\n",phaserr);
	if(phaserr > 0.3) printf("VERY UNCERTAIN PHASE AT THIS EPOCH.\n");

	oprntf("(Helio correction is %f seconds, helio J.D. is %lf)\n",
		tcor,hjd);
			
} 

int set_to_jd(date, use_dst, enter_ut, night_date, stdz, jd, verbose) 

	struct date_time *date;
	short use_dst, enter_ut, night_date;
        double stdz, jd;
        int verbose;  /* 1 = prints, 0 = silent */
{
	/* Takes a jd and loads up the date structure
           to conform to the prevailing conventions for the interpretation
           of times. */

        short dow;
        double jdb, jde, jdloc;

	if((jd < FIRSTJD) || (jd > LASTJD)) {
		oprntf("JD out of range .. allowed %f to %f\n",
			FIRSTJD,LASTJD);
		oprntf("(Corresponds to 1901 --> 2099 calendrical limits.\n");
		oprntf("No action taken ... try again if you want.\n");
		return(-1);
	}

	caldat(jd,date,&dow);   /* first get year (at least) */

        if(enter_ut != 1)  { /* adjust if needed */
           find_dst_bounds(date->y,stdz,use_dst,&jdb,&jde);
	   jd = date_to_jd(*date);
	   jdloc = jd - zone(use_dst,stdz,jd,jdb,jde)/24.;
	   if(use_dst > 0) {  /* north ... daylight savings */
	      if((jde - jd) < 0.041666667 && (jde - jd) > 0.) {
 		oprntf("\n\n IN AMBIGUOUS DST/STD HOUR -- TIME WILL NOT BE CORRECT \n");
	        oprntf("\n   ... use 'g' to set UT input, then do over.\n\n");
	      }
	   }
	   else if(use_dst < 0) {  /* south ... daylight savings */
	      if((jdb - jd) < 0.041666667 && (jdb - jd) > 0.) {
 		oprntf("\n\n IN AMBIGUOUS DST/STD HOUR -- TIME WILL NOT BE CORRECT \n");
	        oprntf("\n   ... use 'g' to set UT input, then do over.\n\n");
	      }
	   }
           caldat(jdloc,date,&dow);
        }
	if(verbose > 0.) {
		oprntf("\nDATE AND TIME SET FROM JD ... to\n");
		print_all(date_to_jd(*date));
		if(enter_ut == 0) oprntf(" local time\n");
        		else oprntf(" Universal time\n"); 
      		if((night_date == 1) && (date->h < 12)) {
			date->d = date->d - 1;
			oprntf("(night-date is on, so internal value of date set to last evening)\n"); 		   
		}
	}
        return(0); /* success */
}

void phaselisting(use_dst, enter_ut, night_date, stdz, lat, 
	  longit, epoch, ra, dec) 

/* prints phase of a repeating phenomenon at this instant. */

short use_dst, enter_ut, night_date;
double stdz, lat, longit, epoch, ra, dec;

{
	double ans, phase, phaserr, hjd, tcor, vcor, curep, curra, curdec;
	double jd, jdstart, jdend, interval;
	struct date_time tempdate;
  	double max_secz, min_alt, max_alt, min_ok_alt, sun_alt, 
		max_ok_sun, secz, par;
  	double sid, ha, altitude, az, rasun, decsun, hasun;
        double ramoon,decmoon,distmoon,moon_alt;
	int valid_date = 0, i = 0, ncycle, ntotry;
        int  nskipped = 0, nprinted = 0, nbanner = 0;
	short dow;

	printf("Current period %f +- %f, epoch %f +- %f\n",
		star_period, star_perr, star_tzero, star_terr);
	printf("Give new period (d), or negative to keep old:");
	getdouble(&ans,-100.,1e6,
		"Answer out of range ... try again:");
	if(ans > 1.0e-5) {
		star_period = ans;
		printf("Give error bar for period:");	
		getdouble(&star_perr,0.,10000.);
	}
	else if(ans > 0. && ans < 1.0e-5) 
	 printf("That period is too short ... retained the old one!\n");
	printf("Give new epoch (full jd), or negative to keep old:");
	getdouble(&ans,-100.,LASTJD,
		"Answer out of range ... try again:");
	if(ans > 1.0e-5) {
		star_tzero = ans;
		printf("Give error bar for epoch:");
		getdouble(&star_terr,0.,10000.);
	}
	else if(ans > 0. && ans < FIRSTJD) 
          printf("That epoch is too early ... retained the old one!\n");

       printf("(Give starting and ending times in current date and time convention ...\n");
       if(enter_ut == 1) printf("Times & dates should be UT, ");
       else {
		printf("Times and dates should be local, zone = %4.1f hr W\n");
	        if(use_dst != 0) printf("Daylight savings applies if current,\n");
                if(night_date == 1) printf("** NIGHT-DATE on.\n");
       } 
       printf("Give starting date for listing, yyyy mm dd, neg year exits:");
       if((valid_date = get_date(&tempdate)) < 0) return;
       printf("Give starting time, h m s:");
       get_time(&tempdate);  
       jdstart = true_jd(tempdate, use_dst, enter_ut, night_date, stdz); 
       printf("Give ending date for listing:");
       get_date(&tempdate);
       printf("Give ending time:");
       get_time(&tempdate);  
       jdend = true_jd(tempdate, use_dst, enter_ut, night_date, stdz); 
       if(jdend < jdstart) {
		printf("Sorry ... end after beginning!\n");
		return;
	}
       printf("Give time interval of phase points in minutes:");
       scanf("%lf",&interval);
       interval /= 1440.;  /* in days ... */
       ntotry = (int) ((jdend - jdstart) / interval);
       printf("Give maximum airmass for listing, negative to print all:");
       scanf("%lf",&max_secz);
       if(max_secz > 0.) min_ok_alt = 90. - DEG_IN_RADIAN * acos(1/max_secz);
       else min_ok_alt = -100.;
       curep = 2000. + (jdstart - J2000)/365.25;    
           /* precess to check observ. */
       cooxform(ra,dec,epoch,curep,&curra,&curdec,
		XFORM_JUSTPRE,XFORM_FROMSTD);
       min_max_alt(lat,curdec,&min_alt,&max_alt);
       if(max_alt < min_ok_alt) {
          oprntf("Object never rises so high -- max alt %5.1f deg -- returning!\n",
           max_alt);
	  return;
       }

       printf("Give maximum allowable sun alt to print (e.g. -12), 90 to print all:");
       scanf("%lf",&max_ok_sun);

		
 	oprntf("\n\nRA & dec: ");
  	put_coords(ra,3,0);
  	oprntf(", ");
  	put_coords(dec,2,1);
  	oprntf(", epoch %6.1f\n", epoch);
	oprntf("Ephemeris: %lf (+- %lf) + \n     %14.9f (+- %12.9f) E (Heliocentric.)\n",
		star_tzero, star_terr, star_period, star_perr);
        oprntf("Periodic listing of phases in interval from ");
  	print_all(jdstart,&dow);
  	oprntf(" UT to ");
  	print_all(jdend,&dow);
  	oprntf(" UT.\n");
  	if(min_ok_alt > 0.) oprntf("Only events at sec(z) < %5.2f are shown.\n",
            max_secz);
  	if(max_ok_sun < 89.) 
	    oprntf("Only events with sun altit. < %6.1f degrees are shown.\n\n",
		max_ok_sun);
       curep = 2000. + (jdstart - J2000) / 365.25;
       cooxform(ra,dec,epoch,curep,&curra,&curdec,XFORM_JUSTPRE,XFORM_FROMSTD);
	/* assume coordinates sensibly constant during calculation ... */
  /*      oprntf("..date..(UT)..time     N    phase    err    ha   secz   sunalt  moonalt\n");
*/
	for(i = 0; i <= ntotry; i++) {
          jd = jdstart + (double) i * interval;  
		/* avoids roundoff accumulation */
          sid = lst(jd,longit);
          ha = adj_time(sid - curra);
          altitude = altit(curdec,ha,lat,&az,&par);
	  secz = secant_z(altitude);
          if(altitude >= min_ok_alt) {
             lpsun(jd,&rasun,&decsun);  /* lpsun plenty good enough */
	     hasun = sid - rasun;
	     sun_alt = altit(decsun,hasun,lat,&az,&par);
             if(sun_alt < max_ok_sun) {
                if(nskipped > 0 || nprinted == 0) {
			if(nbanner > 15 || nprinted == 0)  {
                                oprntf("\n...date..(UT)..time     N    phase");
                                oprntf("     err       ha   secz  sunalt ");
                                oprntf("moonalt\n\n");
                                nbanner = 0;
                        }
                        else oprntf("\n");
		}
                nskipped = 0;
	        print_calendar(jd,&dow);
                oprntf(",");
		print_time(jd,0);
 	        helcor(jd,curra,curdec,0.,0.,0.,&tcor,&vcor);
	         	/* throwing away vcor .... tcor doesn't need ha, lat,
			or elevation. */
	        hjd = jd + tcor/86400.;
	        phase = (hjd - star_tzero) / star_period;
                ncycle = (int) (floor(phase));
	  /* timing error is error in period times number of cycles
		since epoch added in quadrature with error of epoch ... */
	  	phaserr = (phase * star_perr * phase * star_perr) 
		    + (star_terr * star_terr) ; 
	        phaserr = sqrt(phaserr / (star_period * star_period));
                phase = phase - floor(phase);
	        oprntf(" %7d   %5.3f ",ncycle,phase);
	        oprntf(" +- %5.3f  ",phaserr);
		put_coords(ha,0,1);
                oprntf(" ");
		if(fabs(secz) < 100.) 
          	   oprntf("%5.2f ",secz);
		else oprntf("v.low  ");
		if(sun_alt > -18.) oprntf("  %4.0f ",sun_alt);
	        else oprntf("  .... ");
		lpmoon(jd,lat,sid,&ramoon,&decmoon,&distmoon); /* close enuf */
		moon_alt=altit(decmoon,(sid-ramoon),lat,&az,&par);
		if(moon_alt < -2.)  oprntf("  ....\n");
                else oprntf("  %4.0f\n",moon_alt);
                nbanner++;
                nprinted++;
             }
	     else nskipped++;
	  }  
	}			
} 

main()

{
	struct date_time date,dateback;
	struct coord ttime; 
	double jd, jdmid, jdc, jdtmp;
	double jdb, jde, test;  /* jd of begin and end of dst */
	double sid, sss;
	double Dt; /* ephemeris time correction */
	short option, trying;
	short done = 0, optdone = 0, valid_date = 0, nreturns=0, nxret = 0;
	short day, yearday, dst=0, dow;
	int cc, end_in;  /* control character for circumstances loop */
	int cx;    /* control character for extra goodies ... */
        char cxstr[3];
	double objra=0., objdec=0., objepoch=2000.,dectest, ep;     
	char constelname[5];   /* 3-letter constellation abbreviation. */
	double curra, curdec, curep, obj_moon;
	double pra[10],pdec[10];
	double alt, az, ha, secz, jdloc;
	double rasun,decsun,ramoon,decmoon,distmoon;
	short enter_ut = 0; /* are times to be entered as UT? */
	short night_date = 1; /* interperet current date as evening or true? */
	char str[80]; /* dummy string */
	char errprompt[80];
	short nch;
	double glong, glat, eclong, eclat, par;
	int status;
	double mura_sec=0.,mura_arcs=0.,mudec=0.,mura,mudc;  /* proper motions */
	double objra_adj, objdec_adj;           /* equinox of std epoch,
			adjusted for proper motion only */
	short hr_span,i;  /* for table of hour angles */
	double tcor, vcor, vcorlsr; /* time and velocity helio corrections */
	char obs_code;
	double minoffset; /* minutes offset from system clock ... */
        double aberra, aberdec, parra, pardec;  /* aberration and 
		parallax factors ... */

	/* all the site-specific quantities are here:
		longit     = W longitude in decimal hours
		lat        = N latitude in decimal degrees
		stdz       = standard time zone offset, hours
		use_dst    = 1 for USA DST, 2 for Spanish, negative for south,
				 0 to use standard time year round
		zone_name  = name of time zone, e. g. Eastern
		zabr       = single-character abbreviation of time zone
		site_name  = name of site.
	        elevsea    = elevation above sea level, meters
		elev       = elevation above local horizon, meters
	*/

	/* Kitt peak, MDM observatory, is initialized here as a default.
	   User later gets to choose from a menu of possible sites -
	   they're all hard-coded in the routine 'load_site'. */

	char site_name[45];  /* initialized later with
			   strcpy for portability */
	char zabr = 'M';
	char zone_name[25]; /* this too */
	short use_dst = 0;
	double longit = 7.44111;
	double elevsea = 1925.;  /* for MDM, strictly */
	double elev = 500.; /* well, sorta -- height above horizon */
	double horiz = 0.7174;
	double lat = 31.9533;
	double stdz = 7.;
	
	/* and off we go.... */

	strcpy(site_name,"Kitt Peak");
	strcpy(zone_name,"Mountain");

	star_tzero = 0.; 
	star_period = 0.;
	star_terr = 0.;
	star_perr = 0.;  
	
	oprntf("\nAstronomical calculator program, by John Thorstensen.\n\n");

	load_site(&longit,&lat,&stdz,&use_dst,zone_name,&zabr,
			&elevsea,&elev,&horiz,site_name);
	oprntf("You have selected %s\n",site_name);
        oprntf("(This can be changed later using the 's' [site] command.)\n\n");

#if SYS_CLOCK_OK == 1

        if(get_sys_date(&date,use_dst,enter_ut,night_date,stdz,0.) != 0) {
	  date.y = 2000;  /* have to have a default date.*/
	  date.mo = 1;
	  date.d = 1;
	  date.h = 0.;
	  date.mn = 0.;
	  date.s = 0.;        
	  oprntf("SYSTEM CLOCK didn't read. Time & date set arbitrarily to\n");
	  print_all(date_to_jd(date));
	  oprntf("\n");
        }

        else set_zenith(date, use_dst, enter_ut, night_date, stdz, lat, 
	  longit, objepoch, &objra, &objdec); 

#else 
       	  date.y = 2000;  /* have to have a default date.*/
	  date.mo = 1;
	  date.d = 1;
	  date.h = 0.;
	  date.mn = 0.;
	  date.s = 0.;        
	  oprntf("SYSTEM CLOCK options turned off, so \n ");
	  oprntf("time and date set arbitrarily to:\n");
	  print_all(date_to_jd(date));
	  oprntf("\n\n");
#endif 

        oprntf("\nREADY TO COMPUTE.  Use simple commands to set the date, time, RA\n");
        oprntf("dec, and so on; then use other commands to compute and display\n");
        oprntf("observability information.\n\n");
        oprntf("NEW or RUSTY USERS: type 'f' (and return) for FAST guided tour.\n");

	while((cc = getch()) != 'Q')    switch(cc) {
		case '?':    /* print a menu */
			print_menu();
			nreturns=0;
			break;
                case 'f':    /* print a short tutorial */
                        print_tutorial();
			nreturns=0;
			break;
		case 'r':   /* enter the object's right ascension */
			objra = get_coord();
			nreturns=0;         
			break;
		case 'd':   /* enter the object's declination */ 
			/* filter declination in put for 'date' input! */
			dectest = get_coord();
			if(fabs(dectest) <= 90.) { 
				objdec = dectest;
			}
			else {
			    oprntf("REJECTED 'd' INPUT - DECLINATION MUST BE < 90.\n");
			    oprntf("if you want DATE, Use 'y' (yyyy mm dd)\n");
			}
			nreturns=0;         
			break;

		case 'C':  /* enter RA and dec together ... */
			 objra = get_coord();
	                        /* filter declination in put for 'date' input! */
                        dectest = get_coord();
                        if(fabs(dectest) <= 90.) {
                                objdec = dectest;
                        }
                        else {
                            oprntf("REJECTED 'd' INPUT - DECLINATION MUST BE < 90.\n");
                            oprntf("if you want DATE, Use 'y' (yyyy mm dd)\n");
                        }
                        nreturns=0;
                        break;

		case 'p':   /* enter the object's proper motions -- */

                	status = get_pm(objdec,&mura_sec,&mudec);
			nreturns=0;
			break;  
		case 'e':   /* enter the input epoch */
			getdouble(&objepoch,-10000.,10000.,
				"Give input epoch again...\n"); /* liberal lims*/
                        if(objepoch < -5000.) {
				objepoch = 2000.+
	                            (true_jd(date, use_dst, enter_ut, night_date, stdz) 
				     -J2000)/365.25;
				printf("LARGE NEGATIVE EPOCH --- causes input epoch to be set to current!\n");
				printf("set to Julian epoch %9.4f\n",objepoch);
                        }
			break;
		case 't':   /* enter the time ... hours min sec */
			get_time(&date);
			nreturns=0;         
#if SYS_CLOCK_OK == 1
			if(update_on == 1) {
				update_on = 0;
				printf("Automatic updating DISABLED.\n");
			}
#endif
			break;
		case 'T':   /* read system clock -- set date & time to that. */
			;
#if SYS_CLOCK_OK == 1
			printf("Set to how many minutes into the future? :");
			scanf("%lf",&minoffset);
			get_sys_date(&date,use_dst,enter_ut,night_date,
				stdz,minoffset);
			if(update_on == 1) {
				update_on = 0;
				printf("Automatic updating DISABLED.\n");
			}
#else
			printf("Sorry -- system clock options are disabled, probably because of an\n");
			printf("incompatibility between your system and the standard time library\n");
			printf("functions used in the program.  Turning them on would require fixing\n");
			printf("the problem in source code and recompiling.\n");
#endif	               
			nreturns=0;
                        break;
		case 'g':   /* toggle whether times are entered as
				      Greenwich or local */
			if(enter_ut == 1) enter_ut = 0;
			else {
				enter_ut = 1;
				night_date = 0;
			}
			if(enter_ut == 1) 
				oprntf("Dates and times entered are now UT.\n");
else oprntf("Dates & times entered are local, dates are literal (not evening).\n");
			oprntf("TIME IS CHANGED to %d %02d %02d, %02d %02d %02.0f",
				date.y,date.mo,date.d,date.h,date.mn,date.s); 
			if(enter_ut == 1) oprntf(" UNIVERSAL time.\n");
				else oprntf(" LOCAL time.\n");
			nreturns=0;         
			break;
		case 'n': /* toggle whether the current date is to
				  be interpreted as the evening date (for 
				  all night) or the true date .... */
			if(enter_ut == 1) {
oprntf("You're entering times as UT, so 'evening date' makes no sense!....\n"); 
oprntf("No action taken on 'n', first use 'g' first to enable local time input.\n");
			}
			else if(night_date == 1) {
			   night_date = 0;
oprntf("The date in effect will now be interpreted literally, not as evening.\n");
			}
			else {
			   night_date = 1;
oprntf("The date in effect will now be interpreted as the evening date.\n");
			}
			nreturns=0;         
			break;
		case 'a':   
			if(sclogfl != NULL) fprintf(sclogfl,"\n\n"); 
					/* space it */ 
			oprntf("*** Almanac for the currently specified date ***");
			if(night_date != 1) {
oprntf(", but CAUTION!!\nThe 'night date' option is off, so be especially careful\n");
oprntf("this is the correct night.  See 'g' and 'n'...");
			}
			else oprntf(":");
			/* require a trap here to catch bad dates .... */
			if(date.y < 1901 || date.y > 2099) 
				oprntf("Bad date!  Trapped! \n");
			else {
			    print_tonight(date,lat,longit,elevsea,elev,horiz,site_name,stdz,
			       zone_name,zabr,use_dst,&jdb,&jde,2);
			    printf("\nType command, 'f' for fast tour, or '?' for menu:");
			}
			nreturns=0;         
			break;
		case 'y':   /* enter the date, yyyy mm dd */
			get_date(&date);
			nreturns=0;         
#if SYS_CLOCK_OK == 1
			if(update_on == 1) {
				update_on = 0;
				printf("Automatic updating DISABLED.\n");
			}
#endif
			break;

	/* The site parameters must all be changed at the same time; hence
	   user is forced to change them all. */

		case 's':  /* change the site parameters */
			load_site(&longit,&lat,&stdz,&use_dst,
			   zone_name,&zabr,&elevsea,&elev,
			   &horiz,site_name);
			oprntf("New site = %s\n",site_name);
			printf("(Give command, or ? for menu.)\n");
			nreturns=0;             
			break;
		case '=':  /* PRINT CIRCUMSTANCES for current params */
#if SYS_CLOCK_OK == 1
			if(update_on == 1) {
        			if(get_sys_date(&date,use_dst,enter_ut,
					night_date,stdz,update_delta) != 0) 
				       printf("Can't get system date! \n");
			}
#endif
			if(sclogfl != NULL) fprintf(sclogfl,"\n\n*** Instantaneous Circumstances ***\n");
			if(setup_time_place(date,longit,lat,stdz,
			    use_dst,zone_name,zabr, site_name,enter_ut,
			    night_date,&jd,&jdloc,&jdb,&jde,&sid,
			    &curep) < 0) break;
			print_circumstances(objra,objdec,objepoch,jd,
			    curep,mura_arcs,mura_sec,mudec,
				   sid,lat,elevsea,horiz);
			nreturns=0;         
			break;
		case 'm':  /* print positions of major planets */
#if SYS_CLOCK_OK == 1
			if(update_on == 1) {
        			if(get_sys_date(&date,use_dst,enter_ut,
					night_date,stdz,update_delta) != 0) 
				       printf("Can't get system date! \n");
			}
#endif
			if(setup_time_place(date,longit,lat,stdz,
				use_dst,zone_name,zabr, site_name,enter_ut,night_date,
				&jd,&jdloc,&jdb,&jde,&sid,&curep) < 0)
				   break;
			comp_el(jd);
			pposns(jd,lat,sid,1,pra,pdec);
			nreturns=0;         
			break;                                  
		case 'h':  /* print an hourly airmass table */
#if SYS_CLOCK_OK == 1
			if(update_on == 1) {
        			if(get_sys_date(&date,use_dst,enter_ut,
					night_date,stdz,update_delta) != 0) 
				       printf("Can't get system date! \n");
			}
#endif
			hourly_airmass(date,stdz,lat,longit,horiz,
			   use_dst,objra,objdec,objepoch, mura_sec,
			   mura_arcs,mudec);
			nreturns=0;         
			break;
		case 'o':
			if(sclogfl != NULL) fprintf(sclogfl,"\n");
		        obs_season(objra,objdec,objepoch,
			     lat,longit);
			nreturns=0;
                        break;
		case 'c':  /* print galactic and ecliptic coordinates */
			galact(objra,objdec,objepoch,&glat,&glong);
			oprntf("Galactic: l = %5.2f, b = %5.2f\n",
				glat,glong);
			eclipt(objra,objdec,objepoch,date_to_jd(date),
				&curep,&eclong,&eclat);
			oprntf("Ecliptic (equinox %7.2f): long = %5.2f, lat = %5.2f\n",
				curep, eclong, eclat);
                        radec_to_constel(objra, objdec, objepoch, constelname);
			oprntf("Located in contellation %s.\n",constelname);

			nreturns=0;
			break;
		case 'l': 
			print_params(date,enter_ut,night_date,
				stdz,lat,longit,site_name,elevsea,elev,use_dst,
				objra,objdec,objepoch,mura_sec,mura_arcs,mudec);
			nreturns=0;         
			break;
		case 'i':  /* print a short tutorial */
			print_examples();
			nreturns=0;         
			break;
		case 'w':  /* print information about algorithms, acc. */
			print_accuracy();
			nreturns=0;             
			break;
		case 'x':
			nxret = 0;
    /*			printf("(Give xtra goodies subcommand, ? for menu)\n");
			scanf("%s",cxstr);    */
		        while(isspace(cx = getch()) != 0) {
			    nxret++;
			    if(nxret == 3) {
				printf("Give an extra goodies command, or ? for menu!\n"); 
				nxret = 0;
			    }
			}
            /*          cx = cxstr[0];   */
			switch(cx)  {
                           case '?':
			oprntf("Extra goodies commands are:\n");
			oprntf("  x? ... print extra goodies menu.\n");
			oprntf("  xy ... print day of year.\n");
			oprntf("  xc ... give galactic and ecliptic coords & lsr corr'n.\n");
			oprntf("  xa ... print apparent place (nutation and aberration corrected)\n");
			oprntf("  xd ... give rough value of delta T = TDT - UT.\n");
                        oprntf("  xv ... list geocentric times of repeating phenom (Variable star)\n");
			oprntf("  xf ... give phase of repeating phenom.\n");
			oprntf("  xp ... parallax factors and aberration.\n");
                        oprntf("  xb ... precess a bunch of coords the same way.\n");
			oprntf("  xj ... calculate calendar dates given julian dates.\n");
			oprntf("  xJ ... *set* date and time values from Julian date.\n");
			oprntf("  xZ ... *set* RA and dec to Zenith\n");
#if SYS_CLOCK_OK == 1
			oprntf("  xU ... *toggle* automatic system-clock update.\n");
#endif

#if LOG_FILES_OK == 1
			oprntf("LOG-FILE COMMAND:\n");
			oprntf("  xL ... toggles log file open or closed\n");
#endif
			oprntf("COMMANDS FOR FILES OF OBJECTS:\n");
			oprntf("  xR ... read objects from a file, format: name h m s d m s epoch\n");
			oprntf("  xl ... type out (part of) object list.\n");
			oprntf("  xN ... find object by name, set to its coords\n");
			oprntf("  xS ... sort and select object by a rank, set to coords.\n");
			oprntf("(Note that capital letters affect more than one quantity, eg. RA and dec)");
			   break;

                           case 'v':
			      if(sclogfl != NULL) 
				fprintf(sclogfl,"\n\n  *** Ephemeris predictions ***\n\n");
			      ephemgen(objra,objdec,objepoch,lat,longit);
			   break;

			   case 'f':
#if SYS_CLOCK_OK == 1
			      if(update_on == 1) {
        			    if(get_sys_date(&date,use_dst,enter_ut,
					night_date,stdz,update_delta) != 0) 
				       printf("Can't get system date! \n");
			      }
#endif
/*                              print1phase(date, use_dst, enter_ut, 
				night_date, stdz, lat, 
	  			longit, objepoch, objra, objdec); 
*/
				phaselisting(use_dst, enter_ut, night_date, 						stdz, lat, longit, 
					objepoch, objra, objdec);
			   break;	

			   case 'b':
                              mass_precess();
			   break;

			   case 'a':
#if SYS_CLOCK_OK == 1
			      if(update_on == 1) {
        			    if(get_sys_date(&date,use_dst,enter_ut,
					night_date,stdz,update_delta) != 0) 
				       printf("Can't get system date! \n");
			      }
#endif
			      jd = true_jd(date, use_dst, enter_ut, night_date, stdz);
			      print_apparent(objra,objdec,objepoch,mura_sec,mudec,jd,
				lat,longit,elevsea);
			   break;

			   case 'p':
                                 jd = true_jd(date, use_dst, enter_ut, 
					night_date, stdz);
				 parellipse(jd, objra, objdec, objepoch, lat,
					 longit, &parra,
                			 &pardec, &aberra, &aberdec);

			   break;

			   case 'c':
			      oprntf("Equatorial: RA = ");
			      put_coords(objra,3,0);
			      oprntf(", dec = ");
			      put_coords(objdec,2,1);
			      oprntf(" (epoch %6.1f)\n",
                                  objepoch);
			      galact(objra,objdec,objepoch,&glat,&glong);
			      oprntf("  Galactic: l = %5.2f, b = %5.2f\n",
				glat,glong);
			      eclipt(objra,objdec,objepoch,date_to_jd(date),
		           		&curep,&eclong,&eclat);
			      oprntf("Ecliptic (equinox %7.2f): long = %5.2f, lat = %5.2f\n",
				curep, eclong, eclat);
			      lsrcor(objra,objdec,objepoch,&vcorlsr); 
			      oprntf("Rough correction from helio to lsr %5.1f km/s\n",vcorlsr);
			   break;
			case 'y': /* print day of year ... */
#if SYS_CLOCK_OK == 1
				if(update_on == 1) {   
       					if(get_sys_date(&date,use_dst,enter_ut,
						night_date,stdz,update_delta) != 0) 
			       		printf("Can't get system date! \n");
				}
#endif

				jdtmp = true_jd(date, use_dst, enter_ut, night_date, stdz);
				
				printf("\nUT date & time: ");
				print_all(jdtmp);	
				printf("\n");
				printf("UT day of year:%10.5f\n",day_of_year(jdtmp));
				break;

			   case 'd':
			      jd = true_jd(date, use_dst, enter_ut, night_date, stdz);
       	                      Dt = etcorr(jd);
			      oprntf("Delta t = TDT - UT = %5.1f seconds\n",Dt);
			      oprntf("JD %f (UT) --> ",jd);
			      jd = jd + Dt / SEC_IN_DAY;
			      oprntf(" %f (TDT)\n",jd);
			      if(date.y > 2001) 
				oprntf("(Value is an extrapolated guess ... only computable after the fact.)\n");
			      else oprntf("+- 0.5 sec, based on 5-year linear interpolations.\n"); 
		              if(date.y < 1983) oprntf("Before 1983, TDT was preceded by ephemeris time (ET).\n");
			   break;
#if SYS_CLOCK_OK == 1
			   case 'U':
				if(update_on == 0) {
					update_on = 1;
					printf("Automatic update toggled ON.\n");
					printf("Offset (minutes) into future?:");
					scanf("%lf",&update_delta);
				}
			        else {
					update_on = 0;
					printf("Automatic update toggled OFF.\n");
				}
		           break;
#endif

			   case 'j':
			      jdc = 1;
			      oprntf("jd to calendar conversion. \n");
			      while(jdc > 0.) {
				  printf("Give jd to convert, negative to exit:");
				  getdouble(&jdc,-1000000.,3000000.,
	  			    "Give JD to convert, negative to exit");
			   	  oprntf("%f -- > ",jdc);
			          print_all(jdc);
			          oprntf("\n");
			      }
			       oprntf("(Value of date in main program is unaffected.)\n");
	                   break;			      

			   case 'J':
 			      oprntf("Sets date and time from an input JD.\n");
			      printf("Give JD to set to, negative value for no action: ");
				  getdouble(&jdc,-1000000.,LASTJD,
	  			    "Give jd to set to, negative for no action");
			      set_to_jd(&date, use_dst, enter_ut, 
					night_date, stdz, jdc,1); 
#if SYS_CLOCK_OK == 1
				if(update_on == 1) {
				       update_on = 0;
				       printf("Automatic updating DISABLED.\n");
				}
#endif
	                   break;			      
#if LOG_FILES_OK == 1 
			   case 'L':
				if(sclogfl == NULL) {
				    trying = 1;
				    while(sclogfl == NULL && trying == 1) {
				    	printf("Give filename for log file, type NONE to cancel:");
				    	scanf("%s",str);
					if(strcmp(str,"NONE") == 0) 
						trying = 0;
					else {				    
				    		sclogfl = fopen(str,"a");
					}
				    }
				    if(sclogfl != NULL) 
					printf("log file %s is OPEN in append mode.\n",str);
				    else printf("LOG FILE NOT OPENED.\n");
				}
			        else {
				    fclose(sclogfl);
				    sclogfl = NULL;  /* reset it explicitly */
				    printf("Log file has been CLOSED.\n");
				}
				break;
#endif
			   case 'R':
				read_obj_list();
				break;
			   case 'l':
				type_list(date,use_dst,enter_ut,night_date, 
				    stdz,lat,longit);      	
     				break;
			   case 'N':
				find_by_name(&objra,&objdec,objepoch,date, 
				  use_dst,enter_ut,night_date,stdz,lat,longit);
				break;
			   case 'S':
				find_nearest(&objra,&objdec,objepoch,date,
                                  use_dst,enter_ut,night_date,stdz,lat,longit);
				break;
			   case 'Z':
                                set_zenith(date,use_dst,enter_ut,night_date,
 					stdz,lat,longit,objepoch,&objra,
					&objdec);
				break;
			   case 't':  /* test */
				
	test = true_jd(date, use_dst, enter_ut, night_date, stdz);
	printf(".... true_jd gives --> %f\n",test);
	break;
			   case ' ':  ;
			   break;
			   case '\n': ;
    			   break;
			   default: oprntf("Unrecognized character %c ... no action.\n",cx);
                           break;
                        }      /* end of 'extra goodies' menu. */
			nreturns=0;
			printf("\n(eXtra goodies doesn't loop.)\n");
			printf("Back in main commands, 'Q' quits, '?' menu, 'f' fast tour.\n");
			break;
		case '\n': /* ignore carriage returns */        
			nreturns++;  /* but guide the user if they keep 
			       hitting returns .... */
			if(nreturns == 3) {
			  printf("You're repeating carriage returns. There are no prompts.\n");
			  printf("Type 'f' for fast tour, 'i' for instructions, ? for a menu.\n");     
			  nreturns = 0;
			}
			break;
		case ' ':  /* ignore blank spaces */
			break;
		case 'q':  /* prompt if user's trying to quit */
			printf("Type an UPPER CASE Q to quit.\n");
			break;
		default:   /* complain if unrecognizable */
			printf("Unknown command, %c\n",cc);
	}       /* closing switch loop */
	BLUNDER:; /* DUMMY STATEMENT */
	oprntf("Suggestions or comments --> john.thorstensen@dartmouth.edu\n");
	oprntf("Goodbye.\n");
}
