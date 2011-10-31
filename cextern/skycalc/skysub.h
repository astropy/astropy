/* Header file for skycalc c-language routines.

   ANSI function prototypes, plus brief descriptions, of all the 
   routines used in skycalc.  Also includes constant definitions and
   definitions of a few data structures used.  */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
/* #include <stdarg.h> */
#include <string.h>

/* a couple of the system-dependent magic numbers are defined here */

#define SYS_CLOCK_OK 1    /* 1 means ANSI-standard time libraries do work,
   2 means they don't.  This is used by compiler switches in file 5 and
   the main program.  */

#define LOG_FILES_OK 1  /* 1 means that log files are enabled.
			Any other value means they're not.  */

#define MAX_OBJECTS 500
#define MINSHORT -32767   /* min, max short integers and double precision */
#define MAXSHORT 32767
#define MAXDOUBLE 1.0e38
#define MINDOUBLE -1.0e38
#define BUFSIZE 150 

#define XFORM_FROMSTD  1  /* defined quantities for apparent place transforms .. */
#define XFORM_TOSTDEP  -1
#define XFORM_JUSTPRE  1
#define XFORM_DOAPPAR  0
#define XFORM_DOABER   1
#define XFORM_NOABER   0

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
#define  EARTHRAD_IN_AU    23454.7910556298  /* number of earth rad in 1 au */
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
int update_on;  
double update_delta;   
#endif

/*
char zname[40];
char znabr[3];
char sname[40];
*/

struct coord
   {
     int sign;  /* carry sign explicitly since -0 not neg. */
     double hh;
     double mm;
     double ss;
   };

struct date_time
   {
	int y;
	int mo;
	int d;
	int h;
	int mn;
	float s;
   };

FILE *sclogfl;

double star_tzero, star_terr, 
	star_period, star_perr;  /* for ephemeris calculations ... global */

void oprntf(char *fmt, ...);

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

/* elements of K&R hp calculator, basis of commands */

char buf[BUFSIZE];
int bufp;


char getch();
 /* get a (possibly pushed back) character */

void ungetch(int c); /* push character back on input */

/* some functions for getting well-tested input. */

int legal_num_part(char);
int legal_int_part(char);
int legal_command_char(char);
int parsedouble(char *,double *);
	/* return values 0 = ok, with number, 1 = found a valid command,
	   but no number, and -1 = an error of some sort (unexpected char)*/
int getdouble(double *, double, double, char *);
int parseint(char *, int *);
int getint(int *, int, int, char *);
double bab_to_dec(struct coord);
void dec_to_bab (double, struct coord *);
int get_line(char *);
double get_coord();
double roundx(double, int);
void round_coord(struct coord *, struct coord *, int);
void put_hrs(double, int, int, int, int);
void put_coords(double, int, int);
void put_colon_coords(double, int, int);
void fput_hrs(FILE *, double, int, int, int, int);
void fput_coords(FILE *, double, int, int);
void load_site(char *, double *, double *, double *, int *,
	char *, char *, double *, double *, double *, char *);
               
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

double atan_circ(double, double);   /* x, y ... use atan2 instead ... */
void min_max_alt(double, double, double *, double *);
			/* latitude, dec, min, max */
double altit(double, double, double, double *, double *);
	/* dec, ha, lat, azimuth, parallactic ... altitude returned. */
double secant_z(double);
double true_airmass(double);     /* arg is sec-z */
double ha_alt(double, double, double);  /* declination, latitude, altitude.*/
double subtend(double, double, double, double);
	/* ra1, dec1, ra2, dec2 */
int get_pm(double, double *, double *);  /* dec, mu-ra (s/yr), mu-dec (") */ 
int get_date(struct date_time *);
int get_time(struct date_time *date); 
double date_to_jd(struct date_time);
int day_of_week(double);  /* takes jd, returns 0=Mon --> 6 = Sun */
#define IGREG 2299161
void caldat(double, struct date_time *, int *);
	/* jd, date_time, day of week */
double day_of_year(double);
void print_day(int);
void print_all(double);
void print_current(struct date_time, int, int);
	/* date, night_date, enter_ut */
void print_calendar(double, int *);
void print_time(double, int);  /* jd, precision. */
double frac_part(double);
double lst(double, double); /* jd, W longitude in hours. */
double adj_time(double);  /* adjusts decimal time to +-12 */
void lpmoon(double, double, double, double *, double *, double *);
	/* double jd,lat,sid,*ra,*dec,*dist;  */
void lpsun(double, double *, double *);
	/* double jd, *ra, *dec; */
void eclrot(double, double *, double *, double *);
	/* double jd, *x, *y, *z; */
double circulo(double x);
	   /* modulo 360 degrees. */
void geocent(double, double, double, double *, double *, double *);
	/* double geolong (dec. hrs), geolat (deg), height (m), 
		*x_geo, *y_geo, *z_geo; */
double etcorr(double);
	/* et correction for a given jd */
void topocorr(double,double,double,double,double,
        double,double *,double *,double *);
/*  topocentric correction, arguments are 
double ra, dec, dist, geolat, lst,
         elev,double * topora, *topodec, *topodist  */

void accumoon(double, double, double, double, double *, double *,double *,
     double *, double *, double *);
  /*            
	double jd,geolat,lst,elevsea;
     	double *geora,*geodec,*geodist,*topora,*topodec,*topodist; */
void flmoon(int, int, double *); 
   /* lunation, phase 0=new, 3=last, jd of phase */        

float lun_age(double, int *);
   /* for jd in, finds day since last new and lunation of last new moon */
void print_phase(double);
   /* verbal description of moon phase given jd */
double lunskybright (double, double, double, double, double, double); 
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

void accusun(double, double, double, double *, double *, double *,
		double *, double *, double *, double *, double *);
	/* (jd,lst,geolat,ra,dec,dist,topora,topodec,x,y,z)  */
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
double jd_moon_alt(double, double, double, double, double);
    /* alt,jdguess,lat,longit,elevsea 
	 returns jd at which moon is at a given 
	altitude, given jdguess as a starting point. In current version
	uses high-precision moon -- execution time does not seem to be
	excessive on modern hardware.  If it's a problem on your machine,
	you can replace calls to 'accumoon' with 'lpmoon' and remove
	the 'elevsea' argument. */

double jd_sun_alt(double, double, double, double);
	/* alt,jdguess,lat,longit; 
	 returns jd at which sun is at a given 
	altitude, given jdguess as a starting point. Uses
	low-precision sun, which is plenty good enough. */

float ztwilight(double);  
/*	double alt;
	evaluates a polynomial expansion for the approximate brightening
   in magnitudes of the zenith in twilight compared to its 
   value at full night, as function of altitude of the sun (in degrees).
   To get this expression I looked in Meinel, A.,
   & Meinel, M., "Sunsets, Twilight, & Evening Skies", Cambridge U.
   Press, 1983; there's a graph on p. 38 showing the decline of 
   zenith twilight.  I read points off this graph and fit them with a
   polynomial; I don't even know what band there data are for! */
/* Comparison with Ashburn, E. V. 1952, JGR, v.57, p.85 shows that this
   is a good fit to his B-band measurements.  */

void find_dst_bounds(int, double, int, double *, double *);
/*	int yr;
	double stdz;
	int use_dst;
  	double *jdb,*jde; 

	 finds jd's at which daylight savings time begins 
	    and ends.  The parameter use_dst allows for a number
	    of conventions, namely:
		0 = don't use it at all (standard time all the time)
		1 = use USA convention (1st Sun in April to
		     last Sun in Oct after 1986; last Sun in April before)
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

double zone(int, double, double, double, double); 
    /*       
	int use_dst;
	double stdz,jd,jdb,jde; 
           
	 Returns zone time offset when standard time zone is stdz,
	   when daylight time begins (for the year) on jdb, and ends
	   (for the year) on jde.  This is parochial to the northern
	   hemisphere. 
	 Extension -- specifying a negative value of use_dst reverses
	   the logic for the Southern hemisphere; then DST is assumed for
	   the Southern hemisphere summer (which is the end and beginning
	   of the year. */

double true_jd(struct date_time, int, int, int, double);
     /* date, use_dst, enter_ut, night_date, stdz 
       
 takes the values in the date-time structure, the standard time
   zone (in hours west), the prevailing conventions for date and
   time entry, and returns the value of the true julian date. */

void print_tz(double, int, double, double , char);
             
     /*	double jd;
	int use;
	double jdb,jde;
	char zabr; 
	 prints correct time abbreviation, given zabr as the
	   single character abbreviation for the time zone,
	   "D" or "S" depending on daylight or standard (dst 
	    begins at jdb, ends at jde) and current jd. */

void xyz_cel(double, double, double, double *, double *);
     /* x,y,z,ra,dec 
	double x,y,z;   cartesian coordinate triplet 
	double *ra, *dec;   corresponding right ascension and declination,
                returned in decimal hours and decimal degrees. */

void aberrate(double, double *, int);
    /*  epoch, vec, from_std) 
	double epoch,   decimal year ...  
	vec[];   celestial unit vector ...  
        int from_std;   1 = apply aberration, -1 = take aberration out. */

void nutation_params(double, double *, double *);
	/* double date_epoch, *del_psi, *del_ep;

   computes the nutation parameters delta psi and
   delta epsilon at julian epoch (in years) using approximate
   formulae given by Jean Meeus, Astronomical Formulae for
   Calculators, Willman-Bell, 1985, pp. 69-70. Accuracy
   appears to be a few hundredths of an arcsec or better
   and numerics have been checked against his example. 
   Nutation parameters are returned in radians. */

void cooxform(double, double, double, double, double *, double *,
	     int, int, int);
/*	 cooxform(rin, din, std_epoch,
  date_epoch, rout, dout, just_precess, do_aber, from_std)
     	double rin, din;   input ra and dec 
	double std_epoch;
	double date_epoch;        
	double *rout, *dout;   output 
	int just_precess;   flag ... 1 does just precession, 0 
			includes aberration and nutation. 
        int do_aber;    flag ... 1 does aberration, 0 does not. 
	int from_std;     flag ... 1 --> from std to date,
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

double near_hor_refr(double, double);
	/* double app_alt, pressure  */
	/* Almanac 1992, p. B62 -- ignores temperature variation */
	/* formula for near horizon, function-ized for iteration ... */

double refract_size(double, double);
	/* double alt    altitude in degrees 
	double elev;   meters */
 	/* Almanac for 1992, p. B 62.  Ignores variation in temperature
           and just assumes T = 20 celsius.  */
void refract_corr(double *, double *, double, double, double *, int);
	/* ha , dec, lat, elev, size, sense
	double *ha, *dec, lat, eleve, *size;
	int sense; */
/* if sense == 1 , applies refraction to a true ha and dec; if
   == -1, de-corrects already refracted coordinates. Uses elevation of
   observatory above sea level to estimate a mean atmospheric pressure. */

void mass_precess();

void print_apparent(double,double,double,double,double,double,double,double,
		double);
/* (rain,decin,epochin,mura_sec,mudec,jd,lat,longit,elev) 
	double rain, decin, epochin, mura_sec, mudec, jd, lat, longit, elev;
*/

void galact(double, double, double, double *, double *);
           
	/* double ra,dec,epoch,*glong,*glat;
	 Homebrew algorithm for 3-d Euler rotation into galactic. 
	   Perfectly rigorous, and with reasonably accurate input 
	   numbers derived from original IAU definition of galactic
	   pole (12 49, +27.4, 1950) and zero of long (at PA 123 deg
	   from pole.) */

void gal2radec(double, double, double, double *, double *);

	/* inverts the above. */

void eclipt(double, double, double, double, double *, double *, double *);
	/* ra,dec,epoch,jd,curep,eclong,eclat */
	 /* ra in decimal hrs, other coords in dec. deg. */
	/* converts ra and dec to ecliptic coords -- precesses to current
   epoch first (and hands current epoch back for printing.)  */

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

void comp_el(double);
	/* double jd; */
	
void planetxyz(int, double, double *, double *, double *);
             
	/*	int p; 
	double jd, *x, *y, *z;
    produces ecliptic x,y,z coordinates for planet number 'p'
   at date jd. */

void planetvel(int, double, double *, double *, double *);
	/* ( p, jd, vx, vy, vz)
	int p; 
	double jd, *vx, *vy, *vz;
	numerically evaluates planet velocity by brute-force
	numerical differentiation. Very unsophisticated algorithm. */
	/* answer should be in ecliptic coordinates, in AU per day.*/

void xyz2000(double jd, double x, double y, double z); 
	/*	double jd, x, y, z;
simply transforms a vector x, y, and z to 2000 coordinates
   and prints -- for use in diagnostics. */

void earthview(double *x, double *y, double *, int, double *, double *);
	/* double *x, *y, *z;
	int i; 
	double *ra, *dec;
    given computed planet positions for planets 1-10, computes
   ra and dec of i-th planet as viewed from earth (3rd) */

void pposns(double,double, double, int, double *, double *);
     /* (jd,lat,sid,print_option,planra,plandec) 
	double jd,lat,sid;
	int print_option;
	double *planra, *plandec; 
                      
    computes and optionally prints positions for all the planets.
     print_option 1 = print positions, 0 = silent */

void barycor(double, double *, double *, double *, double *, double *, double *);
	/*	jd,x,y,z,xdot,ydot,zdot)

	double jd,*x,*y,*z;
	double *xdot,*ydot,*zdot;
            
   This routine takes the position
   x,y,z and velocity xdot,ydot,zdot, assumed heliocentric,
   and corrects them to the solar system barycenter taking into
   account the nine major planets.  Routine evolved by inserting
   planetary data (given above) into an earlier, very crude
   barycentric correction.  */

void helcor(double, double, double, double, double, double, double *, double*); 

	/* jd,ra,dec,ha,lat,elevsea,tcor,vcor
	double jd,ra,dec,ha;
  	double lat,elevsea,*tcor,*vcor;
           
   finds heliocentric correction for given jd, ra, dec, ha, and lat.
   tcor is time correction in seconds, vcor velocity in km/s, to 
   be added to the observed values.  
   Input ra and dec assumed to be at current epoch */

/* A couple of eclipse predictors.... */

float overlap(double, double, double);
	/* double r1, r2, sepn; 
   for two circles of radii r1 and r2,
   computes the overlap area in terms of the area of r1
   if their centers are separated by
   sepn. */

void solecl(double, double, double);
    /* sun_moon,distmoon,distsun) 
	double sun_moon,distmoon,distsun;   */

int lunecl(double, double, double, double, double, double);

    /* (georamoon,geodecmoon,geodistmoon,rasun,decsun,distsun)
	 quickie lunar eclipse predictor -- makes a number of 
	   minor assumptions, e. g. small angle approximations, plus
	   projects phenomena onto a plane at distance = geocentric 
	   distance of moon . */
 
void planet_alert(double, double, double, double);
                 
	/*	double jd,ra,dec,tolerance;
                 
/* given a jd, ra, and dec, this computes rough positions
   for all the planets, and alerts the user if any of them
   are within less than a settable tolerance of the ra and dec. */


int setup_time_place(struct date_time date, double, double, double, 
	int, char *, char, char *, int, int, double *, double *, double *,
	double *, double *, double *);

  /*  setup_time_place(date,longit,lat,stdz,use_dst,zone_name,
        zabr,site_name,enter_ut,night_date,jdut,jdlocal,jdb,jde,sid,
	curepoch)
  
struct date_time date;
double lat, longit, stdz, *jdut, *jdlocal, *jdb, *jde, *sid, *curepoch;
int use_dst, enter_ut, night_date;
char zabr;
char *site_name;
char *zone_name;
       
 This takes the date (which contains the time), and the site parameters,
   and prints out a banner giving the various dates and times; also
   computes and returns various jd's, the sidereal time, and the epoch. 
   Returns negative number to signal error if date is out of range of
   validity of algorithms, or if you specify a bad time during daylight-time
   change; returns zero if successful.  */


void print_tonight(struct date_time,double,double,double,double,double,
 	char *, double, char *, char, int, double *, double *, int);
/*  print_tonight(date,lat,longit,elevsea,elev,horiz,site_name,stdz,
	zone_name,zabr,use_dst,jdb,jde,int_long)

struct date_time date;
double lat, longit, elevsea, elev, horiz, stdz, *jdb, *jde;
char *site_name, *zone_name, zabr;
int use_dst, int_long;      int_long is a fossil argument which 
                    allows a slightly inter version to be printed. 
 Given site and time information, prints a summary of
   the important phenomena for a single night. 
   The coding in this routine is extremely tortuous ... I even use
   the dreaded goto statement!  It's inelegant, but (a) the logic is
   indeed somewhat complicated and (b) it works.  */

void print_circumstances(double, double, double, double, double,
        double, double, double, double, double, double, double);
	/* (objra,objdec,objepoch,jd,curep,
	mura_arcs,mura_sec,mudec,sid,lat,elevsea,horiz) 

	double objra,objdec,objepoch,curep,mura_arcs,mura_sec,mudec,lat,horiz;
	double jd,sid,elevsea;

 Given object, site, and time information, prints the circumstances
   of an observation.  The heart of the "calculator" mode. */

void hourly_airmass(struct date_time, double, double, double, double, int,
	double, double, double, double, double, double, char *);
/*  hourly_airmass(date,stdz,lat,longit,horiz,use_dst,objra,objdec,
  objepoch, mura_sec,mura_arcs,mudec)

   Given a slew of information, prints a table of hourly airmass, etc.
   for use in scheduling observations.  Also prints sun and moon
   altitude when these are relevant.  Precesses coordinates as well. */


void print_params(struct date_time,int, int, double, double,
    double, char *, double, double, int, double, double, double, 
    double, double, double);
 /*  print_params(date,enter_ut,night_date,stdz,lat,longit,site_name,
    elevsea,elev,use_dst,objra,objdec,objepoch,mura_sec,mura_arcs,
    mudec)

    struct date_time date;
    int enter_ut;
    int night_date;
    double stdz;
    double lat;
    double longit;
    char *site_name;
    double elevsea;
    double elev;
    int use_dst;
    double objra;
    double objdec;
    double objepoch;
    double mura_sec;
    double mura_arcs;
    double mudec;
		 
 This simply prints a nicely formatted list of the *input* parameters
   without doing any computations.  Helpful for the perplexed user, and
   for checking things. */

void print_menu(); 
void print_tutorial();
void print_examples();
void print_accuracy(); 
void print_legalities(); 
void ephemgen(double, double, double, double, double);
  /*  double ra, dec, ep, lat, longit

 Prompts for elements of an ephemeris, then 
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

#define ALT_3  19.47  /* 19.47 degrees altitude => sec z = 3 */
#define ALT_2  30.
#define ALT_15 41.81
#define SID_RATE 1.0027379093  /* sidereal / solar rate */

double hrs_up(double, double, double, double);

	/* 	double jdup, jddown, jdeve, jdmorn;
    If an object comes up past a given point at jdup,
      and goes down at jddown, and evening and morning
      are at jdeve and jdmorn, computes how long
      object is up *and* it's dark.  ... Written as
      function 'cause it's done three times below. */

void print_air(double, int);
    	/* double secz;
	int prec; */
void print_ha_air(double, double, int, int);
	/* double ha, secz;
	int prec1, prec2; */
void obs_season(double, double, double, double, double);
  /*  double ra, dec, epoch, lat, longit;	
   prints a table of observability through an observing
   season.  The idea is to help the observer come up
   with an accurately computed "range of acceptable
   dates", to quote NOAO proposal forms ... */

#if SYS_CLOCK_OK == 1
#include <time.h>
#endif

#if SYS_CLOCK_OK == 1

int get_sys_date(struct date_time *, int, int, int, double, double);
	/* date, use_dst, enter_ut, night_date, stdz, toffset) 
	Reads the system clock; loads up the date structure
           to conform to the prevailing conventions for the interpretation
           of times.  Optionally adds "toffset" minutes to the system
           clock, as in x minutes in the future. */
#endif

void indexx(int, float *, int *);

/* Sort routine from Press et al., "Numerical Recipes in C",
  1st edition, Cambridge University Press. 

	int n,indx[];
	float arrin[];
*/

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

int read_obj_list(); 

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

int find_by_name(double *, double *, double, struct date_time, int, int, 
     int, double, double, double);
  /*  find_by_name(ra, dec, epoch, date, use_dst, enter_ut, night_date, stdz,
		lat, longit)      	
	double *ra, *dec, epoch, stdz, lat, longit;
	struct date_time date;	
	int use_dst, enter_ut, night_date;
	finds object by name in list, and sets ra and dec to
           those coords if found.  Precesses to current value of
           epoch. */

void type_list(struct date_time date, int, int, int, 
      double, double, double);
	/* type_list(date, use_dst, enter_ut, night_date, stdz,
                lat, longit)
        double stdz, lat, longit;
        struct date_time date;
        int use_dst, enter_ut, night_date;
	*/

int find_nearest(double *, double *, double, struct date_time,
      int, int, int, double, double, double);
/* find_nearest(ra, dec, epoch, date, use_dst, enter_ut, night_date, stdz,
		lat, longit)     
    
	double *ra, *dec, epoch, stdz, lat, longit;
	struct date_time date;	
	int use_dst, enter_ut, night_date;

   given ra,dec, and epoch, sorts items in list with 
   respect to arc distance away, and queries user
   whether to accept.  */

void set_zenith(struct date_time, int, int, int, double, double, 
	double, double, double *, double *);
/*  set_zenith(date, use_dst, enter_ut, night_date, stdz, lat, 
	  longit, epoch, ra, dec) 
               
 sets RA and dec to zenith as defined by present time and date;
   coords are set to actual zenith but in currently used epoch.  */

void printephase(struct date_time, int, int, int, 
     double, double, double, double, double, double);
/*  printephase(date, use_dst, enter_ut, night_date, stdz, lat, 
	  longit, epoch, ra, dec) 
   prints phase of a repeating phenomenon at this instant. */

int set_to_jd(struct date_time *, int, int, int, 
     double, double);
	/*  set_to_jd(date, use_dst, enter_ut, night_date, stdz, jd) 
	 Takes a jd and loads up the date structure
           to conform to the prevailing conventions for the interpretation
           of times. */

void radec_to_constel(double, double, double, char *);  
        /* void radec_to_constel(double ra, double dec, double epoch, 
		char *constname)  */
