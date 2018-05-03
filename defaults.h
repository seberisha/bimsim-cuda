#ifndef BIMSIM_DEFAULTS
#define BIMSIM_DEFAULTS

//default sphere parameters
#define DEFAULT_SPHERE_X		0
#define DEFAULT_SPHERE_Y		0
#define DEFAULT_SPHERE_Z		0
#define DEFAULT_SPHERE_A		1

//default near field parameters
#define DEFAULT_LAMBDA          "1"
#define DEFAULT_AMPLITUDE		"1"
#define DEFAULT_MATERIAL        "1.4 0.05"
#define DEFAULT_N               1.4
#define DEFAULT_K               0.5
#define DEFAULT_FOCUS           "0 0 0"
//#define DEFAULT_FOCUS_X         "0"
//#define DEFAULT_FOCUS_Y         "0"
//#define DEFAULT_FOCUS_Z         "0"
//#define DEFAULT_INCIDENT_ORDER	20
#define DEFAULT_STABILITY_PARM	1.4

//optics
//#define DEFAULT_CONDENSER_MIN   "0.0"
//#define DEFAULT_CONDENSER_MAX   "1.0"
#define DEFAULT_CONDENSER       "0 1"
#define DEFAULT_POINTS_RINGS    "1 1"

//#define DEFAULT_OBJECTIVE_MIN   "0"
//#define DEFAULT_OBJECTIVE_MAX   "1"
#define DEFAULT_OBJECTIVE       "0 1"

//incident light direction
#define DEFAULT_K_THETA			0
#define DEFAULT_K_PHI			0

//default flags
#define DEFAULT_PLANEWAVE		false
#define DEFAULT_VECTORSIM		false
#define DEFAULT_APPEND          false
//#define DEFAULT_OUTPUT_POINT	fileoutStruct::imageObjective

#define DEFAULT_PLANE_MIN       "-5 0 -5"
#define DEFAULT_PLANE_MIN_X     -5
#define DEFAULT_PLANE_MIN_Y     0
#define DEFAULT_PLANE_MIN_Z     -5

#define DEFAULT_PLANE_MAX       "5 0 5"
#define DEFAULT_PLANE_MAX_X     5
#define DEFAULT_PLANE_MAX_Y     0
#define DEFAULT_PLANE_MAX_Z     5

#define DEFAULT_PLANE_NORM      "0 1 0"
#define DEFAULT_PLANE_NORM_X    0
#define DEFAULT_PLANE_NORM_Y    1
#define DEFAULT_PLANE_NORM_Z    0

#define DEFAULT_PLANE_SIZE		40
#define DEFAULT_PLANE_POSITION	0


/*
#define DEFAULT_SLICE_MIN_X     -6
#define DEFAULT_SLICE_MIN_Y     -6
#define DEFAULT_SLICE_MIN_Z     1

#define DEFAULT_SLICE_MAX_X     6
#define DEFAULT_SLICE_MAX_Y     6
#define DEFAULT_SLICE_MAX_Z     1

#define DEFAULT_SLICE_NORM_X    0
#define DEFAULT_SLICE_NORM_Y    0
#define DEFAULT_SLICE_NORM_Z    1
*/


#define DEFAULT_FIELD_ORDER     "10"

#define DEFAULT_SAMPLES         "400"

#define DEFAULT_SLICE_RES		"256"

#define DEFAULT_SPHERE_THETA_R  1000

#define DEFAULT_PADDING			"1"
#define DEFAULT_SUPERSAMPLE		"1"

#define DEFAULT_INTENSITY_FILE	    "out_i.bmp"
#define DEFAULT_INCIDENT_FILE	    "out_inc.bmp"

#define DEFAULT_TRANSMITTANCE_FILE	""
#define DEFAULT_ABSORBANCE_FILE	    "out_a.bmp"

#define DEFAULT_ABSORBANCE_SPECTRUM_FILE	    "out_as.bmp"

#define DEFAULT_NEAR_FILE		    "out_n.bmp"
#define DEFAULT_FAR_FILE		    "out_f.bmp"
#define DEFAULT_EXTENDED_SOURCE     ""

#define DEFAULT_INTERPOLATE     "-1"

#define DEFAULT_FIELD_TYPE		    "magnitude"
#define DEFAULT_FORMAT			    fileoutStruct::formatImage
#define DEFAULT_COLORMAP		    "brewer"


#endif
